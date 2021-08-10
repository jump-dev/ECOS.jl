using MathOptInterface
const MOI = MathOptInterface
const CI = MOI.ConstraintIndex
const VI = MOI.VariableIndex

const MOIU = MOI.Utilities

struct Solution
    ret_val::Union{Nothing,Int}
    primal::Vector{Float64}
    dual_eq::Vector{Float64}
    dual_ineq::Vector{Float64}
    slack::Vector{Float64}
    objective_value::Float64
    dual_objective_value::Float64
    objective_constant::Float64
    solve_time::Float64
    iter::Int
end
# The values used by ECOS are from -7 to 10 so -10 should be safe
function Solution()
    return Solution(
        nothing,
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        NaN,
        NaN,
        NaN,
        NaN,
        0,
    )
end

MOIU.@product_of_sets(Zeros, MOI.Zeros)
MOIU.@product_of_sets(PointedCones, MOI.Nonnegatives, MOI.SecondOrderCone, MOI.ExponentialCone)

MOIU.@struct_of_constraints_by_set_types(
    ZerosOrNot,
    MOI.Zeros,
    Union{MOI.Nonnegatives,MOI.SecondOrderCone,MOI.ExponentialCone},
)

const OptimizerCache = MOI.Utilities.GenericModel{
    Cdouble,
    ZerosOrNot{Cdouble}{
        MOIU.MatrixOfConstraints{
            Cdouble,
            MOIU.MutableSparseMatrixCSC{
                Cdouble,
                Clong,
                MOI.Utilities.ZeroBasedIndexing,
            },
            Vector{Cdouble},
            Zeros{Cdouble},
        },
        MOIU.MatrixOfConstraints{
            Cdouble,
            MOIU.MutableSparseMatrixCSC{
                Cdouble,
                Clong,
                MOI.Utilities.ZeroBasedIndexing,
            },
            Vector{Cdouble},
            PointedCones{Cdouble},
        },
    },
}

mutable struct Optimizer <: MOI.AbstractOptimizer
    inner::Union{Nothing,Ptr{Cpwork}}
    maxsense::Bool
    sol::Solution
    silent::Bool
    options::Dict{String,Any}
    function Optimizer(; kwargs...)
        optimizer = new(
            nothing,
            false,
            Solution(),
            false,
            Dict{String,Any}(),
        )
        if length(kwargs) > 0
            @warn(
                "Passing keyword attributes is deprecated. Use " *
                "`set_optimizer_attribute` instead.",
            )
        end
        for (key, value) in kwargs
            MOI.set(optimizer, MOI.RawOptimizerAttribute(String(key)), value)
        end
        return optimizer
    end
end

MOI.get(::Optimizer, ::MOI.SolverName) = "ECOS"

function MOI.supports(::Optimizer, param::MOI.RawOptimizerAttribute)
    return hasfield(Csettings, Symbol(param.name))
end

function MOI.set(optimizer::Optimizer, param::MOI.RawOptimizerAttribute, value)
    optimizer.options[param.name] = value
    return
end

function MOI.get(optimizer::Optimizer, param::MOI.RawOptimizerAttribute)
    # TODO: This gives a poor error message if the name of the parameter is
    # invalid.
    return optimizer.options[param.name]
end

MOI.supports(::Optimizer, ::MOI.Silent) = true
function MOI.set(optimizer::Optimizer, ::MOI.Silent, value::Bool)
    return optimizer.silent = value
end
MOI.get(optimizer::Optimizer, ::MOI.Silent) = optimizer.silent

function MOI.is_empty(optimizer::Optimizer)
    return !optimizer.maxsense && optimizer.inner === nothing
end

function MOI.empty!(optimizer::Optimizer)
    optimizer.maxsense = false
    optimizer.inner = nothing
    optimizer.sol = Solution()
    return
end

function MOI.supports(
    ::Optimizer,
    ::Union{
        MOI.ObjectiveSense,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
    },
)
    return true
end

const AFF = MOI.VectorAffineFunction{Float64}

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{AFF},
    ::Type{
        <:Union{
            MOI.Zeros,
            MOI.Nonnegatives,
            MOI.SecondOrderCone,
            MOI.ExponentialCone,
        },
    },
)
    return true
end

# ECOS orders differently than MOI the second and third dimension of the exponential cone
orderval(val, s) = val
function orderval(val, s::Union{MOI.ExponentialCone,Type{MOI.ExponentialCone}})
    return val[[1, 3, 2]]
end
orderidx(idx, s) = idx
expmap(i) = (1, 3, 2)[i]
function orderidx(idx, s::MOI.ExponentialCone)
    return expmap.(idx)
end

function _copy_to(dest::Optimizer, src::OptimizerCache)
    @assert MOI.is_empty(dest)
    Ab = MOI.Utilities.constraints(src.constraints, AFF, MOI.Zeros)
    A = Ab.coefficients
    Gh = MOI.Utilities.constraints(src.constraints, AFF, MOI.Nonnegatives)
    G = Gh.coefficients
    q = MOI.dimension.(MOI.get.(
        Gh,
        MOI.ConstraintSet(),
        MOI.get(Gh, MOI.ListOfConstraintIndices{AFF,MOI.SecondOrderCone}())
    ))
    @assert A.n == G.n
    obj =
        MOI.get(src, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    c0 = zeros(A.n)
    for term in obj.terms
        c0[term.variable.value] += term.coefficient
    end
    dest.inner = _setup(
        A.n,
        G.m,
        A.m,
        Gh.sets.num_rows[2] - Gh.sets.num_rows[1],
        length(q),
        q,
        MOI.get(Gh, MOI.NumberOfConstraints{AFF,MOI.ExponentialCone}()),
        ECOSMatrix(-G.nzval, G.colptr, G.rowval),
        ECOSMatrix(-A.nzval, A.colptr, A.rowval),
        dest.maxsense ? -c0 : c0,
        Gh.constants,
        Ab.constants;
    )
    dest.maxsense = MOI.get(src, MOI.ObjectiveSense()) == MOI.MAX_SENSE
    return
end

function MOI.copy_to(
    dest::Optimizer,
    src::OptimizerCache;
    copy_names::Bool = false,
)
    _copy_to(dest, src)
    return MOIU.identity_index_map(src)
end

function MOI.copy_to(
    dest::Optimizer,
    src::MOI.Utilities.UniversalFallback{OptimizerCache};
    copy_names::Bool = false,
)
    return MOI.copy_to(dest, src.model)
end

function MOI.copy_to(
    dest::Optimizer,
    src::MOI.ModelLike;
    copy_names::Bool = false,
)
    cache = OptimizerCache()
    index_map = MOI.copy_to(cache, src)
    _copy_to(dest, cache)
    return index_map
end

function MOI.optimize!(optimizer::Optimizer)
    if optimizer.data === nothing
        # optimize! has already been called and no new model has been copied
        return
    end
    options = Dict(Symbol(k) => v for (k, v) in optimizer.options)
    if optimizer.silent
        options[:verbose] = false
    end
    settings(optimizer.inner, dest.options)
    ret_val = ECOS.solve(optimizer.inner)
    ecos_prob = unsafe_load(optimizer.inner)
    stat = unsafe_load(ecos_prob.info)
    solve_time = stat.tsetup + stat.tsolve
    iter = stat.iter
    primal = unsafe_wrap(Array, ecos_prob.x, n)[:]
    dual_eq = unsafe_wrap(Array, ecos_prob.y, cone.f)[:]
    dual_ineq = unsafe_wrap(Array, ecos_prob.z, m)[:]
    slack = unsafe_wrap(Array, ecos_prob.s, m)[:]
    ECOS.cleanup(optimizer.inner, 0)
    objective_value = (optimizer.maxsense ? -1 : 1) * stat.pcost
    dual_objective_value = (optimizer.maxsense ? -1 : 1) * stat.dcost
    optimizer.sol = Solution(
        ret_val,
        primal,
        dual_eq,
        dual_ineq,
        slack,
        objective_value,
        dual_objective_value,
        objective_constant,
        solve_time,
        iter,
    )
    return
end

MOI.get(optimizer::Optimizer, ::MOI.SolveTimeSec) = optimizer.sol.solve_time
MOI.get(optimizer::Optimizer, ::MOI.BarrierIterations) = optimizer.sol.iter
function MOI.get(optimizer::Optimizer, ::MOI.RawStatusString)
    # Strings from https://github.com/ifa-ethz/ecos/blob/master/include/ecos.h
    flag = optimizer.sol.ret_val
    if flag === nothing
        return "Optimize not called"
    elseif flag == ECOS_OPTIMAL
        return "Problem solved to optimality"
    elseif flag == ECOS_OPTIMAL + ECOS_INACC_OFFSET
        return "Problem solved to inaccurate optimality"
    elseif flag == ECOS_PINF
        return "Found certificate of primal infeasibility"
    elseif flag == ECOS_PINF + ECOS_INACC_OFFSET
        return "Found inaccurate certificate of primal infeasibility"
    elseif flag == ECOS_DINF
        return "Found certificate of dual infeasibility"
    elseif flag == ECOS_DINF + ECOS_INACC_OFFSET
        return "Found inaccurate certificate of dual infeasibility"
    elseif flag == ECOS_MAXIT
        return "Maximum number of iterations reached"
    elseif flag == ECOS_NUMERICS
        return "Search direction unreliable"
    elseif flag == ECOS_OUTCONE
        return "s or z got outside the cone, numerics?"
    elseif flag == ECOS_SIGINT
        return "solver interrupted by a signal/ctrl-c"
    else
        @assert flag == ECOS_FATAL
        return "Unknown problem in solver"
    end
end

# Implements getter for result value and statuses
function MOI.get(optimizer::Optimizer, ::MOI.TerminationStatus)
    flag = optimizer.sol.ret_val
    if flag === nothing
        return MOI.OPTIMIZE_NOT_CALLED
    elseif flag == ECOS.ECOS_OPTIMAL
        return MOI.OPTIMAL
    elseif flag == ECOS.ECOS_PINF
        return MOI.INFEASIBLE
    elseif flag == ECOS.ECOS_DINF
        return MOI.DUAL_INFEASIBLE
    elseif flag == ECOS.ECOS_MAXIT
        return MOI.ITERATION_LIMIT
    elseif flag == ECOS.ECOS_NUMERICS || flag == ECOS.ECOS_OUTCONE
        return MOI.NUMERICAL_ERROR
    elseif flag == ECOS.ECOS_SIGINT
        return MOI.INTERRUPTED
    elseif flag == ECOS.ECOS_FATAL
        return MOI.OTHER_ERROR
    elseif flag == ECOS.ECOS_OPTIMAL + ECOS.ECOS_INACC_OFFSET
        return MOI.ALMOST_OPTIMAL
    elseif flag == ECOS.ECOS_PINF + ECOS.ECOS_INACC_OFFSET
        return MOI.ALMOST_INFEASIBLE
    elseif flag == ECOS.ECOS_DINF + ECOS.ECOS_INACC_OFFSET
        return MOI.ALMOST_DUAL_INFEASIBLE
    else
        error("Unrecognized ECOS solve status flag: $flag.")
    end
end

function MOI.get(optimizer::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(optimizer, attr)
    value = optimizer.sol.objective_value
    if !MOIU.is_ray(MOI.get(optimizer, MOI.PrimalStatus()))
        value += optimizer.sol.objective_constant
    end
    return value
end
function MOI.get(optimizer::Optimizer, attr::MOI.DualObjectiveValue)
    MOI.check_result_index_bounds(optimizer, attr)
    value = optimizer.sol.dual_objective_value
    if !MOIU.is_ray(MOI.get(optimizer, MOI.DualStatus()))
        value += optimizer.sol.objective_constant
    end
    return value
end

function MOI.get(optimizer::Optimizer, attr::MOI.PrimalStatus)
    if attr.result_index > MOI.get(optimizer, MOI.ResultCount())
        return MOI.NO_SOLUTION
    end
    flag = optimizer.sol.ret_val
    if flag == ECOS.ECOS_OPTIMAL
        return MOI.FEASIBLE_POINT
    elseif flag == ECOS.ECOS_PINF
        return MOI.INFEASIBLE_POINT
    elseif flag == ECOS.ECOS_DINF
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif flag == ECOS.ECOS_MAXIT
        return MOI.UNKNOWN_RESULT_STATUS
    elseif flag == ECOS.ECOS_OPTIMAL + ECOS.ECOS_INACC_OFFSET
        return MOI.NEARLY_FEASIBLE_POINT
    elseif flag == ECOS.ECOS_PINF + ECOS.ECOS_INACC_OFFSET
        return MOI.INFEASIBLE_POINT
    elseif flag == ECOS.ECOS_DINF + ECOS.ECOS_INACC_OFFSET
        return MOI.NEARLY_INFEASIBILITY_CERTIFICATE
    else
        return MOI.OTHER_RESULT_STATUS
    end
end
# Swapping indices 2 <-> 3 is an involution (it is its own inverse)
const reorderval = orderval
function MOI.get(optimizer::Optimizer, attr::MOI.VariablePrimal, vi::VI)
    MOI.check_result_index_bounds(optimizer, attr)
    return optimizer.sol.primal[vi.value]
end
function MOI.get(optimizer::Optimizer, a::MOI.VariablePrimal, vi::Vector{VI})
    return MOI.get.(optimizer, Ref(a), vi)
end
function MOI.get(
    optimizer::Optimizer,
    attr::MOI.ConstraintPrimal,
    ci::CI{<:MOI.AbstractFunction,MOI.Zeros},
)
    MOI.check_result_index_bounds(optimizer, attr)
    rows = constrrows(optimizer, ci)
    return zeros(length(rows))
end
function MOI.get(
    optimizer::Optimizer,
    attr::MOI.ConstraintPrimal,
    ci::CI{<:MOI.AbstractFunction,S},
) where {S<:MOI.AbstractSet}
    MOI.check_result_index_bounds(optimizer, attr)
    offset = constroffset(optimizer, ci)
    rows = constrrows(optimizer, ci)
    return reorderval(optimizer.sol.slack[offset.+rows], S)
end

function MOI.get(optimizer::Optimizer, attr::MOI.DualStatus)
    if attr.result_index > MOI.get(optimizer, MOI.ResultCount())
        return MOI.NO_SOLUTION
    end
    flag = optimizer.sol.ret_val
    if flag == ECOS.ECOS_OPTIMAL
        return MOI.FEASIBLE_POINT
    elseif flag == ECOS.ECOS_PINF
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif flag == ECOS.ECOS_DINF
        return MOI.INFEASIBLE_POINT
    elseif flag == ECOS.ECOS_MAXIT
        return MOI.UNKNOWN_RESULT_STATUS
    elseif flag == ECOS.ECOS_OPTIMAL + ECOS.ECOS_INACC_OFFSET
        return MOI.NEARLY_FEASIBLE_POINT
    elseif flag == ECOS.ECOS_PINF + ECOS.ECOS_INACC_OFFSET
        return MOI.NEARLY_INFEASIBILITY_CERTIFICATE
    elseif flag == ECOS.ECOS_DINF + ECOS.ECOS_INACC_OFFSET
        return MOI.INFEASIBLE_POINT
    else
        return MOI.OTHER_RESULT_STATUS
    end
end
function _dual(optimizer, ci::CI{<:MOI.AbstractFunction,<:MOI.Zeros})
    return optimizer.sol.dual_eq
end
_dual(optimizer, ci::CI) = optimizer.sol.dual_ineq
function MOI.get(
    optimizer::Optimizer,
    attr::MOI.ConstraintDual,
    ci::CI{<:MOI.AbstractFunction,S},
) where {S<:MOI.AbstractSet}
    MOI.check_result_index_bounds(optimizer, attr)
    offset = constroffset(optimizer, ci)
    rows = constrrows(optimizer, ci)
    return reorderval(_dual(optimizer, ci)[offset.+rows], S)
end

MOI.get(optimizer::Optimizer, ::MOI.ResultCount) = 1
