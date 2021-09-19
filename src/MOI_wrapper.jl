const MOIU = MOI.Utilities

const AFF = MOI.VectorAffineFunction{Float64}

mutable struct Solution
    ret_val::Union{Nothing,Int}
    primal::Vector{Float64}
    dual_eq::Vector{Float64}
    dual_ineq::Vector{Float64}
    slack::Vector{Float64}
    objective_value::Float64
    dual_objective_value::Float64
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
        0,
    )
end

MOIU.@product_of_sets(Zeros, MOI.Zeros)
MOIU.@product_of_sets(
    Cones,
    MOI.Nonnegatives,
    MOI.SecondOrderCone,
    PermutedExponentialCone
)

MOIU.@struct_of_constraints_by_set_types(
    ZerosOrNot,
    MOI.Zeros,
    Union{MOI.Nonnegatives,MOI.SecondOrderCone,PermutedExponentialCone},
)

const OptimizerCache = MOI.Utilities.GenericModel{
    Cdouble,
    MOIU.ObjectiveContainer{Cdouble},
    MOIU.VariablesContainer{Cdouble},
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
            Cones{Cdouble},
        },
    },
}

mutable struct Optimizer <: MOI.AbstractOptimizer
    zeros::Union{Nothing,Zeros{Cdouble}}
    cones::Union{Nothing,Cones{Cdouble}}
    sol::Solution
    silent::Bool
    options::Dict{String,Any}
    function Optimizer(; kwargs...)
        optimizer = new(nothing, nothing, Solution(), false, Dict{String,Any}())
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

function MOI.get(::Optimizer, ::MOI.Bridges.ListOfNonstandardBridges)
    return [PermutedExponentialBridge{Cdouble}]
end

function _rows(optimizer::Optimizer, ci::MOI.ConstraintIndex{AFF,MOI.Zeros})
    return MOIU.rows(optimizer.zeros, ci)
end
function _rows(optimizer::Optimizer, ci::MOI.ConstraintIndex{AFF})
    return MOIU.rows(optimizer.cones, ci)
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
    return optimizer.zeros === nothing && optimizer.cones === nothing
end

function MOI.empty!(optimizer::Optimizer)
    optimizer.zeros = nothing
    optimizer.cones = nothing
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

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{AFF},
    ::Type{
        <:Union{
            MOI.Zeros,
            MOI.Nonnegatives,
            MOI.SecondOrderCone,
            PermutedExponentialCone,
        },
    },
)
    return true
end

function _optimize!(dest::Optimizer, src::OptimizerCache)
    MOI.empty!(dest)
    Ab = MOI.Utilities.constraints(src.constraints, AFF, MOI.Zeros)
    A = Ab.coefficients
    Gh = MOI.Utilities.constraints(src.constraints, AFF, MOI.Nonnegatives)
    G = Gh.coefficients
    q =
        MOI.dimension.(
            MOI.get.(
                Gh,
                MOI.ConstraintSet(),
                MOI.get(
                    Gh,
                    MOI.ListOfConstraintIndices{AFF,MOI.SecondOrderCone}(),
                ),
            ),
        )
    @assert A.n == G.n
    max_sense = MOI.get(src, MOI.ObjectiveSense()) == MOI.MAX_SENSE
    obj =
        MOI.get(src, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    objective_constant = MOI.constant(obj)
    c0 = zeros(A.n)
    for term in obj.terms
        c0[term.variable.value] += term.coefficient
    end

    dest.zeros = deepcopy(Ab.sets) # TODO copy(Ab.sets)
    dest.cones = deepcopy(Gh.sets) # TODO copy(Gh.sets)

    options = Dict(Symbol(k) => v for (k, v) in dest.options)
    if dest.silent
        options[:verbose] = false
    end

    inner = _setup(
        A.n,
        G.m,
        A.m,
        Gh.sets.num_rows[1],
        length(q),
        q,
        MOI.get(Gh, MOI.NumberOfConstraints{AFF,PermutedExponentialCone}()),
        ECOSMatrix(-G.nzval, G.colptr, G.rowval),
        ECOSMatrix(-A.nzval, A.colptr, A.rowval),
        max_sense ? -c0 : c0,
        Gh.constants,
        Ab.constants,
    )
    settings(inner, options)
    ret_val = ECOS.solve(inner)
    ecos_prob = unsafe_load(inner)
    stat = unsafe_load(ecos_prob.info)
    solve_time = stat.tsetup + stat.tsolve
    iter = stat.iter
    primal = unsafe_wrap(Array, ecos_prob.x, ecos_prob.n)[:]
    dual_eq = unsafe_wrap(Array, ecos_prob.y, ecos_prob.p)[:]
    dual_ineq = unsafe_wrap(Array, ecos_prob.z, ecos_prob.m)[:]
    slack = unsafe_wrap(Array, ecos_prob.s, ecos_prob.m)[:]
    ECOS.cleanup(inner, 0)
    objective_value = (max_sense ? -1 : 1) * stat.pcost
    dual_objective_value = (max_sense ? -1 : 1) * stat.dcost
    dest.sol = Solution(
        ret_val,
        primal,
        dual_eq,
        dual_ineq,
        slack,
        objective_value,
        dual_objective_value,
        solve_time,
        iter,
    )
    if !MOIU.is_ray(MOI.get(dest, MOI.PrimalStatus()))
        dest.sol.objective_value += objective_constant
    end
    if !MOIU.is_ray(MOI.get(dest, MOI.DualStatus()))
        dest.sol.dual_objective_value += objective_constant
    end
    return
end

function MOI.optimize!(dest::Optimizer, src::OptimizerCache)
    _optimize!(dest, src)
    return MOIU.identity_index_map(src), false
end

function MOI.optimize!(dest::Optimizer, src::MOI.ModelLike)
    cache = OptimizerCache()
    index_map = MOI.copy_to(cache, src)
    _optimize!(dest, cache)
    return index_map, false
end

MOI.get(optimizer::Optimizer, ::MOI.SolveTimeSec) = optimizer.sol.solve_time
function MOI.get(optimizer::Optimizer, ::MOI.BarrierIterations)
    return Int64(optimizer.sol.iter)
end
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
    return optimizer.sol.objective_value
end
function MOI.get(optimizer::Optimizer, attr::MOI.DualObjectiveValue)
    MOI.check_result_index_bounds(optimizer, attr)
    return optimizer.sol.dual_objective_value
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
function MOI.get(
    optimizer::Optimizer,
    attr::MOI.VariablePrimal,
    vi::MOI.VariableIndex,
)
    MOI.check_result_index_bounds(optimizer, attr)
    return optimizer.sol.primal[vi.value]
end
function MOI.get(
    optimizer::Optimizer,
    attr::MOI.ConstraintPrimal,
    ci::MOI.ConstraintIndex{AFF,MOI.Zeros},
)
    MOI.check_result_index_bounds(optimizer, attr)
    return zeros(length(_rows(optimizer, ci)))
end
function MOI.get(
    optimizer::Optimizer,
    attr::MOI.ConstraintPrimal,
    ci::MOI.ConstraintIndex,
)
    MOI.check_result_index_bounds(optimizer, attr)
    return optimizer.sol.slack[_rows(optimizer, ci)]
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
function _dual(optimizer, ci::MOI.ConstraintIndex{AFF,MOI.Zeros})
    return optimizer.sol.dual_eq
end
_dual(optimizer, ci::MOI.ConstraintIndex) = optimizer.sol.dual_ineq
function MOI.get(
    optimizer::Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex{AFF},
)
    MOI.check_result_index_bounds(optimizer, attr)
    return _dual(optimizer, ci)[_rows(optimizer, ci)]
end

MOI.get(optimizer::Optimizer, ::MOI.ResultCount) = 1
