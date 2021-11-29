using MathOptInterface
const MOI = MathOptInterface

include("permuted_exponential_cone.jl")

MOI.Utilities.@product_of_sets(Zeros, MOI.Zeros)

MOI.Utilities.@product_of_sets(
    Cones,
    MOI.Nonnegatives,
    MOI.SecondOrderCone,
    PermutedExponentialCone
)

MOI.Utilities.@struct_of_constraints_by_set_types(
    ZerosOrNot,
    MOI.Zeros,
    Union{MOI.Nonnegatives,MOI.SecondOrderCone,PermutedExponentialCone},
)

# !!! danger
#     Make sure to use ECOS's definition of pfloat for double and idxint for
#     int. These are _not_ the same across platforms, differ from Julia's
#     Cdouble and Clong, and even differ across patch releases of ECOS :(.

const OptimizerCache = MOI.Utilities.GenericModel{
    pfloat,
    MOI.Utilities.ObjectiveContainer{pfloat},
    MOI.Utilities.VariablesContainer{pfloat},
    ZerosOrNot{pfloat}{
        MOI.Utilities.MatrixOfConstraints{
            pfloat,
            MOI.Utilities.MutableSparseMatrixCSC{
                pfloat,
                idxint,
                MOI.Utilities.ZeroBasedIndexing,
            },
            Vector{pfloat},
            Zeros{pfloat},
        },
        MOI.Utilities.MatrixOfConstraints{
            pfloat,
            MOI.Utilities.MutableSparseMatrixCSC{
                pfloat,
                idxint,
                MOI.Utilities.ZeroBasedIndexing,
            },
            Vector{pfloat},
            Cones{pfloat},
        },
    },
}

mutable struct _Solution
    ret_val::Union{Nothing,idxint}
    primal::Vector{pfloat}
    dual_eq::Vector{pfloat}
    dual_ineq::Vector{pfloat}
    slack::Vector{pfloat}
    objective_value::pfloat
    dual_objective_value::pfloat
    solve_time::pfloat
    iter::idxint
end

# The values used by ECOS are from -7 to 10 so -10 should be safe
function _Solution()
    return _Solution(
        nothing,
        pfloat[],
        pfloat[],
        pfloat[],
        pfloat[],
        NaN,
        NaN,
        NaN,
        0,
    )
end

mutable struct Optimizer <: MOI.AbstractOptimizer
    zeros::Union{Nothing,Zeros{pfloat}}
    cones::Union{Nothing,Cones{pfloat}}
    sol::_Solution
    silent::Bool
    options::Dict{Symbol,Any}
    function Optimizer(; kwargs...)
        return new(nothing, nothing, _Solution(), false, Dict{Symbol,Any}())
    end
end

function MOI.is_empty(optimizer::Optimizer)
    return optimizer.zeros === nothing && optimizer.cones === nothing
end

function MOI.empty!(optimizer::Optimizer)
    optimizer.zeros = nothing
    optimizer.cones = nothing
    optimizer.sol = _Solution()
    return
end

function MOI.get(::Optimizer, ::MOI.Bridges.ListOfNonstandardBridges)
    return [PermutedExponentialBridge{pfloat}]
end

function _rows(
    optimizer::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VectorAffineFunction{pfloat},MOI.Zeros},
)
    return MOI.Utilities.rows(optimizer.zeros, ci)
end

function _rows(
    optimizer::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VectorAffineFunction{pfloat}},
)
    return MOI.Utilities.rows(optimizer.cones, ci)
end

MOI.get(::Optimizer, ::MOI.SolverName) = "ECOS"

MOI.get(::Optimizer, ::MOI.SolverVersion) = unsafe_string(ECOS_ver())

# MOI.RawOptimizerAttribute

function MOI.supports(::Optimizer, param::MOI.RawOptimizerAttribute)
    return hasfield(settings, Symbol(param.name))
end

function MOI.set(optimizer::Optimizer, param::MOI.RawOptimizerAttribute, value)
    optimizer.options[Symbol(param.name)] = value
    return
end

function MOI.get(optimizer::Optimizer, param::MOI.RawOptimizerAttribute)
    # TODO(odow): This gives a poor error message if the name of the parameter
    # is invalid.
    return optimizer.options[Symbol(param.name)]
end

# MOI.Silent

MOI.supports(::Optimizer, ::MOI.Silent) = true

function MOI.set(optimizer::Optimizer, ::MOI.Silent, value::Bool)
    optimizer.silent = value
    return
end

MOI.get(optimizer::Optimizer, ::MOI.Silent) = optimizer.silent

# MOI.supports

function MOI.supports(
    ::Optimizer,
    ::Union{
        MOI.ObjectiveSense,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{pfloat}},
    },
)
    return true
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VectorAffineFunction{pfloat}},
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
    Ab = MOI.Utilities.constraints(
        src.constraints,
        MOI.VectorAffineFunction{pfloat},
        MOI.Zeros,
    )
    A = Ab.coefficients
    Gh = MOI.Utilities.constraints(
        src.constraints,
        MOI.VectorAffineFunction{pfloat},
        MOI.Nonnegatives,
    )
    G = Gh.coefficients
    soc_indices = MOI.get(
        Gh,
        MOI.ListOfConstraintIndices{
            MOI.VectorAffineFunction{pfloat},
            MOI.SecondOrderCone,
        }(),
    )
    q = MOI.dimension.(MOI.get.(Gh, MOI.ConstraintSet(), soc_indices))
    @assert A.n == G.n
    max_sense = MOI.get(src, MOI.ObjectiveSense()) == MOI.MAX_SENSE
    obj =
        MOI.get(src, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{pfloat}}())
    objective_constant = MOI.constant(obj)
    c0 = zeros(A.n)
    for term in obj.terms
        c0[term.variable.value] += term.coefficient
    end
    dest.zeros = deepcopy(Ab.sets) # TODO copy(Ab.sets)
    dest.cones = deepcopy(Gh.sets) # TODO copy(Gh.sets)
    options = copy(dest.options)
    if dest.silent
        options[:verbose] = false
    end
    num_exponential = MOI.get(
        Gh,
        MOI.NumberOfConstraints{
            MOI.VectorAffineFunction{pfloat},
            PermutedExponentialCone,
        }(),
    )
    inner = ECOS_setup(
        A.n,
        G.m,
        A.m,
        Gh.sets.num_rows[1],
        length(q),
        q,
        num_exponential,
        -G.nzval,
        G.colptr,
        G.rowval,
        -A.nzval,
        A.colptr,
        A.rowval,
        max_sense ? -c0 : c0,
        Gh.constants,
        Ab.constants,
    )
    if inner == C_NULL
        error("ECOS failed to construct problem.")
    end
    unsafe_add_settings(inner, options)
    ret_val = ECOS_solve(inner)
    ecos_prob = unsafe_load(inner)::pwork
    stat = unsafe_load(ecos_prob.info)::stats
    dest.sol = _Solution(
        ret_val,
        copy(unsafe_wrap(Array, ecos_prob.x, ecos_prob.n)),
        copy(unsafe_wrap(Array, ecos_prob.y, ecos_prob.p)),
        copy(unsafe_wrap(Array, ecos_prob.z, ecos_prob.m)),
        copy(unsafe_wrap(Array, ecos_prob.s, ecos_prob.m)),
        max_sense ? -stat.pcost : stat.pcost,
        max_sense ? -stat.dcost : stat.dcost,
        stat.tsetup + stat.tsolve,
        stat.iter,
    )
    ECOS_cleanup(inner, 0)
    if !MOI.Utilities.is_ray(MOI.get(dest, MOI.PrimalStatus()))
        dest.sol.objective_value += objective_constant
    end
    if !MOI.Utilities.is_ray(MOI.get(dest, MOI.DualStatus()))
        dest.sol.dual_objective_value += objective_constant
    end
    return
end

function MOI.optimize!(dest::Optimizer, src::OptimizerCache)
    _optimize!(dest, src)
    return MOI.Utilities.identity_index_map(src), false
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
    elseif flag == ECOS_OPTIMAL
        return MOI.OPTIMAL
    elseif flag == ECOS_PINF
        return MOI.INFEASIBLE
    elseif flag == ECOS_DINF
        return MOI.DUAL_INFEASIBLE
    elseif flag == ECOS_MAXIT
        return MOI.ITERATION_LIMIT
    elseif flag == ECOS_NUMERICS || flag == ECOS_OUTCONE
        return MOI.NUMERICAL_ERROR
    elseif flag == ECOS_SIGINT
        return MOI.INTERRUPTED
    elseif flag == ECOS_FATAL
        return MOI.OTHER_ERROR
    elseif flag == ECOS_OPTIMAL + ECOS_INACC_OFFSET
        return MOI.ALMOST_OPTIMAL
    elseif flag == ECOS_PINF + ECOS_INACC_OFFSET
        return MOI.ALMOST_INFEASIBLE
    else
        @assert flag == ECOS_DINF + ECOS_INACC_OFFSET
        return MOI.ALMOST_DUAL_INFEASIBLE
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
    if flag == ECOS_OPTIMAL
        return MOI.FEASIBLE_POINT
    elseif flag == ECOS_PINF
        return MOI.INFEASIBLE_POINT
    elseif flag == ECOS_DINF
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif flag == ECOS_MAXIT
        return MOI.UNKNOWN_RESULT_STATUS
    elseif flag == ECOS_OPTIMAL + ECOS_INACC_OFFSET
        return MOI.NEARLY_FEASIBLE_POINT
    elseif flag == ECOS_PINF + ECOS_INACC_OFFSET
        return MOI.INFEASIBLE_POINT
    elseif flag == ECOS_DINF + ECOS_INACC_OFFSET
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
    ci::MOI.ConstraintIndex{MOI.VectorAffineFunction{pfloat},MOI.Zeros},
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
    if flag == ECOS_OPTIMAL
        return MOI.FEASIBLE_POINT
    elseif flag == ECOS_PINF
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif flag == ECOS_DINF
        return MOI.INFEASIBLE_POINT
    elseif flag == ECOS_MAXIT
        return MOI.UNKNOWN_RESULT_STATUS
    elseif flag == ECOS_OPTIMAL + ECOS_INACC_OFFSET
        return MOI.NEARLY_FEASIBLE_POINT
    elseif flag == ECOS_PINF + ECOS_INACC_OFFSET
        return MOI.NEARLY_INFEASIBILITY_CERTIFICATE
    elseif flag == ECOS_DINF + ECOS_INACC_OFFSET
        return MOI.INFEASIBLE_POINT
    end
    return MOI.OTHER_RESULT_STATUS
end

function MOI.get(
    optimizer::Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex{MOI.VectorAffineFunction{pfloat},MOI.Zeros},
)
    MOI.check_result_index_bounds(optimizer, attr)
    return optimizer.sol.dual_eq[_rows(optimizer, ci)]
end

function MOI.get(
    optimizer::Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex{MOI.VectorAffineFunction{pfloat}},
)
    MOI.check_result_index_bounds(optimizer, attr)
    return optimizer.sol.dual_ineq[_rows(optimizer, ci)]
end

MOI.get(optimizer::Optimizer, ::MOI.ResultCount) = 1