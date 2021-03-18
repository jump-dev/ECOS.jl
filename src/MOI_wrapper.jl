using MathOptInterface
const MOI = MathOptInterface
const CI = MOI.ConstraintIndex
const VI = MOI.VariableIndex

const MOIU = MOI.Utilities

struct Solution
    ret_val::Union{Nothing, Int}
    primal::Vector{Float64}
    dual_eq::Vector{Float64}
    dual_ineq::Vector{Float64}
    slack::Vector{Float64}
    objective_value::Float64
    dual_objective_value::Float64
    objective_constant::Float64
    solve_time::Float64
end
# The values used by ECOS are from -7 to 10 so -10 should be safe
Solution() = Solution(nothing, Float64[], Float64[], Float64[],
                      Float64[], NaN, NaN, NaN, NaN)

# Used to build the data with allocate-load during `copy_to`.
# When `optimize!` is called, a the data is used to build `ECOSMatrix`
# and the `ModelData` struct is discarded
mutable struct ModelData
    m::Int # Number of rows/constraints
    n::Int # Number of cols/variables
    IA::Vector{Int} # List of conic rows
    JA::Vector{Int} # List of conic cols
    VA::Vector{Float64} # List of conic coefficients
    b::Vector{Float64} # List of conic coefficients
    IG::Vector{Int} # List of equality rows
    JG::Vector{Int} # List of equality cols
    VG::Vector{Float64} # List of equality coefficients
    h::Vector{Float64} # List of equality coefficients
    objective_constant::Float64 # The objective is min c'x + objective_constant
    c::Vector{Float64}
end

# This is tied to ECOS's internal representation
mutable struct ConeData
    f::Int # number of linear equality constraints
    l::Int # length of LP cone
    q::Int # length of SOC cone
    qa::Vector{Int} # array of second-order cone constraints
    ep::Int # number of primal exponential cone triples
    # The following four field store model information needed to compute `ConstraintPrimal` and `ConstraintDual`
    eqnrows::Dict{Int, Int}   # The number of rows of Zeros
    ineqnrows::Dict{Int, Int} # The number of rows of each vector sets except Zeros
    function ConeData()
        new(0, 0, 0, Int[], 0,
            Dict{Int, UnitRange{Int}}(),
            Dict{Int, UnitRange{Int}}())
    end
end

mutable struct Optimizer <: MOI.AbstractOptimizer
    cone::ConeData
    maxsense::Bool
    data::Union{Nothing, ModelData} # only non-Nothing between MOI.copy_to and MOI.optimize!
    sol::Solution
    silent::Bool
    options::Dict{Symbol, Any}
    function Optimizer(; kwargs...)
        optimizer = new(ConeData(), false, nothing, Solution(), false, Dict{Symbol, Any}())
        for (key, value) in kwargs
            MOI.set(optimizer, MOI.RawParameter(String(key)), value)
        end
        return optimizer
    end
end

MOI.get(::Optimizer, ::MOI.SolverName) = "ECOS"

function MOI.set(optimizer::Optimizer, param::MOI.RawParameter, value)
    if !(param.name isa String)
        Base.depwarn(
            "passing `$(param.name)` to `MOI.RawParameter` as type " *
            "`$(typeof(param.name))` is deprecated. Use a string instead.",
            Symbol("MOI.set")
        )
    end
    optimizer.options[Symbol(param.name)] = value
end
function MOI.get(optimizer::Optimizer, param::MOI.RawParameter)
    if !(param.name isa String)
        Base.depwarn(
            "passing $(param.name) to `MOI.RawParameter` as type " *
            "$(typeof(param.name)) is deprecated. Use a string instead.",
            Symbol("MOI.get")
        )
    end
    # TODO: This gives a poor error message if the name of the parameter is invalid.
    return optimizer.options[Symbol(param.name)]
end

MOI.supports(::Optimizer, ::MOI.Silent) = true
function MOI.set(optimizer::Optimizer, ::MOI.Silent, value::Bool)
    optimizer.silent = value
end
MOI.get(optimizer::Optimizer, ::MOI.Silent) = optimizer.silent

function MOI.is_empty(optimizer::Optimizer)
    !optimizer.maxsense && optimizer.data === nothing
end

function MOI.empty!(optimizer::Optimizer)
    optimizer.maxsense = false
    optimizer.data = nothing # It should already be nothing except if an error is thrown inside copy_to
    optimizer.sol = Solution()
end

MOIU.supports_allocate_load(::Optimizer, copy_names::Bool) = !copy_names

function MOI.supports(::Optimizer,
                      ::Union{MOI.ObjectiveSense,
                              MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}})
    return true
end

function MOI.supports_constraint(::Optimizer,
                                 ::Type{MOI.VectorAffineFunction{Float64}},
                                 ::Type{<:Union{MOI.Zeros, MOI.Nonnegatives,
                                                MOI.SecondOrderCone,
                                                MOI.ExponentialCone}})
    return true
end

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; kws...)
    return MOIU.automatic_copy_to(dest, src; kws...)
end

using SparseArrays

# Computes cone dimensions
constroffset(cone::ConeData, ci::CI{<:MOI.AbstractFunction, MOI.Zeros}) = ci.value
function _allocate_constraint(cone::ConeData, f, s::MOI.Zeros)
    ci = cone.f
    cone.f += MOI.dimension(s)
    ci
end
constroffset(cone::ConeData, ci::CI{<:MOI.AbstractFunction, MOI.Nonnegatives}) = ci.value
function _allocate_constraint(cone::ConeData, f, s::MOI.Nonnegatives)
    ci = cone.l
    cone.l += MOI.dimension(s)
    ci
end
constroffset(cone::ConeData, ci::CI{<:MOI.AbstractFunction, <:MOI.SecondOrderCone}) = cone.l + ci.value
function _allocate_constraint(cone::ConeData, f, s::MOI.SecondOrderCone)
    push!(cone.qa, s.dimension)
    ci = cone.q
    cone.q += MOI.dimension(s)
    ci
end
constroffset(cone::ConeData, ci::CI{<:MOI.AbstractFunction, <:MOI.ExponentialCone}) = cone.l + cone.q + ci.value
function _allocate_constraint(cone::ConeData, f, s::MOI.ExponentialCone)
    ci = 3cone.ep
    cone.ep += 1
    ci
end
constroffset(optimizer::Optimizer, ci::CI) = constroffset(optimizer.cone, ci::CI)
function MOIU.allocate_constraint(optimizer::Optimizer, f::F, s::S) where {F <: MOI.AbstractFunction, S <: MOI.AbstractSet}
    CI{F, S}(_allocate_constraint(optimizer.cone, f, s))
end

# Build constraint matrix
output_index(t::MOI.VectorAffineTerm) = t.output_index
variable_index_value(t::MOI.ScalarAffineTerm) = t.variable_index.value
variable_index_value(t::MOI.VectorAffineTerm) = variable_index_value(t.scalar_term)
coefficient(t::MOI.ScalarAffineTerm) = t.coefficient
coefficient(t::MOI.VectorAffineTerm) = coefficient(t.scalar_term)
constrrows(s::MOI.AbstractVectorSet) = 1:MOI.dimension(s)
constrrows(optimizer::Optimizer, ci::CI{<:MOI.AbstractVectorFunction, MOI.Zeros}) = 1:optimizer.cone.eqnrows[constroffset(optimizer, ci)]
constrrows(optimizer::Optimizer, ci::CI{<:MOI.AbstractVectorFunction, <:MOI.AbstractVectorSet}) = 1:optimizer.cone.ineqnrows[constroffset(optimizer, ci)]
matrix(data::ModelData, s::MOI.Zeros) = data.b, data.IA, data.JA, data.VA
matrix(data::ModelData, s::Union{MOI.Nonnegatives, MOI.SecondOrderCone, MOI.ExponentialCone}) = data.h, data.IG, data.JG, data.VG
matrix(optimizer::Optimizer, s) = matrix(optimizer.data, s)
# ECOS orders differently than MOI the second and third dimension of the exponential cone
orderval(val, s) = val
function orderval(val, s::Union{MOI.ExponentialCone, Type{MOI.ExponentialCone}})
    val[[1, 3, 2]]
end
orderidx(idx, s) = idx
expmap(i) = (1, 3, 2)[i]
function orderidx(idx, s::MOI.ExponentialCone)
    expmap.(idx)
end
function MOIU.load_constraint(optimizer::Optimizer, ci::MOI.ConstraintIndex, f::MOI.VectorAffineFunction, s::MOI.AbstractVectorSet)
    func = MOIU.canonical(f)
    I = Int[output_index(term) for term in func.terms]
    J = Int[variable_index_value(term) for term in func.terms]
    V = Float64[-coefficient(term) for term in func.terms]
    offset = constroffset(optimizer, ci)
    rows = constrrows(s)
    if s isa MOI.Zeros
        optimizer.cone.eqnrows[offset] = length(rows)
    else
        optimizer.cone.ineqnrows[offset] = length(rows)
    end
    i = offset .+ rows
    # The ECOS format is b - Ax ∈ cone
    # so minus=false for b and minus=true for A
    b, Is, Js, Vs = matrix(optimizer, s)
    b[i] .= orderval(f.constants, s)
    append!(Is, offset .+ orderidx(I, s))
    append!(Js, J)
    append!(Vs, V)
end

function MOIU.allocate_variables(optimizer::Optimizer, nvars::Integer)
    optimizer.cone = ConeData()
    VI.(1:nvars)
end

function MOIU.load_variables(optimizer::Optimizer, nvars::Integer)
    cone = optimizer.cone
    m = cone.l + cone.q + 3cone.ep
    IA = Int[]
    JA = Int[]
    VA = Float64[]
    b = zeros(cone.f)
    IG = Int[]
    JG = Int[]
    VG = Float64[]
    h = zeros(m)
    c = zeros(nvars)
    optimizer.data = ModelData(m, nvars, IA, JA, VA, b, IG, JG, VG, h, 0., c)
end

function MOIU.allocate(optimizer::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    optimizer.maxsense = sense == MOI.MAX_SENSE
end
function MOIU.allocate(::Optimizer, ::MOI.ObjectiveFunction,
                       ::MOI.ScalarAffineFunction{Float64})
end

function MOIU.load(::Optimizer, ::MOI.ObjectiveSense, ::MOI.OptimizationSense)
end
function MOIU.load(optimizer::Optimizer, ::MOI.ObjectiveFunction,
                   f::MOI.ScalarAffineFunction)
    c0 = Vector(sparsevec(variable_index_value.(f.terms), coefficient.(f.terms),
                          optimizer.data.n))
    optimizer.data.objective_constant = f.constant
    optimizer.data.c = optimizer.maxsense ? -c0 : c0
    return nothing
end

function MOI.optimize!(optimizer::Optimizer)
    if optimizer.data === nothing
        # optimize! has already been called and no new model has been copied
        return
    end
    cone = optimizer.cone
    m = optimizer.data.m
    n = optimizer.data.n
    A = ECOS.ECOSMatrix(sparse(optimizer.data.IA, optimizer.data.JA, optimizer.data.VA, cone.f, n))
    b = optimizer.data.b
    G = ECOS.ECOSMatrix(sparse(optimizer.data.IG, optimizer.data.JG, optimizer.data.VG, m, n))
    h = optimizer.data.h
    objective_constant = optimizer.data.objective_constant
    c = optimizer.data.c
    optimizer.data = nothing # Allows GC to free optimizer.data before A is loaded to ECOS
    options = optimizer.options
    if optimizer.silent
        options = copy(options)
        options[:verbose] = false
    end
    ecos_prob_p = setup(
        n,
        m,
        cone.f,
        cone.l,
        length(cone.qa),
        cone.qa,
        cone.ep,
        G,
        A,
        c,
        h,
        b;
        options...,
    )
    ret_val = ECOS_solve(ecos_prob_p)
    ecos_prob = unsafe_load(ecos_prob_p.ptr)
    stat = unsafe_load(ecos_prob.info)
    solve_time = stat.tsetup + stat.tsolve
    primal    = copy(unsafe_wrap(Array, ecos_prob.x, n))
    dual_eq   = copy(unsafe_wrap(Array, ecos_prob.y, cone.f))
    dual_ineq = copy(unsafe_wrap(Array, ecos_prob.z, m))
    slack     = copy(unsafe_wrap(Array, ecos_prob.s, m))
    # ECOS_cleanup(ecos_prob_ptr, 0)
    objective_value = (optimizer.maxsense ? -1 : 1) * stat.pcost
    dual_objective_value = (optimizer.maxsense ? -1 : 1) * stat.dcost
    optimizer.sol = Solution(ret_val, primal, dual_eq, dual_ineq, slack, objective_value,
                            dual_objective_value, objective_constant, solve_time)
    return
end

MOI.get(optimizer::Optimizer, ::MOI.SolveTime) = optimizer.sol.solve_time
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
    if attr.N > MOI.get(optimizer, MOI.ResultCount())
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
    optimizer.sol.primal[vi.value]
end
MOI.get(optimizer::Optimizer, a::MOI.VariablePrimal, vi::Vector{VI}) = MOI.get.(optimizer, Ref(a), vi)
function MOI.get(optimizer::Optimizer, attr::MOI.ConstraintPrimal, ci::CI{<:MOI.AbstractFunction, MOI.Zeros})
    MOI.check_result_index_bounds(optimizer, attr)
    rows = constrrows(optimizer, ci)
    return zeros(length(rows))
end
function MOI.get(optimizer::Optimizer, attr::MOI.ConstraintPrimal, ci::CI{<:MOI.AbstractFunction, S}) where S <: MOI.AbstractSet
    MOI.check_result_index_bounds(optimizer, attr)
    offset = constroffset(optimizer, ci)
    rows = constrrows(optimizer, ci)
    return reorderval(optimizer.sol.slack[offset .+ rows], S)
end

function MOI.get(optimizer::Optimizer, attr::MOI.DualStatus)
    if attr.N > MOI.get(optimizer, MOI.ResultCount())
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
_dual(optimizer, ci::CI{<:MOI.AbstractFunction, <:MOI.Zeros}) = optimizer.sol.dual_eq
_dual(optimizer, ci::CI) = optimizer.sol.dual_ineq
function MOI.get(optimizer::Optimizer, attr::MOI.ConstraintDual, ci::CI{<:MOI.AbstractFunction, S}) where S <: MOI.AbstractSet
    MOI.check_result_index_bounds(optimizer, attr)
    offset = constroffset(optimizer, ci)
    rows = constrrows(optimizer, ci)
    return reorderval(_dual(optimizer, ci)[offset .+ rows], S)
end

MOI.get(optimizer::Optimizer, ::MOI.ResultCount) = 1
