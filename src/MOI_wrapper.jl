using MathOptInterface
const MOI = MathOptInterface
const CI = MOI.ConstraintIndex
const VI = MOI.VariableIndex

const MOIU = MOI.Utilities

const SF = Union{MOI.SingleVariable, MOI.ScalarAffineFunction{Float64}, MOI.VectorOfVariables, MOI.VectorAffineFunction{Float64}}
const SS = Union{MOI.EqualTo{Float64}, MOI.GreaterThan{Float64}, MOI.LessThan{Float64}, MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives, MOI.SecondOrderCone, MOI.ExponentialCone}

struct Solution
    ret_val::Int
    primal::Vector{Float64}
    dual_eq::Vector{Float64}
    dual_ineq::Vector{Float64}
    slack::Vector{Float64}
    objval::Float64
    objbnd::Float64
end
const OPTIMIZE_NOT_CALLED = -1
Solution() = Solution(OPTIMIZE_NOT_CALLED, Float64[], Float64[], Float64[],
                      Float64[], NaN, NaN)

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
    objconstant::Float64 # The objective is min c'x + objconstant
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
    eqsetconstant::Dict{Int, Float64}   # For the constant of EqualTo
    eqnrows::Dict{Int, Int}             # The number of rows of Zeros
    ineqsetconstant::Dict{Int, Float64} # For the constant of LessThan and GreaterThan
    ineqnrows::Dict{Int, Int}           # The number of rows of each vector sets except Zeros
    function ConeData()
        new(0, 0, 0, Int[], 0,
            Dict{Int, Float64}(),
            Dict{Int, UnitRange{Int}}(),
            Dict{Int, Float64}(),
            Dict{Int, UnitRange{Int}}())
    end
end

mutable struct Optimizer <: MOI.AbstractOptimizer
    cone::ConeData
    maxsense::Bool
    data::Union{Nothing, ModelData} # only non-Nothing between MOI.copy_to and MOI.optimize!
    sol::Solution
    options
    function Optimizer(; kwargs...)
        new(ConeData(), false, nothing, Solution(), kwargs)
    end
end

MOI.get(::Optimizer, ::MOI.SolverName) = "ECOS"

function MOI.is_empty(instance::Optimizer)
    !instance.maxsense && instance.data === nothing
end

function MOI.empty!(instance::Optimizer)
    instance.maxsense = false
    instance.data = nothing # It should already be nothing except if an error is thrown inside copy_to
    instance.sol = Solution()
end

MOIU.supports_allocate_load(::Optimizer, copy_names::Bool) = !copy_names

function MOI.supports(::Optimizer,
                      ::Union{MOI.ObjectiveSense,
                              MOI.ObjectiveFunction{MOI.SingleVariable},
                              MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}})
    return true
end

MOI.supports_constraint(::Optimizer, ::Type{<:SF}, ::Type{<:SS}) = true

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; kws...)
    return MOIU.automatic_copy_to(dest, src; kws...)
end

using Compat.SparseArrays

const ZeroCones = Union{MOI.EqualTo, MOI.Zeros}
const LPCones = Union{MOI.GreaterThan, MOI.LessThan, MOI.Nonnegatives, MOI.Nonpositives}

# Computes cone dimensions
constroffset(cone::ConeData, ci::CI{<:MOI.AbstractFunction, <:ZeroCones}) = ci.value
function _allocate_constraint(cone::ConeData, f, s::ZeroCones)
    ci = cone.f
    cone.f += MOI.dimension(s)
    ci
end
constroffset(cone::ConeData, ci::CI{<:MOI.AbstractFunction, <:LPCones}) = ci.value
function _allocate_constraint(cone::ConeData, f, s::LPCones)
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
constroffset(instance::Optimizer, ci::CI) = constroffset(instance.cone, ci::CI)
function MOIU.allocate_constraint(instance::Optimizer, f::F, s::S) where {F <: MOI.AbstractFunction, S <: MOI.AbstractSet}
    CI{F, S}(_allocate_constraint(instance.cone, f, s))
end

# Build constraint matrix
scalecoef(rows, coef, minus, s) = minus ? -coef : coef
scalecoef(rows, coef, minus, s::Union{MOI.LessThan, Type{<:MOI.LessThan}, MOI.Nonpositives, Type{MOI.Nonpositives}}) = minus ? coef : -coef
output_index(t::MOI.VectorAffineTerm) = t.output_index
variable_index_value(t::MOI.ScalarAffineTerm) = t.variable_index.value
variable_index_value(t::MOI.VectorAffineTerm) = variable_index_value(t.scalar_term)
coefficient(t::MOI.ScalarAffineTerm) = t.coefficient
coefficient(t::MOI.VectorAffineTerm) = coefficient(t.scalar_term)
constrrows(::MOI.AbstractScalarSet) = 1
constrrows(s::MOI.AbstractVectorSet) = 1:MOI.dimension(s)
constrrows(instance::Optimizer, ci::CI{<:MOI.AbstractScalarFunction, <:MOI.AbstractScalarSet}) = 1
constrrows(instance::Optimizer, ci::CI{<:MOI.AbstractVectorFunction, MOI.Zeros}) = 1:instance.cone.eqnrows[constroffset(instance, ci)]
constrrows(instance::Optimizer, ci::CI{<:MOI.AbstractVectorFunction, <:MOI.AbstractVectorSet}) = 1:instance.cone.ineqnrows[constroffset(instance, ci)]
matrix(data::ModelData, s::ZeroCones) = data.b, data.IA, data.JA, data.VA
matrix(data::ModelData, s::Union{LPCones, MOI.SecondOrderCone, MOI.ExponentialCone}) = data.h, data.IG, data.JG, data.VG
matrix(instance::Optimizer, s) = matrix(instance.data, s)
MOIU.load_constraint(instance::Optimizer, ci, f::MOI.SingleVariable, s) = MOIU.load_constraint(instance, ci, MOI.ScalarAffineFunction{Float64}(f), s)
function MOIU.load_constraint(instance::Optimizer, ci, f::MOI.ScalarAffineFunction, s::MOI.AbstractScalarSet)
    a = sparsevec(variable_index_value.(f.terms), coefficient.(f.terms))
    # sparsevec combines duplicates with + but does not remove zeros created so we call dropzeros!
    dropzeros!(a)
    offset = constroffset(instance, ci)
    row = constrrows(s)
    i = offset + row
    # The ECOS format is b - Ax ∈ cone
    # so minus=false for b and minus=true for A
    setconstant = MOIU.getconstant(s)
    if s isa MOI.EqualTo
        instance.cone.eqsetconstant[offset] = setconstant
    else
        instance.cone.ineqsetconstant[offset] = setconstant
    end
    constant = f.constant - setconstant
    b, I, J, V = matrix(instance, s)
    b[i] = scalecoef(row, constant, false, s)
    append!(I, fill(i, length(a.nzind)))
    append!(J, a.nzind)
    append!(V, scalecoef(row, a.nzval, true, s))
end
MOIU.load_constraint(instance::Optimizer, ci, f::MOI.VectorOfVariables, s) = MOIU.load_constraint(instance, ci, MOI.VectorAffineFunction{Float64}(f), s)
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
function MOIU.load_constraint(instance::Optimizer, ci, f::MOI.VectorAffineFunction, s::MOI.AbstractVectorSet)
    A = sparse(output_index.(f.terms), variable_index_value.(f.terms), coefficient.(f.terms))
    # sparse combines duplicates with + but does not remove zeros created so we call dropzeros!
    dropzeros!(A)
    I, J, V = findnz(A)
    offset = constroffset(instance, ci)
    rows = constrrows(s)
    if s isa MOI.Zeros
        instance.cone.eqnrows[offset] = length(rows)
    else
        instance.cone.ineqnrows[offset] = length(rows)
    end
    i = offset .+ rows
    # The ECOS format is b - Ax ∈ cone
    # so minus=false for b and minus=true for A
    b, Is, Js, Vs = matrix(instance, s)
    b[i] .= scalecoef(rows, orderval(f.constants, s), false, s)
    append!(Is, offset .+ orderidx(I, s))
    append!(Js, J)
    append!(Vs, scalecoef(I, V, true, s))
end

function MOIU.allocate_variables(instance::Optimizer, nvars::Integer)
    instance.cone = ConeData()
    VI.(1:nvars)
end

function MOIU.load_variables(instance::Optimizer, nvars::Integer)
    cone = instance.cone
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
    instance.data = ModelData(m, nvars, IA, JA, VA, b, IG, JG, VG, h, 0., c)
end

function MOIU.allocate(instance::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    instance.maxsense = sense == MOI.MAX_SENSE
end
function MOIU.allocate(::Optimizer, ::MOI.ObjectiveFunction,
                       ::MOI.Union{MOI.SingleVariable,
                                   MOI.ScalarAffineFunction{Float64}})
end

function MOIU.load(::Optimizer, ::MOI.ObjectiveSense, ::MOI.OptimizationSense)
end
function MOIU.load(optimizer::Optimizer, ::MOI.ObjectiveFunction,
                   f::MOI.SingleVariable)
    MOIU.load(optimizer,
              MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
              MOI.ScalarAffineFunction{Float64}(f))
end
function MOIU.load(instance::Optimizer, ::MOI.ObjectiveFunction,
                   f::MOI.ScalarAffineFunction)
    c0 = Vector(sparsevec(variable_index_value.(f.terms), coefficient.(f.terms),
                          instance.data.n))
    instance.data.objconstant = f.constant
    instance.data.c = instance.maxsense ? -c0 : c0
    return nothing
end

function MOI.optimize!(instance::Optimizer)
    if instance.data === nothing
        # optimize! has already been called and no new model has been copied
        return
    end
    cone = instance.cone
    m = instance.data.m
    n = instance.data.n
    A = ECOS.ECOSMatrix(sparse(instance.data.IA, instance.data.JA, instance.data.VA, cone.f, n))
    b = instance.data.b
    G = ECOS.ECOSMatrix(sparse(instance.data.IG, instance.data.JG, instance.data.VG, m, n))
    h = instance.data.h
    objconstant = instance.data.objconstant
    c = instance.data.c
    instance.data = nothing # Allows GC to free instance.data before A is loaded to ECOS
    ecos_prob_ptr = ECOS.setup(n, m, cone.f, cone.l, length(cone.qa), cone.qa,
                               cone.ep, G, A, c, h, b; instance.options...)
    ret_val = ECOS.solve(ecos_prob_ptr)
    ecos_prob = unsafe_wrap(Array, ecos_prob_ptr, 1)[1]
    primal    = unsafe_wrap(Array, ecos_prob.x, n)[:]
    dual_eq   = unsafe_wrap(Array, ecos_prob.y, cone.f)[:]
    dual_ineq = unsafe_wrap(Array, ecos_prob.z, m)[:]
    slack     = unsafe_wrap(Array, ecos_prob.s, m)[:]
    ECOS.cleanup(ecos_prob_ptr, 0)
    objval = (instance.maxsense ? -1 : 1) * dot(c, primal)
    if ret_val != ECOS.ECOS_DINF
        objval += objconstant
    end
    objbnd = -(dot(b, dual_eq) + dot(h, dual_ineq))
    if ret_val != ECOS.ECOS_PINF
        objbnd += objconstant
    end
    instance.sol = Solution(ret_val, primal, dual_eq, dual_ineq, slack, objval,
                            objbnd)
end

# Implements getter for result value and statuses
function MOI.get(instance::Optimizer, ::MOI.TerminationStatus)
    flag = instance.sol.ret_val
    if flag == OPTIMIZE_NOT_CALLED
        return MOI.OPTIMIZE_NOT_CALLED
    elseif flag == ECOS.ECOS_OPTIMAL
        return MOI.OPTIMAL
    elseif flag == ECOS.ECOS_PINF
        return MOI.INFEASIBLE
    elseif flag == ECOS.ECOS_DINF
        return MOI.DUAL_INFEASIBLE
    elseif flag == ECOS.ECOS_MAXIT
        return MOI.IterationLimit
    elseif flag == ECOS.ECOS_OPTIMAL + ECOS.ECOS_INACC_OFFSET
        return MOI.ALMOST_OPTIMAL
    elseif flag == ECOS.ECOS_PINF + ECOS.ECOS_INACC_OFFSET
        return MOI.ALMOST_INFEASIBLE
    elseif flag == ECOS.ECOS_DINF + ECOS.ECOS_INACC_OFFSET
        return MOI.ALMOST_DUAL_INFEASIBLE
    else
        return MOI.OTHER_ERROR
    end
end

MOI.get(instance::Optimizer, ::MOI.ObjectiveValue) = instance.sol.objval
MOI.get(instance::Optimizer, ::MOI.ObjectiveBound) = instance.sol.objbnd

function MOI.get(instance::Optimizer, ::MOI.PrimalStatus)
    flag = instance.sol.ret_val
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
function MOI.get(instance::Optimizer, ::MOI.VariablePrimal, vi::VI)
    instance.sol.primal[vi.value]
end
MOI.get(instance::Optimizer, a::MOI.VariablePrimal, vi::Vector{VI}) = MOI.get.(instance, Ref(a), vi)
# setconstant: Retrieve set constant stored in `ConeData` during `copy_to`
setconstant(instance::Optimizer, offset, ::CI{<:MOI.AbstractFunction, <:MOI.EqualTo}) = instance.cone.eqsetconstant[offset]
setconstant(instance::Optimizer, offset, ::CI) = instance.cone.ineqsetconstant[offset]
_unshift(instance::Optimizer, offset, value, ::CI) = value
_unshift(instance::Optimizer, offset, value, ci::CI{<:MOI.AbstractScalarFunction, <:MOI.AbstractScalarSet}) = value + setconstant(instance, offset, ci)
function MOI.get(instance::Optimizer, ::MOI.ConstraintPrimal, ci::CI{<:MOI.AbstractFunction, MOI.Zeros})
    rows = constrrows(instance, ci)
    zeros(length(rows))
end
function MOI.get(instance::Optimizer, ::MOI.ConstraintPrimal, ci::CI{<:MOI.AbstractFunction, <:MOI.EqualTo})
    offset = constroffset(instance, ci)
    setconstant(instance, offset, ci)
end
function MOI.get(instance::Optimizer, ::MOI.ConstraintPrimal, ci::CI{<:MOI.AbstractFunction, S}) where S <: MOI.AbstractSet
    offset = constroffset(instance, ci)
    rows = constrrows(instance, ci)
    _unshift(instance, offset, scalecoef(rows, reorderval(instance.sol.slack[offset .+ rows], S), false, S), ci)
end

function MOI.get(instance::Optimizer, ::MOI.DualStatus)
    flag = instance.sol.ret_val
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
_dual(instance, ci::CI{<:MOI.AbstractFunction, <:ZeroCones}) = instance.sol.dual_eq
_dual(instance, ci::CI) = instance.sol.dual_ineq
function MOI.get(instance::Optimizer, ::MOI.ConstraintDual, ci::CI{<:MOI.AbstractFunction, S}) where S <: MOI.AbstractSet
    offset = constroffset(instance, ci)
    rows = constrrows(instance, ci)
    scalecoef(rows, reorderval(_dual(instance, ci)[offset .+ rows], S), false, S)
end

MOI.get(instance::Optimizer, ::MOI.ResultCount) = 1
