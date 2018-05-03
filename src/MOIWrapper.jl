export ECOSOptimizer

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
end
Solution() = Solution(0, Float64[], Float64[], Float64[], Float64[], NaN)

# Used to build the data with allocate-load during `copy!`.
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

mutable struct ECOSOptimizer <: MOI.AbstractOptimizer
    cone::ConeData
    maxsense::Bool
    data::Union{Void, ModelData} # only non-Void between MOI.copy! and MOI.optimize!
    sol::Solution
    function ECOSOptimizer()
        new(ConeData(), false, nothing, Solution())
    end
end

function MOI.isempty(instance::ECOSOptimizer)
    !instance.maxsense && instance.data === nothing
end

function MOI.empty!(instance::ECOSOptimizer)
    instance.maxsense = false
    instance.data = nothing # It should already be nothing except if an error is thrown inside copy!
end

MOI.canaddvariable(instance::ECOSOptimizer) = false
MOIU.needsallocateload(instance::ECOSOptimizer) = true
MOI.supports(::ECOSOptimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}) = true
MOI.supportsconstraint(::ECOSOptimizer, ::Type{<:SF}, ::Type{<:SS}) = true
MOI.copy!(dest::ECOSOptimizer, src::MOI.ModelLike; copynames=true) = MOIU.allocateload!(dest, src, copynames)

using Compat.SparseArrays

const ZeroCones = Union{MOI.EqualTo, MOI.Zeros}
const LPCones = Union{MOI.GreaterThan, MOI.LessThan, MOI.Nonnegatives, MOI.Nonpositives}

_dim(s::MOI.AbstractScalarSet) = 1
_dim(s::MOI.AbstractVectorSet) = MOI.dimension(s)

# Computes cone dimensions
constroffset(cone::ConeData, ci::CI{<:MOI.AbstractFunction, <:ZeroCones}) = ci.value
function _allocateconstraint!(cone::ConeData, f, s::ZeroCones)
    ci = cone.f
    cone.f += _dim(s)
    ci
end
constroffset(cone::ConeData, ci::CI{<:MOI.AbstractFunction, <:LPCones}) = ci.value
function _allocateconstraint!(cone::ConeData, f, s::LPCones)
    ci = cone.l
    cone.l += _dim(s)
    ci
end
constroffset(cone::ConeData, ci::CI{<:MOI.AbstractFunction, <:MOI.SecondOrderCone}) = cone.l + ci.value
function _allocateconstraint!(cone::ConeData, f, s::MOI.SecondOrderCone)
    push!(cone.qa, s.dimension)
    ci = cone.q
    cone.q += _dim(s)
    ci
end
constroffset(cone::ConeData, ci::CI{<:MOI.AbstractFunction, <:MOI.ExponentialCone}) = cone.l + cone.q + ci.value
function _allocateconstraint!(cone::ConeData, f, s::MOI.ExponentialCone)
    ci = 3cone.ep
    cone.ep += 1
    ci
end
constroffset(instance::ECOSOptimizer, ci::CI) = constroffset(instance.cone, ci::CI)
MOIU.canallocateconstraint(::ECOSOptimizer, ::Type{<:SF}, ::Type{<:SS}) = true
function MOIU.allocateconstraint!(instance::ECOSOptimizer, f::F, s::S) where {F <: MOI.AbstractFunction, S <: MOI.AbstractSet}
    CI{F, S}(_allocateconstraint!(instance.cone, f, s))
end

# Build constraint matrix
scalecoef(rows, coef, minus, s) = minus ? -coef : coef
scalecoef(rows, coef, minus, s::Union{MOI.LessThan, Type{<:MOI.LessThan}, MOI.Nonpositives, Type{MOI.Nonpositives}}) = minus ? coef : -coef
_varmap(f) = map(vi -> vi.value, f.variables)
constrrows(::MOI.AbstractScalarSet) = 1
constrrows(s::MOI.AbstractVectorSet) = 1:MOI.dimension(s)
constrrows(instance::ECOSOptimizer, ci::CI{<:MOI.AbstractScalarFunction, <:MOI.AbstractScalarSet}) = 1
constrrows(instance::ECOSOptimizer, ci::CI{<:MOI.AbstractVectorFunction, MOI.Zeros}) = 1:instance.cone.eqnrows[constroffset(instance, ci)]
constrrows(instance::ECOSOptimizer, ci::CI{<:MOI.AbstractVectorFunction, <:MOI.AbstractVectorSet}) = 1:instance.cone.ineqnrows[constroffset(instance, ci)]
matrix(data::ModelData, s::ZeroCones) = data.b, data.IA, data.JA, data.VA
matrix(data::ModelData, s::Union{LPCones, MOI.SecondOrderCone, MOI.ExponentialCone}) = data.h, data.IG, data.JG, data.VG
matrix(instance::ECOSOptimizer, s) = matrix(instance.data, s)
MOIU.canloadconstraint(::ECOSOptimizer, ::Type{<:SF}, ::Type{<:SS}) = true
MOIU.loadconstraint!(instance::ECOSOptimizer, ci, f::MOI.SingleVariable, s) = MOIU.loadconstraint!(instance, ci, MOI.ScalarAffineFunction{Float64}(f), s)
function MOIU.loadconstraint!(instance::ECOSOptimizer, ci, f::MOI.ScalarAffineFunction, s::MOI.AbstractScalarSet)
    a = sparsevec(_varmap(f), f.coefficients)
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
MOIU.loadconstraint!(instance::ECOSOptimizer, ci, f::MOI.VectorOfVariables, s) = MOIU.loadconstraint!(instance, ci, MOI.VectorAffineFunction{Float64}(f), s)
# SCS orders differently than MOI the second and third dimension of the exponential cone
orderval(val, s) = val
function orderval(val, s::Union{MOI.ExponentialCone, Type{MOI.ExponentialCone}})
    val[[1, 3, 2]]
end
orderidx(idx, s) = idx
expmap(i) = (1, 3, 2)[i]
function orderidx(idx, s::MOI.ExponentialCone)
    expmap.(idx)
end
function MOIU.loadconstraint!(instance::ECOSOptimizer, ci, f::MOI.VectorAffineFunction, s::MOI.AbstractVectorSet)
    A = sparse(f.outputindex, _varmap(f), f.coefficients)
    # sparse combines duplicates with + but does not remove zeros created so we call dropzeros!
    dropzeros!(A)
    colval = zeros(Int, length(A.rowval))
    for col in 1:A.n
        colval[A.colptr[col]:(A.colptr[col+1]-1)] = col
    end
    @assert !any(iszero.(colval))
    offset = constroffset(instance, ci)
    rows = constrrows(s)
    if s isa MOI.Zeros
        instance.cone.eqnrows[offset] = length(rows)
    else
        instance.cone.ineqnrows[offset] = length(rows)
    end
    i = offset + rows
    # The ECOS format is b - Ax ∈ cone
    # so minus=false for b and minus=true for A
    b, I, J, V = matrix(instance, s)
    b[i] = scalecoef(rows, orderval(f.constant, s), false, s)
    append!(I, offset + orderidx(A.rowval, s))
    append!(J, colval)
    append!(V, scalecoef(A.rowval, A.nzval, true, s))
end

function MOIU.allocatevariables!(instance::ECOSOptimizer, nvars::Integer)
    instance.cone = ConeData()
    VI.(1:nvars)
end

function MOIU.loadvariables!(instance::ECOSOptimizer, nvars::Integer)
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

MOIU.canallocate(::ECOSOptimizer, ::MOI.ObjectiveSense) = true
function MOIU.allocate!(instance::ECOSOptimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    instance.maxsense = sense == MOI.MaxSense
end
MOIU.canallocate(::ECOSOptimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}) = true
function MOIU.allocate!(::ECOSOptimizer, ::MOI.ObjectiveFunction, ::MOI.ScalarAffineFunction) end

MOIU.canload(::ECOSOptimizer, ::MOI.ObjectiveSense) = true
function MOIU.load!(::ECOSOptimizer, ::MOI.ObjectiveSense, ::MOI.OptimizationSense) end
MOIU.canload(::ECOSOptimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}) = true
function MOIU.load!(instance::ECOSOptimizer, ::MOI.ObjectiveFunction, f::MOI.ScalarAffineFunction)
    c0 = full(sparsevec(_varmap(f), f.coefficients, instance.data.n))
    instance.data.objconstant = f.constant
    instance.data.c = instance.maxsense ? -c0 : c0
end

function MOI.optimize!(instance::ECOSOptimizer)
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
    ecos_prob_ptr = ECOS.setup(n, m, cone.f, cone.l, length(cone.qa), cone.qa, cone.ep, G, A, c, h, b)
    ret_val = ECOS.solve(ecos_prob_ptr)
    ecos_prob = unsafe_wrap(Array, ecos_prob_ptr, 1)[1]
    primal    = unsafe_wrap(Array, ecos_prob.x, n)[:]
    dual_eq   = unsafe_wrap(Array, ecos_prob.y, cone.f)[:]
    dual_ineq = unsafe_wrap(Array, ecos_prob.z, m)[:]
    slack     = unsafe_wrap(Array, ecos_prob.s, m)[:]
    ECOS.cleanup(ecos_prob_ptr, 0)
    objval = (instance.maxsense ? -1 : 1) * dot(c, primal) + objconstant
    instance.sol = Solution(ret_val, primal, dual_eq, dual_ineq, slack, objval)
end

# Implements getter for result value and statuses
MOI.canget(instance::ECOSOptimizer, ::MOI.TerminationStatus) = true
function MOI.get(instance::ECOSOptimizer, ::MOI.TerminationStatus)
    flag = instance.sol.ret_val
    if flag == ECOS.ECOS_OPTIMAL
        MOI.Success
    elseif flag == ECOS.ECOS_PINF
        MOI.Success
    elseif flag == ECOS.ECOS_DINF  # Dual infeasible = primal unbounded, probably
        MOI.Success
    elseif flag == ECOS.ECOS_MAXIT
        MOI.IterationLimit
    elseif flag == ECOS.ECOS_OPTIMAL + ECOS.ECOS_INACC_OFFSET
        m.solve_stat = MOI.AlmostSuccess
    else
        m.solve_stat = MOI.OtherError
    end
end

MOI.canget(instance::ECOSOptimizer, ::MOI.ObjectiveValue) = true
MOI.get(instance::ECOSOptimizer, ::MOI.ObjectiveValue) = instance.sol.objval

function MOI.canget(instance::ECOSOptimizer, ::MOI.PrimalStatus)
    instance.sol.ret_val != ECOS.ECOS_PINF
end
function MOI.get(instance::ECOSOptimizer, ::MOI.PrimalStatus)
    flag = instance.sol.ret_val
    if flag == ECOS.ECOS_OPTIMAL
        MOI.FeasiblePoint
    elseif flag == ECOS.ECOS_PINF
        MOI.InfeasiblePoint
    elseif flag == ECOS.ECOS_DINF  # Dual infeasible = primal unbounded, probably
        MOI.InfeasibilityCertificate
    elseif flag == ECOS.ECOS_MAXIT
        MOI.UnknownResultStatus
    elseif flag == ECOS.ECOS_OPTIMAL + ECOS.ECOS_INACC_OFFSET
        m.solve_stat = MOI.NearlyFeasiblePoint
    else
        m.solve_stat = MOI.OtherResultStatus
    end
end
# Swapping indices 2 <-> 3 is an involution (it is its own inverse)
const reorderval = orderval
function MOI.canget(instance::ECOSOptimizer, ::Union{MOI.VariablePrimal, MOI.ConstraintPrimal}, ::Type{<:MOI.Index})
    instance.sol.ret_val != ECOS.ECOS_PINF
end
function MOI.get(instance::ECOSOptimizer, ::MOI.VariablePrimal, vi::VI)
    instance.sol.primal[vi.value]
end
MOI.get(instance::ECOSOptimizer, a::MOI.VariablePrimal, vi::Vector{VI}) = MOI.get.(instance, a, vi)
# setconstant: Retrieve set constant stored in `ConeData` during `copy!`
setconstant(instance::ECOSOptimizer, offset, ::CI{<:MOI.AbstractFunction, <:MOI.EqualTo}) = instance.cone.eqsetconstant[offset]
setconstant(instance::ECOSOptimizer, offset, ::CI) = instance.cone.ineqsetconstant[offset]
_unshift(instance::ECOSOptimizer, offset, value, ::CI) = value
_unshift(instance::ECOSOptimizer, offset, value, ci::CI{<:MOI.AbstractScalarFunction, <:MOI.AbstractScalarSet}) = value + setconstant(instance, offset, ci)
function MOI.get(instance::ECOSOptimizer, ::MOI.ConstraintPrimal, ci::CI{<:MOI.AbstractFunction, MOI.Zeros})
    rows = constrrows(instance, ci)
    zeros(length(rows))
end
function MOI.get(instance::ECOSOptimizer, ::MOI.ConstraintPrimal, ci::CI{<:MOI.AbstractFunction, <:MOI.EqualTo})
    offset = constroffset(instance, ci)
    setconstant(instance, offset, ci)
end
function MOI.get(instance::ECOSOptimizer, ::MOI.ConstraintPrimal, ci::CI{<:MOI.AbstractFunction, S}) where S <: MOI.AbstractSet
    offset = constroffset(instance, ci)
    rows = constrrows(instance, ci)
    _unshift(instance, offset, scalecoef(rows, reorderval(instance.sol.slack[offset + rows], S), false, S), ci)
end

function MOI.canget(instance::ECOSOptimizer, ::MOI.DualStatus)
    instance.sol.ret_val != ECOS.ECOS_DINF
end
function MOI.get(instance::ECOSOptimizer, ::MOI.DualStatus)
    flag = instance.sol.ret_val
    if flag == ECOS.ECOS_OPTIMAL
        MOI.FeasiblePoint
    elseif flag == ECOS.ECOS_PINF
        MOI.InfeasibilityCertificate
    elseif flag == ECOS.ECOS_DINF  # Dual infeasible = primal unbounded, probably
        MOI.InfeasiblePoint
    elseif flag == ECOS.ECOS_MAXIT
        MOI.UnknownResultStatus
    elseif flag == ECOS.ECOS_OPTIMAL + ECOS.ECOS_INACC_OFFSET
        m.solve_stat = MOI.NearlyFeasiblePoint
    else
        m.solve_stat = MOI.OtherResultStatus
    end
end
function MOI.canget(instance::ECOSOptimizer, ::MOI.ConstraintDual, ::Type{<:CI})
    instance.sol.ret_val != ECOS.ECOS_DINF
end
_dual(instance, ci::CI{<:MOI.AbstractFunction, <:ZeroCones}) = instance.sol.dual_eq
_dual(instance, ci::CI) = instance.sol.dual_ineq
function MOI.get(instance::ECOSOptimizer, ::MOI.ConstraintDual, ci::CI{<:MOI.AbstractFunction, S}) where S <: MOI.AbstractSet
    offset = constroffset(instance, ci)
    rows = constrrows(instance, ci)
    scalecoef(rows, reorderval(_dual(instance, ci)[offset + rows], S), false, S)
end

MOI.canget(instance::ECOSOptimizer, ::MOI.ResultCount) = true
MOI.get(instance::ECOSOptimizer, ::MOI.ResultCount) = 1
