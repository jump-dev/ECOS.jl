#############################################################################
# ECOS.jl
# Wrapper around the ECOS solver https://github.com/ifa-ethz/ecos
# See http://github.com/jump-dev/ECOS.jl
#############################################################################
# ECOS.jl
# Contains the wrapper itself
#############################################################################

module ECOS

using SparseArrays
using LinearAlgebra

# Try to load the binary dependency
if VERSION < v"1.3"
    _DEPS_FILE = joinpath(dirname(@__FILE__),"..","deps","deps.jl")
    if isfile(_DEPS_FILE)
        include(_DEPS_FILE)
    else
        error("ECOS not properly installed. Please run Pkg.build(\"ECOS\")")
    end
else
    import ECOS_jll
    const ecos = ECOS_jll.libecos
end

using CEnum
include("gen/ctypes.jl")
include("gen/libecos_common.jl")
include("gen/libecos_api.jl")

# ver  [not exported]
# Returns the version of ECOS in use
#
# Note: Avoid calling this function at the top-level as doing so breaks creating system
# images which include this package
# https://github.com/jump-dev/ECOS.jl/issues/106
ver() = unsafe_string(ECOS_ver())

function __init__()
    libecos_version = VersionNumber(ver())
    if libecos_version < v"2.0.5"
        error("Current ECOS version installed is $(ver()), but we require minimum version of 2.0.8")
    end
    return
end

# A julia container for matrices passed to ECOS.
# This makes it easy for the Julia GC to track these arrays.
# The format matches ECOS's spmat: 0-based compressed sparse column (CSC)
struct ECOSMatrix
    pr::Vector{Cdouble} # nzval
    jc::Vector{Clong}   # colptr
    ir::Vector{Clong}   # rowval
end

function ECOSMatrix(mat::AbstractMatrix)
    sparsemat = sparse(mat)
    return ECOSMatrix(sparsemat.nzval, sparsemat.colptr .- 1, sparsemat.rowval .- 1)
end

mutable struct _ECOS_setup_struct
    ptr::Ptr{pwork}
    q::Union{Vector{Cint},Nothing}
    G::ECOSMatrix
    A::Union{ECOSMatrix,Nothing}
    c::Vector{Float64}
    h::Vector{Float64}
    b::Union{Vector{Float64},Nothing}
end

Base.cconvert(::Type{Ptr{pwork}}, x::_ECOS_setup_struct) = x
Base.unsafe_convert(::Type{Ptr{pwork}}, x::_ECOS_setup_struct) = x.ptr

"""
    setup(args...)

Return a pointer to the ECOS pwork structure (pwork in ECOS.jl) containing a
problem in the form:

    min c'x
    st  A*x = b
        G*x in_K h,
        or equivalently  h - G*x in K
        or equivalently  G*x + s = h, s in K

where K is the product of R+ and multiple second-order cones, with the positive
orthant first.

## Inputs

  n       Number of variables == length(x)
  m       Number of inequality constraints == length(h)
  p       Number of equality constraints == length(b)  (can be 0)
  l       Dimension of the positive orthant,
              i.e. the first l elements of s are >= 0
  ncones  Number of second-order cones present in K
  q       ncones-long vector of the length of the index sets of the SOCs
              e.g. cone 1 is indices 4:6, cone 2 is indices 7:10
                   ->  q = [3, 4]
  e       Number of exponential cones present in problem
  G       The G matrix in ECOSMatrix format.
  A       The A matrix in ECOSMatrix format.
          Can be all nothing if no equalities are present.
  c       Objective coefficients, length(c) == n
  h       RHS for inequality constraints, length(h) == m
  b       RHS for equality constraints, length(b) == b (can be nothing)

## Note

If you call this method, you do _NOT_ need to call ECOS_cleanup.
"""
function setup(
    n::Int,
    m::Int,
    p::Int,
    l::Int,
    ncones::Int,
    q::Union{Vector{Int},Nothing},
    e::Int,
    G::ECOSMatrix,
    A::Union{ECOSMatrix,Nothing},
    c::Vector{Float64},
    h::Vector{Float64},
    b::Union{Vector{Float64},Nothing};
    kwargs...,
)
    problem_ptr = ECOS_setup(
        n,
        m,
        p,
        l,
        ncones,
        q === nothing ? C_NULL : q,
        e,
        G.pr,
        G.jc,
        G.ir,
        A === nothing ? C_NULL : A.pr,
        A === nothing ? C_NULL : A.jc,
        A === nothing ? C_NULL : A.ir,
        c,
        h,
        b === nothing ? C_NULL : b,
    )
    if problem_ptr == C_NULL
        error("ECOS failed to construct problem.")
    end
    if !isempty(kwargs)
        problem = unsafe_load(problem_ptr)
        if problem.stgs == C_NULL
            error("ECOS returned a malformed settings struct.")
        end
        old_settings = unsafe_load(problem.stgs)
        options = Dict{Symbol,Any}(k => v for (k, v) in kwargs)
        new_settings = settings(
            [
                convert(
                    fieldtype(settings, key),
                    get(options, key, getfield(old_settings, key)),
                )
                for key in fieldnames(settings)
            ]...
        )
        unsafe_store!(problem.stgs, new_settings)
    end
    model = _ECOS_setup_struct(
        problem_ptr,
        q,
        G,
        A,
        c,
        h,
        b,
    )
    finalizer(model) do m
        ECOS_cleanup(m, 0)
    end
    return model
end

include("MPB_wrapper.jl")
include("MOI_wrapper.jl")

end # module
