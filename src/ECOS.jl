#############################################################################
# ECOS.jl
# Wrapper around the ECOS solver https://github.com/ifa-ethz/ecos
# See http://github.com/JuliaOpt/ECOS.jl
#############################################################################
# ECOS.jl
# Contains the wrapper itself
#############################################################################

module ECOS

using SparseArrays
using LinearAlgebra

# Try to load the binary dependency
if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("ECOS not properly installed. Please run Pkg.build(\"ECOS\")")
end

# ver  [not exported]
# Returns the version of ECOS in use
function ver()
    ver_ptr = ccall((:ECOS_ver, ECOS.ecos), Ptr{UInt8}, ())
    return unsafe_string(ver_ptr)
end

macro ecos_ccall(func, args...)
    args = map(esc,args)
    f = "ECOS_$(func)"
    quote
        @unix_only ret = ccall(($f,ECOS.ecos), $(args...))
        @windows_only ret = ccall(($f,ECOS.ecos), stdcall, $(args...))
        ret
    end
end

function __init__()
    vnum = VersionNumber(ver())
    if !(vnum.major == 2 && vnum.minor == 0)
        depsdir = realpath(joinpath(dirname(@__FILE__),"..","deps"))
        error("Current ECOS version installed is $(ver()), but we require version 2.0.*. On Linux and Windows, delete the contents of the `$depsdir` directory except for the files `build.jl` and `.gitignore`, then rerun Pkg.build(\"ECOS\"). On OS X, run `using Homebrew; Homebrew.update()` in Julia.")
    end
end

include("types.jl")  # All the types and constants defined in ecos.h

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


# setup  (direct interface)
# Provide ECOS with a problem in the form
# min c'x
# st  A*x = b
#     G*x in_K h,
#       or equivalently  h - G*x in K
#       or equivalently  G*x + s = h, s in K
# K is the product of R+ and multiple second-order cones, with the
# positive orthant first.
# Inputs:
#   n       Number of variables == length(x)
#   m       Number of inequality constraints == length(h)
#   p       Number of equality constraints == length(b)  (can be 0)
#   l       Dimension of the positive orthant,
#               i.e. the first l elements of s are >= 0
#   ncones  Number of second-order cones present in K
#   q       ncones-long vector of the length of the index sets of the SOCs
#               e.g. cone 1 is indices 4:6, cone 2 is indices 7:10
#                    ->  q = [3, 4]
#   e       Number of exponential cones present in problem
#   G       The G matrix in ECOSMatrix format.
#   A       The A matrix in ECOSMatrix format.
#           Can be all nothing if no equalities are present.
#   c       Objective coefficients, length(c) == n
#   h       RHS for inequality constraints, length(h) == m
#   b       RHS for equality constraints, length(b) == b (can be nothing)
# Returns a pointer to the ECOS pwork structure (Cpwork in ECOS.jl). See
# types.jl for more information.
# **NOTE**: ECOS retains references to the problem data passed in here.
# You *must* ensure that G, A, c, h, and b are not freed until after cleanup(), otherwise
# memory corruption will occur.
function setup(n::Int, m::Int, p::Int, l::Int, ncones::Int, q::Union{Vector{Int},Nothing}, e::Int,
               G::ECOSMatrix, A::Union{ECOSMatrix, Nothing},
               c::Vector{Float64}, h::Vector{Float64}, b::Union{Vector{Float64},Nothing}; kwargs...)
    # Convert to canonical forms
    q = (q == nothing) ? convert(Ptr{Clong}, C_NULL) : convert(Vector{Clong},q)
    Apr = (A == nothing) ? convert(Ptr{Cdouble}, C_NULL) : A.pr
    Ajc = (A == nothing) ? convert(Ptr{Cdouble}, C_NULL) : A.jc
    Air = (A == nothing) ? convert(Ptr{Cdouble}, C_NULL) : A.ir
    b = (b == nothing) ? convert(Ptr{Cdouble}, C_NULL) : b
    problem_ptr = ccall((:ECOS_setup, ECOS.ecos), Ptr{Cpwork},
        (Clong, Clong, Clong, Clong, Clong, Ptr{Clong}, Clong,
         Ptr{Cdouble}, Ptr{Clong}, Ptr{Clong},
         Ptr{Cdouble}, Ptr{Clong}, Ptr{Clong},
         Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
        n, m, p, l, ncones, q, e,
        G.pr, G.jc, G.ir,
        Apr, Ajc, Air,
        c, h, b)

    problem_ptr != C_NULL || error("ECOS failed to construct problem.")

    if !isempty(kwargs)
        problem = unsafe_load(problem_ptr)
        problem.stgs != C_NULL || error("ECOS returned a malformed settings struct.")
        settings = unsafe_load(problem.stgs)

        options = Dict{Symbol, Any}()
        for (k,v) in kwargs
            options[k] = v
        end

        new_settings = Csettings([
            setting in keys(options) ?
            convert(fieldtype(typeof(settings), setting), options[setting]) :
            getfield(settings, setting)
            for setting in fieldnames(typeof(settings))]...)

        unsafe_store!(problem.stgs, new_settings)
    end

    problem_ptr
end

# solve
# Solves the provided problem. Results are stored inside the structure,
# but currently there is no convenient interface-provided way to access
# this - use MathProgBase interface.
function solve(problem::Ptr{Cpwork})
    exitflag = ccall((:ECOS_solve, ECOS.ecos), Clong, (Ptr{Cpwork},), problem)
end

# cleanup
# Frees memory allocated by ECOS for the problem.
# The optional keepvars argument is number of variables to NOT free.
function cleanup(problem::Ptr{Cpwork}, keepvars::Int = 0)
    ccall((:ECOS_cleanup, ECOS.ecos), Cvoid, (Ptr{Cpwork}, Clong), problem, keepvars)
end

include("MPB_wrapper.jl")
include("MOI_wrapper.jl")

end # module
