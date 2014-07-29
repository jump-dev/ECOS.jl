#############################################################################
# ECOS.jl
# Wrapper around the ECOS solver https://github.com/ifa-ethz/ecos
# See http://github.com/JuliaOpt/ECOS.jl
#############################################################################
# ECOS.jl
# Contains the wrapper itself
#############################################################################

module ECOS

# Try to load the binary dependency
if isfile(joinpath(Pkg.dir("ECOS"),"deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("ECOS not properly installed. Please run Pkg.build(\"ECOS\")")
end

macro ecos_ccall(func, args...)
    f = "ECOS_$(func)"
    quote
        @unix_only ret = ccall(($f,ECOS.ecos), $(args...))
        @windows_only ret = ccall(($f,ECOS.ecos), stdcall, $(args...))
        ret
    end
end


include("ECOSSolverInterface.jl")  # MathProgBase interface
include("types.jl")  # All the types and constants defined in ecos.h

export setup, solve, cleanup

# parameters
# n is the number of variables,
# m is the number of inequality constraints (dimension 1 of the matrix G and the length of the vector h),
# p is the number of equality constraints (can be 0)
# l is the dimension of the positive orthant, i.e. in Gx+s=h, s in K, the first l elements of s are >=0
# ncones is the number of second-order cones present in K
# q is an array of integers of length ncones, where each element defines the dimension of the cone
# Gpr, Gjc, Gir are the the data, the column index, and the row index arrays, respectively,
# for the matrix G represented in column compressed storage (CCS) format (Google it if you need more
# information on this format, it is one of the standard sparse matrix representations)
# Apr, Ajc, Air is the CCS representation of the matrix A (can be all NULL if no equalities are present)
# c is an array of type pfloat of size n
# h is an array of type pfloat of size m
# b is an array of type pfloat of size p (can be NULL if no equalities are present)
# The setup function returns a struct of type pwork, which you need to define first
# This is the straightforward translation from the C API.
function setup(n::Int64, m::Int64, p::Int64, l::Int64, ncones::Int64, q::Vector{Int64},
        Gpr::Vector{Float64}, Gjc::Vector{Int64}, Gir::Vector{Int64}, Apr::Vector{Float64},
        Ajc::Vector{Int64}, Air::Vector{Int64}, c::Vector{Float64}, h::Vector{Float64},
        b::Vector{Float64})
    problem = ccall((:ECOS_setup, ECOS.ecos), Ptr{Cpwork}, (Clong, Clong, Clong, Clong, Clong, Ptr{Clong}, Ptr{Cdouble}, Ptr{Clong}, Ptr{Clong}, Ptr{Cdouble}, Ptr{Clong}, Ptr{Clong}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), n, m, p, l, ncones, q, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b)
end

VecOrMatOrSparseOrNothing = Union(Vector, Matrix, SparseMatrixCSC, Nothing)
ArrayFloat64OrNothing = Union(Vector{Float64}, Nothing)

function setup(;n::Int64=nothing, m::Int64=nothing, p::Int64=0, l::Int64=0,
        ncones::Int64=0, q::Vector{Int64}=Int64[], G::VecOrMatOrSparseOrNothing=nothing,
        A::VecOrMatOrSparseOrNothing=nothing, c::Vector{Float64}=nothing,
        h::ArrayFloat64OrNothing=nothing, b::ArrayFloat64OrNothing=nothing)
    
    if length(q) == 0
        q = convert(Ptr{Int64}, C_NULL)
    end

    if A == nothing
        Apr = convert(Ptr{Float64}, C_NULL)
        Ajc = convert(Ptr{Int64}, C_NULL)
        Air = convert(Ptr{Int64}, C_NULL)
    else
        sparseA = sparse(A)
        # Hack to make it a float, find a better way
        Apr = sparseA.nzval * 1.0
        # -1 since the C language is 0 indexed
        Ajc = sparseA.colptr - 1
        Air = sparseA.rowval - 1
    end

    if G == nothing
        Gpr = convert(Ptr{Float64}, C_NULL)
        Gjc = convert(Ptr{Int64}, C_NULL)
        Gir = convert(Ptr{Int64}, C_NULL)
    else
        sparseG = sparse(G)
        # Hack to make it a float, find a better way
        Gpr = sparseG.nzval * 1.0
        # -1 since the C language is 0 indexed
        Gjc = sparseG.colptr - 1
        Gir = sparseG.rowval - 1
    end

    if b == nothing
        b = convert(Ptr{Float64}, C_NULL)
    end

    if h == nothing
        h = convert(Ptr{Float64}, C_NULL)
    end

    ccall((:ECOS_setup, ECOS.ecos), Ptr{Cpwork}, (Clong, Clong, Clong, Clong, Clong,
            Ptr{Clong}, Ptr{Cdouble}, Ptr{Clong}, Ptr{Clong}, Ptr{Cdouble}, Ptr{Clong},
            Ptr{Clong}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), n, m, p, l, ncones,
            q, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b)
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
    ccall((:ECOS_cleanup, ECOS.ecos), Void, (Ptr{Cpwork}, Clong), problem, keepvars)
end

# ver  [not exported]
# Returns the version of ECOS in use
function ver()
    ver_ptr = ccall((:ECOS_ver, ECOS.ecos), Ptr{Uint8}, ())
    return bytestring(pointer_to_array(ver_ptr, 7))
end

end # module
