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

# @fieldtype
# for 0.3 compatibility
macro fieldtype(instance, field)
    if VERSION < v"0.4-"
        return :(fieldtype($instance, $field))
    else
        return :(fieldtype(typeof($instance), $field))
    end
end


include("ECOSSolverInterface.jl")  # MathProgBase interface
include("types.jl")  # All the types and constants defined in ecos.h

export setup, solve, cleanup

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
#   Gpr, Gjc, Gir
#           Non-zeros, column indices, and the row index arrays, respectively,
#           for the matrix G represented in column compressed storage (CCS) format
#           which is equivalent to Julia's SparseMatrixCSC
#   Apr, Ajc, Air
#           Equivalent to above for the matrix A.
#           Can be all nothing if no equalities are present.
#   c       Objective coefficients, length(c) == n
#   h       RHS for inequality constraints, length(h) == m
#   b       RHS for equality constraints, length(b) == b (can be nothing)
# Returns a pointer to the ECOS pwork structure (Cpwork in ECOS.jl). See
# types.jl for more information.
function setup(n::Int, m::Int, p::Int, l::Int, ncones::Int, q::Union(Vector{Int},Nothing),
                Gpr::Vector{Float64}, Gjc::Vector{Int}, Gir::Vector{Int},
                Apr::Union(Vector{Float64},Nothing), Ajc::Union(Vector{Int},Nothing), Air::Union(Vector{Int},Nothing),
                c::Vector{Float64}, h::Vector{Float64}, b::Union(Vector{Float64},Nothing); kwargs...)
    # Convert to canonical forms
    q = (q == nothing) ? convert(Ptr{Clong}, C_NULL) : convert(Vector{Clong},q)
    Apr = (Apr == nothing) ? convert(Ptr{Cdouble}, C_NULL) : Apr
    Ajc = (Ajc == nothing) ? convert(Ptr{Cdouble}, C_NULL) : convert(Vector{Clong},Ajc)
    Air = (Air == nothing) ? convert(Ptr{Cdouble}, C_NULL) : convert(Vector{Clong},Air)
    b = (b == nothing) ? convert(Ptr{Cdouble}, C_NULL) : b
    problem_ptr = ccall((:ECOS_setup, ECOS.ecos), Ptr{Cpwork},
        (Clong, Clong, Clong, Clong, Clong, Ptr{Clong},
         Ptr{Cdouble}, Ptr{Clong}, Ptr{Clong},
         Ptr{Cdouble}, Ptr{Clong}, Ptr{Clong},
         Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
        n, m, p, l, ncones, q,
        Gpr, convert(Vector{Clong},Gjc), convert(Vector{Clong},Gir),
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
            convert(@fieldtype(settings, setting), options[setting]) :
            getfield(settings, setting)
            for setting in names(settings)]...)

        unsafe_store!(problem.stgs, new_settings)
    end

    problem_ptr
end

# setup  (more general interface)
# A more tolerant version of the above method that doesn't require
# user to fiddle with internals of the sparse matrix format
# User can pass nothing as argument for A, b, and q
function setup(n, m, p, l, ncones, q, G, A, c, h, b; options...)
    if A == nothing
        if b != nothing
            @assert length(b) == 0
            b = nothing
        end
        @assert p == 0
        Apr = nothing
        Ajc = nothing
        Air = nothing
    else
        numrow, numcol = size(A)
        @assert numcol == n
        @assert numrow == length(b)
        sparseA = sparse(A)
        Apr = convert(Vector{Float64}, sparseA.nzval)
        Ajc = sparseA.colptr - 1  # C is 0-based
        Air = sparseA.rowval - 1
    end

    sparseG = sparse(G)
    Gpr = convert(Vector{Float64}, sparseG.nzval)
    Gjc = sparseG.colptr - 1  # C is 0-based
    Gir = sparseG.rowval - 1

    setup(  n, m, p, l, ncones, q, 
            Gpr, Gjc, Gir,
            Apr, Ajc, Air,
            c, h, b; options...)
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
