module ECOS

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

# ECOS types: 
# pfloat -> Cdouble
# idxint -> SuiteSparse_long -> Clong

# parameters
# n is the number of variables,
# m is the number of inequality constraints (dimension 1 of the matrix G and the length of the vector h),
# p is the number of equality constraints (can be 0)
# l is the dimension of the positive orthant, i.e. in Gx+s=h, s in K, the first l elements of s are >=0
# ncones is the number of second-order cones present in K
# q is an array of integers of length ncones, where each element defines the dimension of the cone
# Gpr, Gjc, Gir are the the data, the column index, and the row index arrays, respectively, for the matrix G represented in column compressed storage (CCS) format (Google it if you need more information on this format, it is one of the standard sparse matrix representations)
# Apr, Ajc, Air is the CCS representation of the matrix A (can be all NULL if no equalities are present)
# c is an array of type pfloat of size n
# h is an array of type pfloat of size m

# b is an array of type pfloat of size p (can be NULL if no equalities are present) The setup function returns a struct of type pwork, which you need to define first


# This is the straightforward translation from the C API.
# TODO: allow the matrix to be a Julia SparseMatrixCSC instead
function setup(n::Int64, m::Int64, p::Int64, l::Int64, ncones::Int64, q::Array{Int64}, Gpr::Array{Float64}, Gjc::Array{Int64}, Gir::Array{Int64}, Apr::Array{Float64}, Ajc::Array{Int64}, Air::Array{Int64}, c::Array{Float64}, h::Array{Float64}, b::Array{Float64})

problem = ccall((:ECOS_setup, ECOS.ecos), Ptr{Void}, (Clong, Clong, Clong, Clong, Clong, Ptr{Clong}, Ptr{Cdouble}, Ptr{Clong}, Ptr{Clong}, Ptr{Cdouble}, Ptr{Clong}, Ptr{Clong}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), n, m, p, l, ncones, q, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b)

end

# ECOS_setup(idxint n, idxint m, idxint p, idxint l, idxint ncones, idxint* q,
#                    pfloat* Gpr, idxint* Gjc, idxint* Gir,
#                    pfloat* Apr, idxint* Ajc, idxint* Air,
#                    pfloat* c, pfloat* h, pfloat* b);

end # module
