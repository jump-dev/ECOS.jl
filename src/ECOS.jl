#############################################################################
# ECOS.jl
# Wrapper around the ECOS solver https://github.com/ifa-ethz/ecos
# See http://github.com/JuliaOpt/ECOS.jl
#############################################################################
# ECOS.jl
# Contains the wrapper itself
#############################################################################

module ECOS

if isfile(joinpath(Pkg.dir("ECOS"),"deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("ECOS not properly installed. Please run Pkg.build(\"ECOS\")")
end

include("ECOSSolverInterface.jl")


macro ecos_ccall(func, args...)
    f = "ECOS_$(func)"
    quote
        @unix_only ret = ccall(($f,ECOS.ecos), $(args...))
        @windows_only ret = ccall(($f,ECOS.ecos), stdcall, $(args...))
        ret
    end
end

immutable Clpcone
    p::Clong
    w::Ptr{Cdouble}
    v::Ptr{Cdouble}
    kkt_idx::Ptr{Clong}
end

immutable Csocone
    p::Clong
    skbar::Ptr{Cdouble}
    zkbar::Ptr{Cdouble}
    a::Clong
    d1::Clong
    w::Cdouble
    eta::Cdouble
    eta_square::Cdouble
    q::Ptr{Cdouble}
    Didx::Ptr{Clong}
    u0::Cdouble
    u1::Cdouble
    v1::Cdouble
end

immutable Ccone
    lpcone::Ptr{Clpcone}
    socone::Ptr{Csocone}
    nsoc::Clong
end

immutable Cspmat
    jc::Ptr{Clong}
    ir::Ptr{Clong}
    pr::Ptr{Cdouble}
    n::Clong
    m::Clong
    nnz::Clong
end

immutable Ckkt
    PKPt::Ptr{Cspmat}
    L::Ptr{Cspmat}
    D::Ptr{Cdouble}
    work1::Ptr{Cdouble}
    work2::Ptr{Cdouble}
    work3::Ptr{Cdouble}
    work4::Ptr{Cdouble}
    work5::Ptr{Cdouble}
    work6::Ptr{Cdouble}
    RHS1::Ptr{Cdouble}
    RHS2::Ptr{Cdouble}
    dx1::Ptr{Cdouble}
    dx2::Ptr{Cdouble}
    dy1::Ptr{Cdouble}
    dy2::Ptr{Cdouble}
    dz1::Ptr{Cdouble}
    dz2::Ptr{Cdouble}
    P::Ptr{Clong}
    Pinv::Ptr{Clong}
    PK::Ptr{Clong}
    Parent::Ptr{Clong}
    Sign::Ptr{Clong}
    Pattern::Ptr{Clong}
    Flag::Ptr{Clong}
    Lnz::Ptr{Clong}
    delta::Cdouble
end

immutable Cstats
    pcost::Cdouble
    dcost::Cdouble
    pres::Cdouble
    dres::Cdouble
    pinf::Cdouble
    dinf::Cdouble
    pinfres::Cdouble
    dinfres::Cdouble
    gap::Cdouble
    relgap::Cdouble
    sigma::Cdouble
    mu::Cdouble
    step::Cdouble
    step_aff::Cdouble
    kapovert::Cdouble
    iter::Clong
    nitref1::Clong
    nitref2::Clong
    nitref3::Clong
    tsetup::Cdouble
    tsolve::Cdouble
    tfactor::Cdouble
    tkktsolve::Cdouble
    torder::Cdouble
    tkktcreate::Cdouble
    ttranspose::Cdouble
    tperm::Cdouble
    tfactor_t1::Cdouble
    tfactor_t2::Cdouble
end

immutable Csettings
    gamma::Cdouble
    delta::Cdouble
    eps::Cdouble
    feastol::Cdouble
    abstol::Cdouble
    reltol::Cdouble
    feastol_inacc::Cdouble
    abstol_inacc::Cdouble
    reltol_inacc::Cdouble
    nitref::Clong
    maxit::Clong
    verbose::Clong
end

immutable Cpwork
    # Dimensions
    n::Clong
    m::Clong
    p::Clong
    D::Clong

    # Variables
    x::Ptr{Cdouble}
    y::Ptr{Cdouble}
    z::Ptr{Cdouble}
    s::Ptr{Cdouble}
    lambda::Ptr{Cdouble}
    kap::Cdouble
    tau::Cdouble

    # Best iterates seen so far
    best_x::Ptr{Cdouble}
    best_y::Ptr{Cdouble}
    best_z::Ptr{Cdouble}
    best_s::Ptr{Cdouble}
    best_kap::Cdouble
    best_tau::Cdouble
    best_cx::Cdouble
    best_by::Cdouble
    best_hz::Cdouble
    best_info::Ptr{Cstats}

    # Temporary variables
    dsaff::Ptr{Cdouble}
    dzaff::Ptr{Cdouble}
    W_times_dzaff::Ptr{Cdouble}
    dsaff_by_W::Ptr{Cdouble}
    saff::Ptr{Cdouble}
    zaff::Ptr{Cdouble}

    # Cone
    C::Ptr{Ccone}
    A::Ptr{Cspmat}
    G::Ptr{Cspmat}
    c::Ptr{Cdouble}
    b::Ptr{Cdouble}
    h::Ptr{Cdouble}

    # equilibration vector
    xequil::Ptr{Cdouble}
    Aequil::Ptr{Cdouble}
    Gequil::Ptr{Cdouble}

    # scalings of problem data
    resx0::Cdouble
    resy0::Cdouble
    resz0::Cdouble

    # residuals
    rx::Ptr{Cdouble}
    ry::Ptr{Cdouble}
    rz::Ptr{Cdouble}
    rt::Cdouble
    hresx::Cdouble
    hresy::Cdouble
    hresz::Cdouble

    # temporary storage
    cx::Cdouble
    by::Cdouble
    hz::Cdouble
    sz::Cdouble

    # KKT System
    KKT::Ptr{Ckkt}

    # info struct
    info::Ptr{Cstats}

    # settings struct
    stgs::Ptr{Csettings}
end


# Solve status flags
# From: https://github.com/ifa-ethz/ecos/blob/master/include/ecos.h
const ECOS_OPTIMAL      = 0   # Problem solved to optimality
const ECOS_PINF         = 1   # Found certificate of primal infeasibility
const ECOS_DINF         = 2   # Found certificate of dual infeasibility
const ECOS_INACC_OFFSET = 10  # Offset exitflag at inaccurate results
const ECOS_MAXIT        = -1  # Maximum number of iterations reached
const ECOS_NUMERICS     = -2  # Search direction unreliable
const ECOS_OUTCONE      = -3  # s or z got outside the cone, numerics?
const ECOS_FATAL        = -7  # Unknown problem in solver


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

function solve(problem::Ptr{Cpwork})
    exitflag = ccall((:ECOS_solve, ECOS.ecos), Clong, (Ptr{Cpwork},), problem)
end

function cleanup(problem::Ptr{Cpwork}, keepvars::Clong)
    ccall((:ECOS_cleanup, ECOS.ecos), Void, (Ptr{Cpwork}, Clong), problem, keepvars)
end

end # module
