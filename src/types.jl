#############################################################################
# ECOS.jl
# Wrapper around the ECOS solver https://github.com/ifa-ethz/ecos
# See http://github.com/JuliaOpt/ECOS.jl
#############################################################################
# types.jl
# Julia implementations of the types defined in ecos.h
# Types are mapped as follows
# pfloat -> Cdouble
# idxint -> SuiteSparse_long -> Clong
#############################################################################

# Solve status flags
# From: https://github.com/ifa-ethz/ecos/blob/master/include/ecos.h
const ECOS_OPTIMAL      =  0  # Problem solved to optimality
const ECOS_PINF         =  1  # Found certificate of primal infeasibility
const ECOS_DINF         =  2  # Found certificate of dual infeasibility
const ECOS_INACC_OFFSET = 10  # Offset exitflag at inaccurate results
const ECOS_MAXIT        = -1  # Maximum number of iterations reached
const ECOS_NUMERICS     = -2  # Search direction unreliable
const ECOS_OUTCONE      = -3  # s or z got outside the cone, numerics?
const ECOS_SIGINT       = -4  # solver interrupted by a signal/ctrl-c
const ECOS_FATAL        = -7  # Unknown problem in solver


struct Clpcone
    p::Clong
    w::Ptr{Cdouble}
    v::Ptr{Cdouble}
    kkt_idx::Ptr{Clong}
end

struct Csocone
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

struct Ccone
    lpcone::Ptr{Clpcone}
    socone::Ptr{Csocone}
    nsoc::Clong
end

struct Cspmat
    jc::Ptr{Clong}
    ir::Ptr{Clong}
    pr::Ptr{Cdouble}
    n::Clong
    m::Clong
    nnz::Clong
end

struct Ckkt
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

struct Cstats
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

struct Csettings
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
    max_bk_iter::Clong
    bk_scale::Cdouble
    centrality::Cdouble
end

if VersionNumber(ver()) >= v"2.0.5"
    @eval struct Cpwork
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

        # indices that map entries of A and G to the KKT matrix
        AtoK::Ptr{Clong}
        GtoK::Ptr{Clong}

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

        # norm iterates
        nx::Cdouble
        ny::Cdouble
        nz::Cdouble
        ns::Cdouble

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
else
    @eval struct Cpwork
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

        # norm iterates
        nx::Cdouble
        ny::Cdouble
        nz::Cdouble
        ns::Cdouble

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
end
