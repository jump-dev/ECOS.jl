# Automatically generated using Clang.jl


const CONEMODE = 0
const INSIDE_CONE = 0
const OUTSIDE_CONE = 1
const idxint = Clong
const pfloat = Cdouble

struct lpcone
    p::idxint
    w::Ptr{pfloat}
    v::Ptr{pfloat}
    kkt_idx::Ptr{idxint}
end

struct socone
    p::idxint
    skbar::Ptr{pfloat}
    zkbar::Ptr{pfloat}
    a::pfloat
    d1::pfloat
    w::pfloat
    eta::pfloat
    eta_square::pfloat
    q::Ptr{pfloat}
    Didx::Ptr{idxint}
    u0::pfloat
    u1::pfloat
    v1::pfloat
end

struct expcone
    colstart::NTuple{3, idxint}
    v::NTuple{6, pfloat}
    g::NTuple{3, pfloat}
end

struct cone
    lpc::Ptr{lpcone}
    soc::Ptr{socone}
    nsoc::idxint
    expc::Ptr{expcone}
    nexc::idxint
    fexv::idxint
end

# Skipping MacroDefinition: init_ctrlc ( )
# Skipping MacroDefinition: remove_ctrlc ( )
# Skipping MacroDefinition: check_ctrlc ( ) 0

const ECOS_VERSION = "2.0.8"
const MAXIT = 100
const FEASTOL = 1.0e-8
const ABSTOL = 1.0e-8
const RELTOL = 1.0e-8
const FTOL_INACC = 0.0001
const ATOL_INACC = 5.0e-5
const RTOL_INACC = 5.0e-5
const GAMMA = 0.99
const STATICREG = 1
const DELTASTAT = 7.0e-8
const DELTA = 2.0e-7
const EPS = 1.0e-13
const VERBOSE = 1
const NITREF = 9
const IRERRFACT = 6
const LINSYSACC = 1.0e-14
const SIGMAMIN = 0.0001
const SIGMAMAX = 1.0
const STEPMIN = 1.0e-6
const STEPMAX = 0.999
const SAFEGUARD = 500
const MAX_BK = 90
const BK_SCALE = 0.8
const MIN_DISTANCE = 0.1
const CENTRALITY = 1
const EQUILIBRATE = 1
const EQUIL_ITERS = 3
const ECOS_OPTIMAL = 0
const ECOS_PINF = 1
const ECOS_DINF = 2
const ECOS_INACC_OFFSET = 10
const ECOS_MAXIT = -1
const ECOS_NUMERICS = -2
const ECOS_OUTCONE = -3
const ECOS_SIGINT = -4
const ECOS_FATAL = -7

# Skipping MacroDefinition: MAX ( X , Y ) ( ( X ) < ( Y ) ? ( Y ) : ( X ) )
# Skipping MacroDefinition: SAFEDIV_POS ( X , Y ) ( ( Y ) < EPS ? ( ( X ) / EPS ) : ( X ) / ( Y ) )

struct settings
    gamma::pfloat
    delta::pfloat
    eps::pfloat
    feastol::pfloat
    abstol::pfloat
    reltol::pfloat
    feastol_inacc::pfloat
    abstol_inacc::pfloat
    reltol_inacc::pfloat
    nitref::idxint
    maxit::idxint
    verbose::idxint
    max_bk_iter::idxint
    bk_scale::pfloat
    centrality::pfloat
end

struct stats
    pcost::pfloat
    dcost::pfloat
    pres::pfloat
    dres::pfloat
    pinf::pfloat
    dinf::pfloat
    pinfres::pfloat
    dinfres::pfloat
    gap::pfloat
    relgap::pfloat
    sigma::pfloat
    mu::pfloat
    step::pfloat
    step_aff::pfloat
    kapovert::pfloat
    iter::idxint
    nitref1::idxint
    nitref2::idxint
    nitref3::idxint
    tsetup::pfloat
    tsolve::pfloat
    pob::idxint
    cb::idxint
    cob::idxint
    pb::idxint
    db::idxint
    affBack::idxint
    cmbBack::idxint
    centrality::pfloat
end

struct spmat
    jc::Ptr{idxint}
    ir::Ptr{idxint}
    pr::Ptr{pfloat}
    n::idxint
    m::idxint
    nnz::idxint
end

struct kkt
    PKPt::Ptr{spmat}
    L::Ptr{spmat}
    D::Ptr{pfloat}
    work1::Ptr{pfloat}
    work2::Ptr{pfloat}
    work3::Ptr{pfloat}
    work4::Ptr{pfloat}
    work5::Ptr{pfloat}
    work6::Ptr{pfloat}
    RHS1::Ptr{pfloat}
    RHS2::Ptr{pfloat}
    dx1::Ptr{pfloat}
    dx2::Ptr{pfloat}
    dy1::Ptr{pfloat}
    dy2::Ptr{pfloat}
    dz1::Ptr{pfloat}
    dz2::Ptr{pfloat}
    P::Ptr{idxint}
    Pinv::Ptr{idxint}
    PK::Ptr{idxint}
    Parent::Ptr{idxint}
    Sign::Ptr{idxint}
    Pattern::Ptr{idxint}
    Flag::Ptr{idxint}
    Lnz::Ptr{idxint}
    delta::pfloat
end

struct pwork
    n::idxint
    m::idxint
    p::idxint
    D::idxint
    x::Ptr{pfloat}
    y::Ptr{pfloat}
    z::Ptr{pfloat}
    s::Ptr{pfloat}
    lambda::Ptr{pfloat}
    kap::pfloat
    tau::pfloat
    best_x::Ptr{pfloat}
    best_y::Ptr{pfloat}
    best_z::Ptr{pfloat}
    best_s::Ptr{pfloat}
    best_kap::pfloat
    best_tau::pfloat
    best_cx::pfloat
    best_by::pfloat
    best_hz::pfloat
    best_info::Ptr{stats}
    dsaff::Ptr{pfloat}
    dzaff::Ptr{pfloat}
    W_times_dzaff::Ptr{pfloat}
    dsaff_by_W::Ptr{pfloat}
    saff::Ptr{pfloat}
    zaff::Ptr{pfloat}
    C::Ptr{cone}
    A::Ptr{spmat}
    G::Ptr{spmat}
    c::Ptr{pfloat}
    b::Ptr{pfloat}
    h::Ptr{pfloat}
    AtoK::Ptr{idxint}
    GtoK::Ptr{idxint}
    xequil::Ptr{pfloat}
    Aequil::Ptr{pfloat}
    Gequil::Ptr{pfloat}
    resx0::pfloat
    resy0::pfloat
    resz0::pfloat
    rx::Ptr{pfloat}
    ry::Ptr{pfloat}
    rz::Ptr{pfloat}
    rt::pfloat
    hresx::pfloat
    hresy::pfloat
    hresz::pfloat
    nx::pfloat
    ny::pfloat
    nz::pfloat
    ns::pfloat
    cx::pfloat
    by::pfloat
    hz::pfloat
    sz::pfloat
    KKT::Ptr{kkt}
    info::Ptr{stats}
    stgs::Ptr{settings}
end

const MI_PRINTLEVEL = 1
const MI_ABS_EPS = 1.0e-6
const MI_REL_EPS = 0.001
const MI_MAXITER = 1000
const MI_INT_TOL = FTOL_INACC
const MI_SOLVED_NON_BRANCHABLE = 3
const MI_SOLVED_BRANCHABLE = 2
const MI_NOT_SOLVED = 1
const MI_FREE = 0
const MI_ONE = 1
const MI_ZERO = 0
const MI_STAR = -1
const MI_OPTIMAL_SOLN = ECOS_OPTIMAL
const MI_INFEASIBLE = ECOS_PINF
const MI_UNBOUNDED = ECOS_DINF
const MI_MAXITER_FEASIBLE_SOLN = ECOS_OPTIMAL + ECOS_INACC_OFFSET
const MI_MAXITER_NO_SOLN = ECOS_PINF + ECOS_INACC_OFFSET
const MI_MAXITER_UNBOUNDED = ECOS_DINF + ECOS_INACC_OFFSET
const MAX_FLOAT_INT = 8388608

@cenum BRANCHING_STRATEGY::UInt32 begin
    BRANCHING_STRATEGY_MOST_INFEASIBLE = 0
    BRANCHING_STRATEGY_STRONG_BRANCHING = 1
    BRANCHING_STRATEGY_PSEUDOCOST_BRANCHING = 2
    BRANCHING_STRATEGY_RELIABILITY = 3
    BRANCHING_STRATEGY_RANDOM = 4
end

@cenum NODE_SELECTION_METHOD::UInt32 begin
    BREADTH_FIRST = 0
    DIVE_LOWER_NODE = 1
    DIVE_UPPER_NODE = 2
end


struct settings_bb
    maxit::idxint
    verbose::idxint
    abs_tol_gap::pfloat
    rel_tol_gap::pfloat
    integer_tol::pfloat
    branching_strategy::BRANCHING_STRATEGY
    reliable_eta::idxint
    node_selection_method::NODE_SELECTION_METHOD
end

struct node
    status::UInt8
    L::pfloat
    U::pfloat
    relaxation::pfloat
    split_idx::idxint
    split_val::pfloat
    prev_split_idx::idxint
    prev_split_val::pfloat
    prev_relaxation::pfloat
    up_branch_node::Cint
end

struct ecos_bb_pwork
    num_bool_vars::idxint
    num_int_vars::idxint
    nodes::Ptr{node}
    bool_node_ids::Cstring
    int_node_ids::Ptr{pfloat}
    bool_vars_idx::Ptr{idxint}
    int_vars_idx::Ptr{idxint}
    ecos_prob::Ptr{pwork}
    A::Ptr{spmat}
    G::Ptr{spmat}
    c::Ptr{pfloat}
    b::Ptr{pfloat}
    h::Ptr{pfloat}
    x::Ptr{pfloat}
    y::Ptr{pfloat}
    z::Ptr{pfloat}
    s::Ptr{pfloat}
    kap::pfloat
    tau::pfloat
    info::Ptr{stats}
    global_U::pfloat
    global_L::pfloat
    tmp_bool_node_id::Cstring
    tmp_int_node_id::Ptr{pfloat}
    iter::idxint
    dive_node_id::idxint
    tmp_branching_bool_node_id::Cstring
    tmp_branching_int_node_id::Ptr{pfloat}
    pseudocost_bin_sum::Ptr{pfloat}
    pseudocost_int_sum::Ptr{pfloat}
    pseudocost_bin_cnt::Ptr{idxint}
    pseudocost_int_cnt::Ptr{idxint}
    Gpr_new::Ptr{pfloat}
    Gjc_new::Ptr{idxint}
    Gir_new::Ptr{idxint}
    h_new::Ptr{pfloat}
    ecos_stgs::Ptr{settings}
    stgs::Ptr{settings_bb}
    default_settings::idxint
end

const PRINTLEVEL = 2
const PROFILING = 1
const DEBUG = 0
const ECOS_INFINITY = Inf
const ECOS_NAN = NaN
# const PRINTTEXT = printf
# const MALLOC = malloc
# const FREE = free
const KKT_PROBLEM = 0
const KKT_OK = 1
