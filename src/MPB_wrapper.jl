#############################################################################
# ECOS.jl
# Wrapper around the ECOS solver https://github.com/ifa-ethz/ecos
# See http://github.com/jump-dev/ECOS.jl
#############################################################################
# ECOSSolverInterface.jl
# MathProgBase.jl interface for the ECOS.jl solver wrapper
#############################################################################

import MathProgBase
const MPB = MathProgBase

#############################################################################
# Define the MPB Solver and Model objects
export ECOSSolver
struct ECOSSolver <: MPB.AbstractMathProgSolver
    options
end
ECOSSolver(;kwargs...) = ECOSSolver(kwargs)

mutable struct ECOSMathProgModel <: MPB.AbstractConicModel
    nvar::Int                           # Number of variables
    nconstr::Int                        # Number of rows in the MPB-form A matrix
    nineq::Int                          # Number of inequalities Gx <=_K h
    neq::Int                            # Number of equalities Ax = b
    npos::Int                           # Number of positive orthant cones
    ncones::Int                         # Number of SO cones
    conedims::Vector{Int}               # Dimension of each SO cone
    nexp_cones::Int                     # Number of exponential cones
    G::ECOSMatrix                       # The G matrix (inequalties)
    A::ECOSMatrix                       # The A matrix (equalities)
    c::Vector{Float64}                  # The objective coeffs (always min)
    orig_sense::Symbol                  # Original objective sense
    h::Vector{Float64}                  # RHS for inequality
    b::Vector{Float64}                  # RHS for equality
    # Post-solve
    solve_stat::Symbol
    solve_time::Float64
    obj_val::Float64
    primal_sol::Vector{Float64}
    dual_sol_eq::Vector{Float64}
    dual_sol_ineq::Vector{Float64}
    # Maps x∈K to ECOS duals
    # .._ind maps a column to a row index
    # .._type  is the original cone type of each column
    col_map_ind::Vector{Int}
    col_map_type::Vector{Symbol}
    # Maps b-Ax∈K to ECOS duals
    # .._ind maps a row to an index
    # .._type  is the original cone type of each row
    row_map_ind::Vector{Int}
    row_map_type::Vector{Symbol}
    # To reorder solution to MPB form
    fwd_map::Vector{Int}
    options
end
ECOSMathProgModel(;kwargs...) = ECOSMathProgModel(0,0,0,0,0,0,
                                        Int[],0,
                                        ECOSMatrix(spzeros(0,0)),
                                        ECOSMatrix(spzeros(0,0)),
                                        Float64[], :Min,
                                        Float64[], Float64[],
                                        :NotSolved, 0.0, 0.0,
                                        Float64[],
                                        Float64[], Float64[],
                                        Int[], Symbol[],
                                        Int[], Symbol[], Int[],
                                        kwargs)

#############################################################################
# Begin implementation of the MPB low-level interface
# Implements
# - model
# - optimize!
# - status
# - numvar
# - numconstr
# http://mathprogbasejl.readthedocs.org/en/latest/lowlevel.html

MPB.ConicModel(s::ECOSSolver) = ECOSMathProgModel(;s.options...)
MPB.LinearQuadraticModel(s::ECOSSolver) = MPB.ConicToLPQPBridge(MPB.ConicModel(s))


function MPB.optimize!(m::ECOSMathProgModel)
    ecos_prob_p = setup(
        m.nvar, m.nineq, m.neq,
        m.npos, m.ncones, m.conedims, m.nexp_cones,
        m.G, m.A,
        m.c,
        m.h, m.b; m.options...)
    # Note: ECOS modifies problem data in setup() and restores it on cleanup()

    t = time()
    flag = ECOS_solve(ecos_prob_p)
    m.solve_time = time() - t
    if flag == ECOS_OPTIMAL
        m.solve_stat = :Optimal
    elseif flag == ECOS_PINF
        m.solve_stat = :Infeasible
    elseif flag == ECOS_DINF  # Dual infeasible = primal unbounded, probably
        m.solve_stat = :Unbounded
    elseif flag == ECOS_MAXIT
        m.solve_stat = :UserLimit
    elseif flag == ECOS_OPTIMAL + ECOS_INACC_OFFSET
        m.solve_stat = :Suboptimal
    else
        m.solve_stat = :Error
    end
    # Extract solution
    ecos_prob = unsafe_load(ecos_prob_p.ptr)
    m.primal_sol = copy(unsafe_wrap(Array, ecos_prob.x, m.nvar))
    m.dual_sol_eq   = copy(unsafe_wrap(Array, ecos_prob.y, m.neq))
    m.dual_sol_ineq = copy(unsafe_wrap(Array, ecos_prob.z, m.nineq))
    # Compute this value after cleanup with restored c.
    m.obj_val = dot(m.c, m.primal_sol) * (m.orig_sense == :Max ? -1 : +1)
    return
end

MPB.status(m::ECOSMathProgModel) = m.solve_stat
MPB.getobjval(m::ECOSMathProgModel) = m.obj_val
MPB.getsolution(m::ECOSMathProgModel) = m.primal_sol[m.fwd_map]

MPB.numvar(m::ECOSMathProgModel) = m.nvar
MPB.numconstr(m::ECOSMathProgModel) = m.nineq + m.neq

#############################################################################
# Begin implementation of the MPB conic interface
# Implements
# - supportedcones
# - loadconicproblem!
# - getconicdual
# http://mathprogbasejl.readthedocs.org/en/latest/conic.html

MPB.supportedcones(m::ECOSSolver) = [:Free,:Zero,:NonNeg,:NonPos,:SOC,:ExpPrimal]

function MPB.loadproblem!(m::ECOSMathProgModel, c, A, b, constr_cones, var_cones)
    if size(A,2) == 0
        Base.warn("Input matrix has no columns. ECOS is known to crash in this corner case.")
    end
    m.nconstr = size(A,1)
    # If A is sparse, we should use an appropriate "zeros"
    zeromat = isa(A,SparseMatrixCSC) ? spzeros : zeros

    # We don't support SOCRotated, SDP, or ExpDual
    bad_cones = [:SOCRotated, :SDP, :ExpDual]
    for cone_vars in constr_cones
        cone_vars[1] in bad_cones && error("Cone type $(cone_vars[1]) not supported")
    end
    for cone_vars in var_cones
        cone_vars[1] in bad_cones && error("Cone type $(cone_vars[1]) not supported")
    end

    # MathProgBase form             ECOS form
    # min  c'x                      min  c'x
    #  st b-Ax ∈ K_1                 st   Ax = b
    #        x ∈ K_2                    h-Gx ∈ K
    #
    # Mapping:
    # * For the constaints (K_1)
    #   * If :Zero cone, then treat as equality constraint in ECOS form
    #   * Otherewise trivially maps to h-Gx in ECOS form
    # * For the variables (K_2)
    #   * If :Free, do nothing
    #   * If :Zero, put in as equality constraint
    #   * If rest, stick in h-Gx
    #
    # Approach:
    # 0. Figure out a mapping between MPB and ECOS
    # 1. Map non-SOC/Exp variables to ECOS form
    # 2. Map non-SOC/Exp constraints to ECOS form
    # 3. Map SOC variables and constraints to ECOS form
    # 4. Map Exp variables and constraints to ECOS form

    # Allocate space for the ECOS variables
    num_vars = length(c)
    fwd_map = Array{Int}(undef, num_vars)  # Will be used for SOCs
    rev_map = Array{Int}(undef, num_vars)  # Need to restore sol. vec.
    idxcone = Array{Symbol}(undef, num_vars)  # We'll use this for non-SOC/Exp
    m.col_map_type = Array{Symbol}(undef, num_vars) # In MPB indices

    # Now build the mapping between MPB variables and ECOS variables
    pos = 1
    for (cone, idxs) in var_cones
        for i in idxs
            fwd_map[i]   = pos   # fwd_map = MPB idx -> ECOS idx
            rev_map[pos] = i     # rev_map = ECOS idx -> MPB idx
            idxcone[pos] = cone
            m.col_map_type[i] = cone
            pos += 1
        end
    end

    # Rearrange data into the internal ordering, make copy
    ecos_c = c[rev_map]
    ecos_A = zeromat(0,num_vars)
    ecos_b = Float64[]

    # Mapping for duals
    m.col_map_ind = zeros(Int,num_vars)
    m.row_map_ind = zeros(Int, length(b))
    m.row_map_type  = Array{Symbol}(undef, length(b))
    function update_row_map(cone_type, cur_ind)
        for (cone,idxs) in constr_cones
            if cone == cone_type
                for idx in idxs
                    m.row_map_ind[idx] = cur_ind
                    m.row_map_type[idx] = cone_type
                    cur_ind += 1
                end
            end
        end
        cur_ind
    end

    ###################################################################
    # PHASE ONE  -  MAP x ∈ K_2 to ECOS form, except SOC/Exp

    # If a variable is in :Zero cone, fix at 0 with equality constraint.
    for j = 1:num_vars
        idxcone[j] != :Zero && continue

        new_row    = zeromat(1,num_vars)
        new_row[j] = 1.0
        ecos_A     = vcat(ecos_A, new_row)
        ecos_b     = vcat(ecos_b, 0.0)
        m.col_map_ind[rev_map[j]] = length(ecos_b)
    end

    # G matrix:
    # * 1 row ∀ :NonNeg & :NonPos cones
    # * 1 row ∀ variable in SOC/Exp cones
    # We will only handle the first case here, the rest in phase 3.
    num_G_row_negpos = 0
    for j = 1:num_vars
        !(idxcone[j] == :NonNeg || idxcone[j] == :NonPos) && continue
        num_G_row_negpos += 1
        m.col_map_ind[rev_map[j]] = num_G_row_negpos
    end
    ecos_G = zeromat(num_G_row_negpos,num_vars)
    ecos_h =   zeros(num_G_row_negpos)

    # Handle the :NonNeg, :NonPos cases
    num_pos_orth = 0
    G_row = 1
    for j = 1:num_vars
        if idxcone[j] == :NonNeg
            ecos_G[G_row,j] = -1.0
            G_row += 1
            num_pos_orth += 1
        elseif idxcone[j] == :NonPos
            ecos_G[G_row,j] = +1.0
            G_row += 1
            num_pos_orth += 1
        end
    end
    @assert G_row == num_pos_orth + 1

    ###################################################################
    # PHASE TWO  -  MAP b-Ax ∈ K_1 to ECOS form, except SOC/Exp

    # Zero rows for Ax=b, NonNegPos rows to append to G,h
     eq_rows = Int[]
    pos_rows = Int[]
    neg_rows = Int[]
    for (cone,idxs) in constr_cones
        cone == :Free && error("Free cone constraints not handled")
        (cone == :SOC || cone == :ExpPrimal) && continue  # Handle later
        idxset = vec(collect(idxs))
        if cone == :Zero
            append!(eq_rows, idxset)
            continue
        end
        cone == :NonNeg && append!(pos_rows, idxset)
        cone == :NonPos && append!(neg_rows, idxset)
    end
    # Update mappings - eq, nonneg, nonpos
     eq_cur_ind = length(ecos_b) + 1
     eq_cur_ind = update_row_map(:Zero, eq_cur_ind)
    neq_cur_ind = length(ecos_h) + 1
    neq_cur_ind = update_row_map(:NonNeg, neq_cur_ind)
    neq_cur_ind = update_row_map(:NonPos, neq_cur_ind)
    # Equality constraints / Zero cones
    ecos_A = vcat(ecos_A,  A[eq_rows,rev_map])
    ecos_b = vcat(ecos_b,  b[eq_rows])
    ecos_G = vcat(ecos_G,  A[pos_rows,rev_map])
    ecos_h = vcat(ecos_h,  b[pos_rows])
    ecos_G = vcat(ecos_G, -A[neg_rows,rev_map])  # b-a'x <= 0 - flip sign,
    ecos_h = vcat(ecos_h, -b[neg_rows])          # then maps to a row in h-Gx
    G_row        += length(pos_rows) + length(neg_rows)
    num_pos_orth += length(pos_rows) + length(neg_rows)

    ###################################################################
    # PHASE THREE  -  MAP SOC variables and constraints to ECOS form

    # Handle the SOC variable cones
    # MPB  form: vector of var (y,x) is in the SOC ||x|| <= y
    # ECOS form: h - Gx ∈ Q  -->  0 - Ix ∈ Q
    num_G_row_soc = 0
    num_G_row_exp = 0
    for j = 1:num_vars
        if idxcone[j] == :SOC
            num_G_row_soc += 1
        elseif idxcone[j] == :ExpPrimal
            num_G_row_exp += 1
        end
    end
    ecos_G = vcat(ecos_G, zeromat(num_G_row_soc,num_vars))
    ecos_h = vcat(ecos_h,   zeros(num_G_row_soc))

    num_SOC_cones = 0
    SOC_conedims  = Int[]
    for (cone, idxs) in var_cones
        cone != :SOC && continue
        # Found a new SOC
        num_SOC_cones += 1
        push!(SOC_conedims, length(idxs))
        # Add the entries (carrying on from pos. orthant)
        for j in idxs
            ecos_G[G_row,fwd_map[j]] = -1.0
            m.col_map_ind[j] = G_row
            G_row += 1
        end
    end
    @assert G_row == num_pos_orth + num_G_row_soc + 1


    # Handle the SOC constraint cones
    # Collect all the rows we'll be appending to G,h
    all_rows = Int[]
    for (cone,idxs) in constr_cones
        if cone == :SOC
            num_SOC_cones += 1
            push!(SOC_conedims, length(idxs))
            idx_list   = collect(idxs)
            all_rows   = vcat(all_rows,   idx_list)
        end
    end
    neq_cur_ind = update_row_map(:SOC, neq_cur_ind)
    ecos_G = vcat(ecos_G, A[all_rows,rev_map])
    ecos_h = vcat(ecos_h, b[all_rows])

    ###################################################################
    # PHASE FOUR  -  MAP Exp variables and constraints to ECOS form

    # Note: The MPB definition of the ExpPrimal cone and the ECOS
    # defintions swap the 2nd and 3rd components, so we have to
    # manually permute the indices.

    num_exp_cones = 0
    ecos_G = vcat(ecos_G, zeromat(num_G_row_exp,num_vars))
    ecos_h = vcat(ecos_h,   zeros(num_G_row_exp))

    G_row += sum(all_rows)
    for (cone, idxs) in var_cones
        cone != :ExpPrimal && continue
        # Found a new Exp
        num_exp_cones += 1
        @assert length(idxs) == 3
        # Add the entries in permuted order
        ecos_G[G_row,fwd_map[idxs[1]]] = -1.0
        ecos_G[G_row+2,fwd_map[idxs[2]]] = -1.0
        ecos_G[G_row+1,fwd_map[idxs[3]]] = -1.0
        m.col_map_ind[idxs[1]] = G_row
        m.col_map_ind[idxs[2]] = G_row+2
        m.col_map_ind[idxs[3]] = G_row+1
        G_row += 3
    end
    @assert G_row == num_pos_orth + num_G_row_soc + sum(all_rows) + num_G_row_exp + 1

    exp_rows = Int[]
    for (cone,idxs) in constr_cones
        if cone == :ExpPrimal
            num_exp_cones += 1
            @assert length(idxs) == 3
            append!(exp_rows, [idxs[1],idxs[3],idxs[2]])

            # override update_row_map to handle the permutation
            m.row_map_type[idxs] .= :ExpPrimal
            m.row_map_ind[idxs[1]] = neq_cur_ind
            m.row_map_ind[idxs[3]] = neq_cur_ind+1
            m.row_map_ind[idxs[2]] = neq_cur_ind+2
            neq_cur_ind += 3
        end
    end
    ecos_G = vcat(ecos_G, A[exp_rows,rev_map])
    ecos_h = vcat(ecos_h, b[exp_rows])

    ###################################################################
    # Store in the ECOS structure
    m.nvar          = num_vars          # Num variable
    m.nineq         = size(ecos_G,1)    # Num inequality constraints
    m.neq           = length(ecos_b)    # Num equality constraints
    m.npos          = num_pos_orth      # Num ineq. constr. in +ve orthant
    m.ncones        = num_SOC_cones     # Num second-order cones
    m.conedims      = SOC_conedims      # Num contr. in each SOC
    m.nexp_cones    = num_exp_cones     # Num exponential cones
    m.G             = ECOSMatrix(ecos_G)
    m.A             = ECOSMatrix(ecos_A)
    m.c             = ecos_c
    m.orig_sense    = :Min
    m.h             = ecos_h
    m.b             = ecos_b
    m.fwd_map       = fwd_map           # Used to return solution
end


function MPB.getdual(m::ECOSMathProgModel)
    duals = zeros(length(m.row_map_ind))
    for (mpb_row,ecos_row) in enumerate(m.row_map_ind)
        cone = m.row_map_type[mpb_row]
        if cone == :Zero
            # This MPB constraint ended up in ECOS equality block
            duals[mpb_row] = m.dual_sol_eq[ecos_row]
        else
            # Ended up in ECOS inequality block
            if cone == :NonPos
                duals[mpb_row] = -m.dual_sol_ineq[ecos_row]
            else
                duals[mpb_row] = m.dual_sol_ineq[ecos_row]
            end
        end
    end
    return duals
end


function MPB.getvardual(m::ECOSMathProgModel)
    duals = zeros(length(m.col_map_ind))
    for (mpb_col,ecos_row) in enumerate(m.col_map_ind)
        cone = m.col_map_type[mpb_col]
        cone == :Free && continue # dual is zero
        if cone == :Zero
            # This MPB zero var ended up in ECOS equality block
            duals[mpb_col] = m.dual_sol_eq[ecos_row]
        else
            # Ended up in ECOS inequality block
            if cone == :NonPos
                duals[mpb_col] = -m.dual_sol_ineq[ecos_row]
            else
                duals[mpb_col] = m.dual_sol_ineq[ecos_row]
            end
        end
    end
    return duals
end

MPB.getsolvetime(m::ECOSMathProgModel) = m.solve_time
