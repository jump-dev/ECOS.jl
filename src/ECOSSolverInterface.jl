#############################################################################
# ECOS.jl
# Wrapper around the ECOS solver https://github.com/ifa-ethz/ecos
# See http://github.com/JuliaOpt/ECOS.jl
#############################################################################
# ECOSSolverInterface.jl
# MathProgBase.jl interface for the ECOS.jl solver wrapper
#############################################################################

importall MathProgBase.MathProgSolverInterface

#############################################################################
# Define the MPB Solver and Model objects
export ECOSSolver
immutable ECOSSolver <: AbstractMathProgSolver
    options
end
ECOSSolver(;kwargs...) = ECOSSolver(kwargs)

type ECOSMathProgModel <: AbstractMathProgModel
    nvar::Int                           # Number of variables
    nineq::Int                          # Number of inequalities Gx <=_K h
    neq::Int                            # Number of equalities Ax = b
    npos::Int                           # Number of positive orthant cones
    ncones::Int                         # Number of SO cones
    conedims::Vector{Int}               # Dimension of each SO cone
    G::SparseMatrixCSC{Float64,Int}     # The G matrix (inequalties)
    A::SparseMatrixCSC{Float64,Int}     # The A matrix (equalities)
    c::Vector{Float64}                  # The objective coeffs (always min)
    orig_sense::Symbol                  # Original objective sense
    h::Vector{Float64}                  # RHS for inequality 
    b::Vector{Float64}                  # RHS for equality
    # Post-solve
    solve_stat::Symbol
    obj_val::Float64
    primal_sol::Vector{Float64}
    dual_sol_eq::Vector{Float64}
    dual_sol_ineq::Vector{Float64}
    # Maps b-Ax∈K to ECOS duals
    # .._ind maps a row to an index
    # .._type  is the original cone type of each row
    row_map_ind::Vector{Int}
    row_map_type::Vector{Symbol}
    # To reorder solution if we solved using the conic interface
    fwd_map::Vector{Int}
    options
end
ECOSMathProgModel(;kwargs...) = ECOSMathProgModel(0,0,0,0,0,
                                        Int[],
                                        spzeros(0,0),
                                        spzeros(0,0),
                                        Float64[], :Min,
                                        Float64[], Float64[],
                                        :NotSolved, 0.0, 
                                        Float64[],
                                        Float64[], Float64[], 
                                        Int[], Symbol[], Int[],
                                        kwargs)

#############################################################################
# Begin implementation of the MPB low-level interface 
# Implements
# - model
# - loadproblem!
# - optimize!
# - status
# http://mathprogbasejl.readthedocs.org/en/latest/lowlevel.html

model(s::ECOSSolver) = ECOSMathProgModel(;s.options...)

# Loads the provided problem data to set up the linear programming problem:
# min c'x
# st  lb <= Ax <= ub
#      l <=  x <= u
# where sense = :Min or :Max
function loadproblem!(m::ECOSMathProgModel, A, collb, colub, obj, rowlb, rowub, sense)
    (nvar = length(collb)) == length(colub) || error("Unequal lengths for column bounds")
    (nrow = length(rowlb)) == length(rowub) || error("Unequal lengths for row bounds")
    
    # Turn variable bounds into constraints
    # Inefficient, because keeps allocating memory!
    # Would need to batch, get tricky...
    for j = 1:nvar
        if collb[j] != -Inf
            # Variable has lower bound
            newrow = zeros(1, nvar)
            newrow[j] = -1.0
            A = vcat(A, newrow)
            rowlb = vcat(rowlb, -Inf)
            rowub = vcat(rowub, -collb[j])
            nrow += 1
        end
        if colub[j] != +Inf
            # Variable has upper bound
            newrow = zeros(1, nvar)
            newrow[j] = 1.0
            A = vcat(A, newrow)
            rowlb = vcat(rowlb, -Inf)
            rowub = vcat(rowub, colub[j])
            nrow += 1
        end
    end

    eqidx   = Int[]      # Equality row indicies
    ineqidx = Int[]      # Inequality row indicies
    eqbnd   = Float64[]  # Bounds for equality rows
    ineqbnd = Float64[]  # Bounds for inequality row
    for it in 1:nrow
        # Equality constraint
        if rowlb[it] == rowub[it]
            push!(eqidx, it)
            push!(eqbnd, rowlb[it])
        # Range constraint - not supported
        elseif rowlb[it] != -Inf && rowub[it] != Inf
            error("Ranged constraints unsupported!")
        # Less-than constraint
        elseif rowlb[it] == -Inf
            push!(ineqidx, it)
            push!(ineqbnd, rowub[it])
        # Greater-than constraint - flip sign so only have <= constraints
        else
            push!(ineqidx, it)
            push!(ineqbnd, -rowlb[it])
            A[it,:] *= -1 # flip signs so we have Ax<=b
        end
    end

    m.nvar      = nvar                  # Number of variables
    m.nineq     = length(ineqidx)       # Number of inequalities Gx <=_K h
    m.neq       = length(eqidx)         # Number of equalities Ax = b
    m.npos      = length(ineqidx)       # Number of positive orthant cone
    m.ncones    = 0                     # Number of SO cones
    m.conedims  = Int[]                 # Dimenions of SO cones
    m.G         = sparse(A[ineqidx,:])  # The G matrix (inequalties)
    m.A         = sparse(A[eqidx,:])    # The A matrix (equalities)
    m.c         = (sense == :Max) ? obj * -1 : obj[:] 
                                        # The objective coeffs (always min)
    m.orig_sense = sense                # Original objective sense
    m.h         = ineqbnd               # RHS for inequality 
    m.b         = eqbnd                 # RHS for equality
    m.fwd_map   = [1:nvar]              # Identity mapping
end

function optimize!(m::ECOSMathProgModel)
    ecos_prob_ptr = setup(
        m.nvar, m.nineq, m.neq,
        m.npos, m.ncones, m.conedims,
        m.G, m.A,
        m.c[:],  # Seems to modify this
        m.h, m.b; m.options...)

    flag = solve(ecos_prob_ptr)
    if flag == ECOS_OPTIMAL
        m.solve_stat = :Optimal
    elseif flag == ECOS_PINF
        m.solve_stat = :Infeasible
    elseif flag == ECOS_DINF  # Dual infeasible = primal unbounded, probably
        m.solve_stat = :Unbounded
    elseif flag == ECOS_MAXIT
        m.solve_stat = :UserLimit
    else
        m.solve_stat = :Error
    end
    # Extract solution
    ecos_prob = pointer_to_array(ecos_prob_ptr, 1)[1]
    m.primal_sol = pointer_to_array(ecos_prob.x, m.nvar)[:]
    m.dual_sol_eq   = pointer_to_array(ecos_prob.y, m.neq)[:]
    m.dual_sol_ineq = pointer_to_array(ecos_prob.z, m.nineq)[:]
    m.obj_val = dot(m.c, m.primal_sol) * (m.orig_sense == :Max ? -1 : +1)  
    cleanup(ecos_prob_ptr, 0)
end

status(m::ECOSMathProgModel) = m.solve_stat
getobjval(m::ECOSMathProgModel) = m.obj_val
getsolution(m::ECOSMathProgModel) = m.primal_sol[m.fwd_map]

#############################################################################
# Begin implementation of the MPB conic interface 
# Implements
# - loadconicproblem!
# http://mathprogbasejl.readthedocs.org/en/latest/conic.html

function loadconicproblem!(m::ECOSMathProgModel, c, A, b, constr_cones, var_cones)
    # If A is sparse, we should use an appropriate "zeros"
    const zeromat = isa(A,SparseMatrixCSC) ? spzeros : zeros

    # We don't support SOCRotated, SDP, or Exp*
    bad_cones = [:SOCRotated, :SDP, :ExpPrimal, :ExpDual]
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
    # 1. Map non-SOC variables to ECOS form
    # 2. Map non-SOC constraints to ECOS form
    # 3. Map SOC variables to ECOS form
    # 4. Map SOC constraints to ECOS form

    # Allocate space for the ECOS variables
    num_vars = length(c)
    fwd_map = Array(Int,    num_vars)  # Will be used for SOCs
    rev_map = Array(Int,    num_vars)  # Need to restore sol. vec.
    idxcone = Array(Symbol, num_vars)  # We'll uses this for non-SOCs

    # Now build the mapping between MPB variables and ECOS variables
    pos = 1
    for (cone, idxs) in var_cones
        for i in idxs
            fwd_map[i]   = pos   # fwd_map = MPB idx -> ECOS idx
            rev_map[pos] = i     # rev_map = ECOS idx -> MPB idx
            idxcone[pos] = cone
            pos += 1
        end
    end

    # Rearrange data into the internal ordering, make copy
    ecos_c = c[rev_map]
    ecos_A = zeromat(0,num_vars)
    ecos_b = Float64[]

    # Mapping for duals
    m.row_map_ind = zeros( Int, length(b))
    m.row_map_type  = Array(Symbol, length(b))
    function update_map(cone_type, cur_ind)
        for (cone,idxs) in constr_cones
            if cone == cone_type
                for idx in idxs
                    m.row_map_ind[idx] = cur_ind
                    m.row_map_type[idx]  = cone_type
                    cur_ind += 1
                end
            end
        end
        cur_ind
    end

    ###################################################################
    # PHASE ONE  -  MAP x ∈ K_2 to ECOS form, except SOC
    
    # If a variable is in :Zero cone, fix at 0 with equality constraint.
    for j = 1:num_vars
        idxcone[j] != :Zero && continue

        new_row    = zeromat(1,num_vars)
        new_row[j] = 1.0
        ecos_A     = vcat(ecos_A, new_row)
        ecos_b     = vcat(ecos_b, 0.0)
    end

    # G matrix:
    # * 1 row ∀ :NonNeg & :NonPos cones
    # * 1 row ∀ variable in :SOC cone
    # We will only handle the first case here, the rest in phase 3.
    num_G_row_negpos = 0
    for j = 1:num_vars
        !(idxcone[j] == :NonNeg || idxcone[j] == :NonPos) && continue
        num_G_row_negpos += 1
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
    # PHASE TWO  -  MAP b-Ax ∈ K_1 to ECOS form, except SOC
    
    # Zero rows for Ax=b, NonNegPos rows to append to G,h
     eq_rows = Int[]
    pos_rows = Int[]
    neg_rows = Int[]
    for (cone,idxs) in constr_cones
        cone == :Free && error("Free cone constraints not handled")
        cone == :SOC  && continue  # Handle later
        idxset = Int[idxs]
        if cone == :Zero
            append!(eq_rows, idxset)
            continue
        end
        cone == :NonNeg && append!(pos_rows, idxset)
        cone == :NonPos && append!(neg_rows, idxset)
    end
    # Update mappings - eq, nonneg, nonpos
     eq_cur_ind = length(ecos_b) + 1
     eq_cur_ind = update_map(:Zero, eq_cur_ind)
    neq_cur_ind = length(ecos_h) + 1
    neq_cur_ind = update_map(:NonNeg, neq_cur_ind)
    neq_cur_ind = update_map(:NonPos, neq_cur_ind)
    # Equality constraints / Zero cones
    ecos_A = vcat(ecos_A,  A[eq_rows,rev_map])
    ecos_b = vcat(ecos_b,  b[eq_rows])
    ecos_G = vcat(ecos_G,  A[pos_rows,rev_map])
    ecos_h = vcat(ecos_h,  b[pos_rows])
    ecos_G = vcat(ecos_G, -A[neg_rows,rev_map])  # b-a'x <= 0 - flip sign, 
    ecos_h = vcat(ecos_h, -b[neg_rows])          # then maps to a row in h-Gx
    G_row         += length(pos_rows) + length(neg_rows)
    num_pos_orth  += length(pos_rows) + length(neg_rows)

    ###################################################################
    # PHASE THREE  -  MAP x ∈ SOC to ECOS form
    
    # Handle the :SOC variable cones
    # MPB  form: vector of var (y,x) is in the SOC ||x|| <= y
    # ECOS form: h - Gx ∈ Q  -->  0 - Ix ∈ Q
    num_G_row_soc = 0
    for j = 1:num_vars
        (idxcone[j] != :SOC) && continue
        num_G_row_soc += 1
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
            G_row += 1
        end
    end
    @assert G_row == num_pos_orth + num_G_row_soc + 1

    ###################################################################
    # PHASE FOUR  -  MAP b-Ax ∈ SOC to ECOS form
    
    # Collect all the rows we'll be appending to G,h
    all_rows = Int[]
    for (cone,idxs) in constr_cones
        if cone == :SOC
            num_SOC_cones += 1
            push!(SOC_conedims, length(idxs))
            idx_list   = Int[idxs]
            all_rows   = vcat(all_rows,   idx_list)
        end
    end
    update_map(:SOC, neq_cur_ind)
    ecos_G = vcat(ecos_G, A[all_rows,rev_map])
    ecos_h = vcat(ecos_h, b[all_rows])

    ###################################################################
    # Store in the ECOS structure
    m.nvar          = num_vars          # Num variable
    m.nineq         = size(ecos_G,1)    # Num inequality constraints
    m.neq           = length(ecos_b)    # Num equality constraints
    m.npos          = num_pos_orth      # Num ineq. constr. in +ve orthant
    m.ncones        = num_SOC_cones     # Num second-order cones
    m.conedims      = SOC_conedims      # Num contr. in each SOC
    m.G             = ecos_G
    m.A             = ecos_A
    m.c             = ecos_c
    m.orig_sense    = :Min
    m.h             = ecos_h
    m.b             = ecos_b
    m.fwd_map       = fwd_map           # Used to return solution
end


function getconicdual(m::ECOSMathProgModel)
    #@show m.row_map_ind
    #@show m.row_map_eq
    duals = zeros(length(m.row_map_ind))
    for (mpb_row,ecos_row) in enumerate(m.row_map_ind)
        cone = m.row_map_type[mpb_row]
        # Not well understood why we need to flip signs
        # of Zero and SOC constraints.
        if cone == :Zero
            # This MPB constraint ended up in ECOS equality block
            duals[mpb_row] = -m.dual_sol_eq[ecos_row]
        else
            # Ended up in ECOS inequality block
            if cone == :NonPos || cone == :SOC
                duals[mpb_row] = -m.dual_sol_ineq[ecos_row]
            else
                duals[mpb_row] = m.dual_sol_ineq[ecos_row]
            end
        end
    end
    return duals
end
