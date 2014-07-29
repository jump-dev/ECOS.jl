#############################################################################
# ECOS.jl
# Wrapper around the ECOS solver https://github.com/ifa-ethz/ecos
# See http://github.com/JuliaOpt/ECOS.jl
#############################################################################
# ECOSSolverInterface.jl
# MathProgBase.jl interface for the ECOS.jl solver wrapper
#############################################################################

require(joinpath(Pkg.dir("MathProgBase"),"src","MathProgSolverInterface.jl"))
importall MathProgSolverInterface

#############################################################################
# Define the MPB Solver and Model objects
export ECOSSolver
immutable ECOSSolver <: AbstractMathProgSolver
end

type ECOSMathProgModel <: AbstractMathProgModel
    nvar::Int                           # Number of variables
    nineq::Int                          # Number of inequalities Gx <=_K h
    neq::Int                            # Number of equalities Ax = b
    npos::Int                           # Number of ???
    ncones::Int                         # Number of SO cones
    conedims::Vector{Int}               # ?
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
    fwd_map::Vector{Int}                # To reorder solution if we solved
end                                     # using the conic interface
ECOSMathProgModel() = ECOSMathProgModel(0,0,0,0,0,
                                        Int[],
                                        spzeros(0,0),
                                        spzeros(0,0),
                                        Float64[], :Min,
                                        Float64[], Float64[],
                                        :NotSolved, 0.0, Float64[], Int[])

#############################################################################
# Begin implementation of the MPB low-level interface 
# Implements
# - model
# - loadproblem!
# - optimize!
# - status
# http://mathprogbasejl.readthedocs.org/en/latest/lowlevel.html

model(s::ECOSSolver) = ECOSMathProgModel()

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
    m.npos      = length(ineqidx)       # Number of ???
    m.ncones    = 0                     # Number of SO cones
    m.conedims  = Int[]                 # ???
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
        n       = m.nvar,
        m       = m.nineq,
        p       = m.neq,
        l       = m.npos,
        ncones  = m.ncones,
        q       = m.conedims,
        G       = m.G,
        A       = m.A,
        c       = m.c[:],  # Seems to modify this
        h       = m.h,
        b       = m.b)

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

function loadconicproblem!(m::ECOSMathProgModel, c, A, b, cones)
    # TODO (if it matters): make this more efficient for sparse A

    # We don't support SOCRotated, SDP, or Exp*
    bad_cones = [:SOCRotated, :SDP, :ExpPrimal, :ExpDual]
    for cone_vars in cones
        cone_vars[1] in bad_cones && error("Cone type $(cone_vars[1]) not supported")
    end

    # MathProgBase form             ECOS form
    # min c'x                       min c'x
    # st  A x = b                   st  A x = b
    #       x in K                      h - Gx in K

    # Expand out the cones info
    # The cones can come in any order, so we need to build a mapping
    # from the variable indices in the input to the internal ordering
    # we will use.
    
    # In the first past we'll just count up the number of variables 
    # of each type.
    num_vars = 0
    for (cone_type, idxs) in cones
        num_vars += length(idxs)
    end
    fwd_map = Array(Int,    num_vars)  # Will be used for SOCs
    rev_map = Array(Int,    num_vars)  # Need to restore sol. vec.
    idxcone = Array(Symbol, num_vars)  # We'll uses this for non-SOCs

    # Now build the mapping
    pos = 1
    for (cone, idxs) in cones
        for i in idxs
            fwd_map[i]   = pos   # fwd_map = orig idx -> internal idx
            rev_map[pos] = i     # rev_map = internal idx -> orig idx
            idxcone[pos] = cone
            pos += 1
        end
    end

    # Rearrange data into the internal ordering
    ecos_c = c[rev_map]
    ecos_A = A[:,rev_map]
    ecos_b = b[:]

    # For all variables in the :Zero cone, fix at 0 with an
    # equality constraint. TODO: Don't even include them

    for j = 1:num_vars
        idxcone[j] != :Zero && continue

        new_row    = zeros(1,num_vars)
        new_row[j] = 1.0
        ecos_A     = vcat(ecos_A, new_row)
        ecos_b     = vcat(ecos_b, 0.0)
    end

    # Build G matrix
    # There will be one row for every :NonNeg and :NonPos cone
    # and an additional row for every variable in a :SOC cone
    # Or in other words, everything that isn't a :Free or :Zero
    # gets a row in G and h
    num_G_row = 0
    for j = 1:num_vars
        idxcone[j] == :Free && continue
        idxcone[j] == :Zero && continue
        num_G_row += 1
    end
    ecos_G = zeros(num_G_row,num_vars)
    ecos_h = zeros(num_G_row)

    # First, handle the :NonNeg, :NonPos cases
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
    # Now handle the SOCs
    # The MPB unput form is basically just says a vector of
    # variables (y,x) lives in the SOC  || x || <= y
    # ECOS wants somethings in the form h - Gx in Q so we
    # will prove 0 - Ix \in Q
    num_SOC_cones = 0
    SOC_conedims  = Int[]
    for (cone, idxs) in cones
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
    @assert G_row == num_G_row + 1

    # Store in the ECOS structure
    m.nvar          = num_vars
    m.nineq         = num_G_row
    m.neq           = length(ecos_b)
    m.npos          = num_pos_orth
    m.ncones        = num_SOC_cones
    m.conedims      = SOC_conedims
    m.G             = ecos_G
    m.A             = ecos_A
    m.c             = ecos_c
    m.orig_sense    = :Min
    m.h             = ecos_h
    m.b             = ecos_b
    m.fwd_map       = fwd_map
end