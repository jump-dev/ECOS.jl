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
    original_sense::Symbol              # Original objective sense
    h::Vector{Float64}                  # RHS for inequality 
    b::Vector{Float64}                  # RHS for equality
    # Post-solve
    solve_stat::Symbol
    obj_val::Float64
    primal_sol::Vector{Float64}
end
ECOSMathProgModel() = ECOSMathProgModel(0,0,0,0,0,
                                        Int[],
                                        spzeros(0,0),
                                        spzeros(0,0),
                                        Float64[], :Min,
                                        Float64[], Float64[],
                                        :NotSolved, 0.0, Float64[])

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
    m.original_sense = sense            # Original objective sense
    m.h         = ineqbnd               # RHS for inequality 
    m.b         = eqbnd                 # RHS for equality
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
    m.obj_val = dot(m.c, m.primal_sol) * (m.original_sense == :Max ? -1 : +1)  
    cleanup(ecos_prob_ptr, 0)
end

status(m::ECOSMathProgModel) = m.solve_stat
getobjval(m::ECOSMathProgModel) = m.obj_val
getsolution(m::ECOSMathProgModel) = m.primal_sol

#############################################################################
# Begin implementation of the MPB conic interface 
# Implements
# - loadconicproblem!
# http://mathprogbasejl.readthedocs.org/en/latest/conic.html

function loadconicproblem!(m::ECOSMathProgModel, c, A, b, cones)
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
    # TODO: Don't assume sorted
    expand_cones = Symbol[]
    for cone_vars in cones
        cone_type, idxs = cone_vars
        if isa(idxs, Int)
            push!(expand_cones, cone_type)
        else 
            # Assume its iterable
            for i in idxs
                push!(expand_cones, cone_type)
            end
        end
    end
    nvar = length(expand_cones)

    # Start with the data provided
    ecos_c = copy(c)
    ecos_A = copy(A)
    ecos_b = copy(b)

    # For all variables in the :Zero cone, fix at 0 with an
    # equality constraint.
    # TODO: Don't even include them
    for j = 1:nvar
        if expand_cones[j] == :Zero
            new_row = zeros(1,nvar)
            new_row[j] = 1.0
            ecos_A = vcat(ecos_A, new_row)
            ecos_b = vcat(ecos_b, 0.0)
        end
    end

    # Build G matrix
    ecos_G = Array(Float64,0,nvar)
    ecos_h = Array(Float64,0)
    # First, handle the :NonNeg, :NonPos cases
    npos = 0
    for j = 1:nvar
        if expand_cones[j] == :NonNeg
            new_row = zeros(1,nvar)
            new_row[j] = -1.0
            ecos_G = vcat(ecos_G, new_row)
            ecos_h = vcat(ecos_h, 0.0)
            npos += 1
        elseif expand_cones[j] == :NonPos
            new_row = zeros(1,nvar)
            new_row[j] = +1.0
            ecos_G = vcat(ecos_G, new_row)
            ecos_h = vcat(ecos_h, 0.0)
            npos += 1
        end
    end
    # Now handle the SOCs
    # Input form is basically just says
    # || x || <= y
    # and we pass into ECOS
    # 0 - Ix \in Q
    # which maps back to the same thing
    ncones = 0
    conedims = Int[]
    for cone_vars in cones
        cone_type, idxs = cone_vars
        cone_type != :SOC && continue
        ncones += 1
        push!(conedims, length(idxs))
        new_rows = zeros(length(idxs),nvar)
        row = 1
        for j in idxs
            new_rows[row,j] = -1.0
            row += 1
            ecos_h = vcat(ecos_h, 0.0)
        end
        ecos_G = vcat(ecos_G, new_rows)
    end

    m.nvar      = nvar
    m.nineq     = length(ecos_h)
    m.neq       = length(ecos_b)
    m.npos      = npos
    m.ncones    = ncones
    m.conedims  = conedims
    m.G         = ecos_G
    m.A         = ecos_A
    m.c         = ecos_c
    m.original_sense = :Min
    m.h         = ecos_h
    m.b         = ecos_b
end