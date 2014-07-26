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
# Begin implementation of the MPB interface
# Implements
# - model
# - loadproblem!
# - optimize!
# - status
# http://mathprogbasejl.readthedocs.org/en/latest/lowlevel.html
# http://mathprogbasejl.readthedocs.org/en/latest/conic.html

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
            error("Not yet support for ranged constraints")
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