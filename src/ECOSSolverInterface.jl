export ECOSSolver

type ECOSMathProgModel <: AbstractMathProgModel
    nvar::Int
    nineq::Int
    neq::Int
    npos::Int
    ncones::Int
    conedims::Vector{Int}
    G::SparseMatrixCSC{Float64,Int}
    A::SparseMatrixCSC{Float64,Int}
    c::Vector{Float64}
    sense::Symbol
    h::Vector{Float64}
    b::Vector{Float64}
    solve_stat::Symbol
end

ECOSMathProgModel() = ECOSMathProgModel(0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        Int[],
                                        spzeros(0),
                                        spzeros(0),
                                        Float64[],
                                        Float64[],
                                        :NotSolved)

immutable ECOSSolver <: AbstractMathProgSolver

function loadproblem!(m::ECOSMathProgModel, A, collb, colub, obj, rowlb, rowub, sense)
    (nvar = length(collb)) == length(colub) || error("Unequal lengths for column bounds")
    (nrow = length(rowlb)) == length(rowub) || "Unequal lengths for row bounds"
    sense == :Max && (obj *= -1)
    eqidx = Int[]
    ineqidx = Int[]
    eqbnd = Float64[]
    ineqbnd = Float64[]
    for it in 1:nrow
        if rowlb[it] == rowub[it]
            push!(eqidx, it)
            push!(eqbnd, rowlb[it])
        elseif rowlb[it] != -Inf && rowub != Inf
            error("Not yet support for ranged constraints")
        elseif rowlb[it] == -Inf
            push!(ineqidx, it)
            push!(ineqbnd, rowub[it])
        else
            push!(ineqidx, it)
            push!(ineqbnd, -rowlb[it])
            A[it,:] *= -1 # flip signs so we have Ax<=b
        end
    end
    ncones = 0
    conedims = Int[]
    m = ECOSMathProgModel(nvar,
                          length(ineqidx),
                          length(eqidx),
                          length(ineqidx), # not sure about this...
                          0,
                          Int[],
                          sparse(A[eqidx,:]), # much less efficient than desirable
                          sparse(A[ineqidx,:]),
                          obj,
                          sense,
                          ineqbnd,
                          eqbnd,
                          :NotSolved)
    nothing
end

function optimize!(m::ECOSMathProgModel)
    A = CSCtoCCS(m.A)
    G = CSCtoCCS(m.G)
    problem = setup(m.nvar,
                    m.nineq,
                    m.neq,
                    m.npos,
                    m.ncones,
                    m.conedims,
                    G,
                    A,
                    m.c,
                    m.h,
                    m.b)

    m.solve_stat = solve(problem)
    cleanup(problem, 1)
    nothing
end