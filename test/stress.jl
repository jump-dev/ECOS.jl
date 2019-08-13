# Generates a random cone problem to test the MPB interface perf
using SparseArrays

rows = 10000
cols = 10000

srand(1988)
A = sprand(rows, cols, 0.01)
b = rand(rows)
c = rand(cols)


coneset = [:NonNeg, :NonPos, :SOC]
var_cones = {}
vars_coned = 0
while true
    (vars_coned >= cols) && break
    if vars_coned == cols - 1
        push!(var_cones, (:NonNeg, cols))
        break
    end
    cone = coneset[rand(1:3)]
    n_var = rand(2:5)
    end_var = min(vars_coned+n_var, cols)
    push!(var_cones, (cone, vars_coned+1:end_var) )
    vars_coned += n_var
end

coneset = [:Zero, :NonNeg, :NonPos, :SOC]
row_cones = {}
rows_coned = 0
while true
    (rows_coned >= rows) && break
    if rows_coned == rows - 1
        push!(row_cones, (:NonNeg, rows))
        break
    end
    cone = coneset[rand(1:4)]
    n_row = rand(2:5)
    end_row = min(rows_coned+n_row, rows)
    push!(row_cones, (cone, rows_coned+1:end_row) )
    rows_coned += n_row
end


using ECOS
m = ECOS.ECOSMathProgModel()

println("First run")
@time    ECOS.loadconicproblem!(m, c, A, b, row_cones, var_cones)
Base.gc()

println("Second run")
@profile ECOS.loadconicproblem!(m, c, A, b, row_cones, var_cones)
Profile.print(format=:flat)
Base.gc()

println("Third run")
@time    ECOS.loadconicproblem!(m, c, A, b, row_cones, var_cones)
