#############################################################################
# ECOS.jl
# Wrapper around the ECOS solver https://github.com/ifa-ethz/ecos
# See http://github.com/jump-dev/ECOS.jl
#############################################################################
# test/MPB_linear.jl
# Test the MathProgBase.jl interface for the ECOS.jl solver wrapper
#############################################################################

using Test
using LinearAlgebra
using MathProgBase
using SparseArrays
using ECOS

solver = ECOSSolver()
objtol = 1e-7
primaltol = 1e-6

# min -x
# s.t. 2x + y <= 1.5
# x,y >= 0
# solution is (0.75,0) with objval -0.75
sol = linprog([-1,0],[2 1],'<',1.5,solver)
@test sol.status == :Optimal
@test isapprox(sol.objval, -0.75, atol=objtol)
@test isapprox(norm(sol.sol - [0.75,0.0]), 0, atol=primaltol)

sol = linprog([-1,0],sparse([2 1]),'<',1.5,solver)
@test sol.status == :Optimal
@test isapprox(sol.objval, -0.75, atol=objtol)
@test isapprox(norm(sol.sol - [0.75,0.0]), 0, atol=primaltol)

# test infeasible problem:
# min x
# s.t. 2x+y <= -1
# x,y >= 0
sol = linprog([1.0,0.0],[2.0 1.0],'<',-1.0,solver)
@test sol.status == :Infeasible

# test unbounded problem:
# min -x-y
# s.t. -x+2y <= 0
# x,y >= 0
sol = linprog([-1.,-1.],[-1. 2.],'<',[0.0],solver)
@test sol.status == :Unbounded

# More complicated through loadproblem!
A = [1.0 1.0 0.0 0.0 0.0;
     0.0 1.0 0.0 0.0 1.0;
     0.0 0.0 1.0 1.0 1.0]
collb = [0.0, 0.0, 0.0, 0.0, 0.0]
colub = [Inf, Inf, Inf, Inf, Inf]
obj   = [3.0, 4.0, 4.0, 9.0, 5.0]
sense = :Max
rowlb = [-Inf, -Inf, -Inf]
rowub = [ 5.0,  3.0,  9.0]
s = ECOSSolver()
m = MathProgBase.LinearQuadraticModel(s)
MathProgBase.loadproblem!(m, A, collb, colub, obj, rowlb, rowub, sense)
MathProgBase.optimize!(m)
@test isapprox(MathProgBase.getobjval(m), 99.0, atol=1e-5)
x = MathProgBase.getsolution(m)
@test isapprox(x[1], 2.0, atol=1e-5)
@test isapprox(x[2], 3.0, atol=1e-5)
@test isapprox(x[3], 0.0, atol=1e-5)
@test isapprox(x[4], 9.0, atol=1e-5)
@test isapprox(x[5], 0.0, atol=1e-5)
