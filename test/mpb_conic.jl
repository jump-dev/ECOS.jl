#############################################################################
# ECOS.jl
# Wrapper around the ECOS solver https://github.com/ifa-ethz/ecos
# See http://github.com/JuliaOpt/ECOS.jl
#############################################################################
# test/mpb.jl
# Test the MathProgBase.jl interface for the ECOS.jl solver wrapper
#############################################################################

using Base.Test
import MathProgBase
using ECOS

s = ECOSSolver()

# Problem 1
# min -3x - 2y - 4z
# st    x +  y +  z <= 3
#            y +  z <= 2
#       x>=0 y>=0 z>=0
# Opt solution = -11
# x = 1, y = 0, z = 2
m = MathProgBase.model(s)
ECOS.loadconicproblem!(m, [-3.0, -2.0, -4.0],
                     [ 1.0   1.0   1.0;
                       0.0   1.0   1.0],
                     [ 3.0,  2.0],
                     [:NonNeg, :NonNeg, :NonNeg])
MathProgBase.optimize!(m)
@test MathProgBase.status(m) == :Optimal
@test_approx_eq_eps MathProgBase.getobjval(m) -11 1e-6
@test_approx_eq_eps MathProgBase.getsolution(m)[1] 1.0 1e-6
@test_approx_eq_eps MathProgBase.getsolution(m)[2] 0.0 1e-6
@test_approx_eq_eps MathProgBase.getsolution(m)[3] 2.0 1e-6
