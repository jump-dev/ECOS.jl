#############################################################################
# ECOS.jl
# Wrapper around the ECOS solver https://github.com/ifa-ethz/ecos
# See http://github.com/JuliaOpt/ECOS.jl
#############################################################################
# test/MPBWrapper.jl
# Test the MathProgBase.jl interface for the ECOS.jl solver wrapper
#############################################################################

@testset "Test the MPB wrapper with linprog" begin
    include("MPB_linear.jl")
end

import ECOS

import MathProgBase
@static if VERSION >= v"0.7-"
    const MPB_test_path = joinpath(dirname(pathof(MathProgBase)), "..", "test")
else
    const MPB_test_path = joinpath(Pkg.dir("MathProgBase"), "test")
end

@testset "Run the conic interface test from MathProgBase.jl" begin
    include(joinpath(MPB_test_path, "conicinterface.jl"))
    coniclineartest(ECOS.ECOSSolver(), duals=true)
    conicSOCtest(ECOS.ECOSSolver(), duals=true)
    conicEXPtest(ECOS.ECOSSolver(), duals=true)
    include(joinpath(MPB_test_path, "quadprog.jl"))
    socptest(ECOS.ECOSSolver())
end
