#############################################################################
# ECOS.jl
# Wrapper around the ECOS solver https://github.com/ifa-ethz/ecos
# See http://github.com/JuliaOpt/ECOS.jl
#############################################################################
# test/MPBWrapper.jl
# Test the MathProgBase.jl interface for the ECOS.jl solver wrapper
#############################################################################

@testset "Test the MPB wrapper with linprog" begin
    include("mpb_linear.jl")
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
    @static if !Compat.Sys.iswindows()
        # Test EXP3 fails  on Windows 32 and 64 bits because the windows
        # binaries are out of date. There failure is:
        # Expression: (-(y[2]) * log(-(y[2]) / y[4]) + y[2]) - y[3] ≤ tol
        # Evaluated: 0.39942722775671957 ≤ 1.0e-6
        # See https://github.com/JuliaOpt/ECOS.jl/issues/47
        conicEXPtest(ECOS.ECOSSolver(), duals=true)
    end

    # See https://github.com/JuliaOpt/MathProgBase.jl/pull/218
    if VERSION < v"0.7-"
        include(joinpath(MPB_test_path, "quadprog.jl"))
        socptest(ECOS.ECOSSolver())
    end
end
