using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIB = MOI.Bridges

const MOIU = MOI.Utilities
MOIU.@model ECOSModelData () (EqualTo, GreaterThan, LessThan) (Zeros, Nonnegatives, Nonpositives, SecondOrderCone, ExponentialCone) () (SingleVariable,) (ScalarAffineFunction,) (VectorOfVariables,) (VectorAffineFunction,)
const optimizer = MOIU.CachingOptimizer(ECOSModelData{Float64}(), ECOS.Optimizer(verbose=false))

# SOC2 requires 1e-4
const config = MOIT.TestConfig(atol=1e-4, rtol=1e-4)

@testset "Continuous linear problems" begin
    MOIT.contlineartest(MOIB.SplitInterval{Float64}(optimizer), config)
end

@testset "Continuous conic problems" begin
    exclude = ["sdp", "rootdet", "logdet"]
    @static if Compat.Sys.iswindows()
        # Test fails on Windows 32 and 64 bits with
        # Expression: (-(y[2]) * log(-(y[2]) / y[4]) + y[2]) - y[3] ≤ tol
        # Evaluated: 0.39942722775671957 ≤ 1.0e-6
        # We do not know the reason yet
        push!(exclude, "exp3")
    end
    MOIT.contconictest(MOIB.GeoMean{Float64}(MOIB.RSOC{Float64}(optimizer)),
                       config, exclude)
end
