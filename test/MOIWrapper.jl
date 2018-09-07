using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIB = MOI.Bridges

const MOIU = MOI.Utilities
MOIU.@model(ECOSModelData,
            (),
            (MOI.EqualTo, MOI.GreaterThan, MOI.LessThan),
            (MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives, MOI.SecondOrderCone,
             MOI.ExponentialCone),
            (),
            (MOI.SingleVariable,),
            (MOI.ScalarAffineFunction,),
            (MOI.VectorOfVariables,),
            (MOI.VectorAffineFunction,))
const optimizer = MOIU.CachingOptimizer(ECOSModelData{Float64}(), ECOS.Optimizer(verbose=false))

# SOC2 requires 1e-4
const config = MOIT.TestConfig(atol=1e-4, rtol=1e-4)

@testset "Continuous linear problems" begin
    MOIT.contlineartest(MOIB.SplitInterval{Float64}(optimizer), config)
end

@testset "Continuous conic problems" begin
    exclude = ["sdp", "rootdet", "logdet"]
    @static if Compat.Sys.iswindows()
        # Test exp3 fails  on Windows 32 and 64 bits because the windows
        # binaries are out of date just like EXP3 fails with the MPB wrapper
        # See https://github.com/JuliaOpt/ECOS.jl/issues/47
        push!(exclude, "exp3")
    end

    MOIT.contconictest(MOIB.GeoMean{Float64}(MOIB.RSOC{Float64}(optimizer)),
                       config, exclude)
end
