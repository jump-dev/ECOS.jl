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

@testset "Unit" begin
    MOIT.unittest(MOIB.SplitInterval{Float64}(optimizer), config,
                  [# Quadratic functions are not supported
                   "solve_qcp_edge_cases", "solve_qp_edge_cases",
                   # Integer and ZeroOne sets are not supported
                   "solve_integer_edge_cases", "solve_objbound_edge_cases"])
end

@testset "Continuous linear problems" begin
    MOIT.contlineartest(MOIB.SplitInterval{Float64}(optimizer), config)
end

@testset "Continuous conic problems" begin
    exclude = ["sdp", "rootdet", "logdet"]
    MOIT.contconictest(MOIB.GeoMean{Float64}(MOIB.RSOC{Float64}(optimizer)),
                       config, exclude)
end

@testset "SolverName" begin
    @test MOI.get(optimizer, MOI.SolverName()) == "ECOS"
end

@testset "supports_allocate_load" begin
    @test MOIU.supports_allocate_load(optimizer.optimizer, false)
    @test !MOIU.supports_allocate_load(optimizer.optimizer, true)
end
