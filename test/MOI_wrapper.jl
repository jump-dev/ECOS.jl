using Compat
using Compat.Test

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

import ECOS
const optimizer = ECOS.Optimizer(verbose=false)

@testset "SolverName" begin
    @test MOI.get(optimizer, MOI.SolverName()) == "ECOS"
end

@testset "supports_allocate_load" begin
    @test MOIU.supports_allocate_load(optimizer, false)
    @test !MOIU.supports_allocate_load(optimizer, true)
end

MOIU.@model(ECOSModelData,
            (), (), # No scalar functions
            (MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives, MOI.SecondOrderCone,
             MOI.ExponentialCone),
            (),
            (), (), # No scalar sets
            (), (MOI.VectorAffineFunction,))
# UniversalFallback is needed for starting values, even if they are ignored by ECOS
const cache = MOIU.UniversalFallback(ECOSModelData{Float64}())
const cached = MOIU.CachingOptimizer(cache, optimizer)

# Essential bridges that are needed for all tests
const bridged = MOIB.full_bridge_optimizer(cached, Float64)

# SOC2 requires 1e-4
const config = MOIT.TestConfig(atol=1e-4, rtol=1e-4)

@testset "Unit" begin
    MOIT.unittest(bridged,
                  config,
                  [# Need https://github.com/JuliaOpt/MathOptInterface.jl/issues/529
                   "solve_qp_edge_cases",
                   # Integer and ZeroOne sets are not supported
                   "solve_integer_edge_cases", "solve_objbound_edge_cases"])
end

@testset "Continuous linear problems" begin
    MOIT.contlineartest(bridged, config)
end

@testset "Continuous quadratic problems" begin
    MOIT.qcptest(bridged, config)
end

@testset "Continuous conic problems" begin
    exclude = ["sdp", "rootdet", "logdet"]
    MOIT.contconictest(bridged, config, exclude)
end
