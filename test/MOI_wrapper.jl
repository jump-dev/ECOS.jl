using Test

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

import ECOS
const optimizer = ECOS.Optimizer()
MOI.set(optimizer, MOI.Silent(), true)

@testset "SolverName" begin
    @test MOI.get(optimizer, MOI.SolverName()) == "ECOS"
end

@testset "supports_allocate_load" begin
    @test MOIU.supports_allocate_load(optimizer, false)
    @test !MOIU.supports_allocate_load(optimizer, true)
end

# UniversalFallback is needed for starting values, even if they are ignored by ECOS
const cache = MOIU.UniversalFallback(MOIU.Model{Float64}())
const cached = MOIU.CachingOptimizer(cache, optimizer)

const bridged = MOIB.full_bridge_optimizer(cached, Float64)

# SOC2 requires 1e-4
const config = MOIT.TestConfig(atol=1e-4, rtol=1e-4)

@testset "Unit" begin
    MOIT.unittest(bridged,
                  config,
                  [
        # `NumberOfThreads` not supported.
        "number_threads",
        # `TimeLimitSec` not supported.
        "time_limit_sec",
        # Need https://github.com/jump-dev/MathOptInterface.jl/issues/529
        "solve_qp_edge_cases",
        # Integer and ZeroOne sets are not supported
        "solve_integer_edge_cases", "solve_objbound_edge_cases",
        "solve_zero_one_with_bounds_1",
        "solve_zero_one_with_bounds_2",
        "solve_zero_one_with_bounds_3"])
end

@testset "Continuous linear problems" begin
    MOIT.contlineartest(bridged, config)
end

@testset "Continuous quadratic problems" begin
    MOIT.qcptest(bridged, config)
end

@testset "Continuous conic problems" begin
    exclude = ["dualexp", "pow", "dualpow", "sdp", "rootdet", "logdet", "normnuc", "normspec"]
    MOIT.contconictest(bridged, config, exclude)
end

@testset "MOI.RawParameter" begin
    model = ECOS.Optimizer()
    # TODO: remove symbol cases when deprecation is removed.
    MOI.set(model, MOI.RawParameter(:abstol), 1e-5)
    @test MOI.get(model, MOI.RawParameter(:abstol)) ≈ 1e-5
    @test MOI.get(model, MOI.RawParameter("abstol")) ≈ 1e-5
    MOI.set(model, MOI.RawParameter("abstol"), 2e-5)
    @test MOI.get(model, MOI.RawParameter(:abstol)) ≈ 2e-5
    @test MOI.get(model, MOI.RawParameter("abstol")) == 2e-5
end

@testset "Iteration Limit" begin
    # Problem data
    v = [5.0, 3.0, 1.0]
    w = [2.0, 1.5, 0.3]

    MOI.empty!(bridged)

    MOI.set(bridged, MOI.RawParameter("maxit"), 1)
    MOI.set(bridged, MOI.Silent(), true)

    x = MOI.add_variables(bridged, 3)
    for xj in x
        MOI.add_constraint(bridged, MOI.SingleVariable(xj), MOI.Interval(0.0, 1.0))
    end

    # Constraint
    MOI.add_constraint(bridged, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(w, x), 0.0), MOI.LessThan(3.0))

    # Set the objective
    MOI.set(bridged,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(v, x), 0.0)
    )
    MOI.set(bridged, MOI.ObjectiveSense(), MOI.MAX_SENSE)

    MOI.optimize!(bridged)

    @test MOI.get(bridged, MOI.TerminationStatus()) == MOI.ITERATION_LIMIT
end
