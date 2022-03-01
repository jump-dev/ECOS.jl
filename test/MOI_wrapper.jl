module TestECOS

using Test
import ECOS
import MathOptInterface

const MOI = MathOptInterface

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function test_runtests()
    model = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        MOI.instantiate(ECOS.Optimizer; with_bridge_type = Float64),
    )
    @test model.optimizer.model.model_cache isa
          MOI.Utilities.UniversalFallback{ECOS.OptimizerCache}
    MOI.set(model, MOI.Silent(), true)
    exclude = String[
        # ZerosBridge does not support ConstraintDual. These are tested below in
        # test_runtests_ZerosBridge
        "test_conic_RotatedSecondOrderCone_INFEASIBLE_2",
        "test_conic_linear_VectorOfVariables_2",
        "test_linear_integration",
        "test_quadratic_constraint_GreaterThan",
        "test_quadratic_constraint_LessThan",
    ]
    if Sys.WORD_SIZE == 32
        # These tests fail on x86 Linux, returning ITERATION_LIMIT instead of
        # proving {primal,dual}_INFEASIBLE.
        push!(exclude, "test_conic_linear_INFEASIBLE")
        push!(exclude, "test_solve_TerminationStatus_DUAL_INFEASIBLE")
    end
    MOI.Test.runtests(
        model,
        MOI.Test.Config(
            atol = 1e-3,
            rtol = 1e-3,
            exclude = Any[
                MOI.ConstraintBasisStatus,
                MOI.VariableBasisStatus,
                MOI.ObjectiveBound,
            ],
        ),
        exclude = exclude,
    )
    return
end

function test_runtests_ZerosBridge()
    optimizer = MOI.instantiate(ECOS.Optimizer; with_bridge_type = Float64)
    MOI.Bridges.remove_bridge(
        optimizer,
        MOI.Bridges.Variable.ZerosBridge{Float64},
    )
    model = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        optimizer,
    )
    MOI.Test.runtests(
        model,
        MOI.Test.Config(
            atol = 1e-3,
            rtol = 1e-3,
            exclude = Any[
                MOI.ConstraintBasisStatus,
                MOI.VariableBasisStatus,
                MOI.ObjectiveBound,
            ],
        ),
        include = String[
            # ZerosBridge does not support ConstraintDual
            "test_conic_RotatedSecondOrderCone_INFEASIBLE_2",
            "test_conic_linear_VectorOfVariables_2",
            "test_linear_integration",
            "test_quadratic_constraint_GreaterThan",
            "test_quadratic_constraint_LessThan",
        ],
    )
    return
end

function test_RawOptimizerAttribute()
    model = ECOS.Optimizer()
    MOI.set(model, MOI.RawOptimizerAttribute("abstol"), 1e-5)
    @test MOI.get(model, MOI.RawOptimizerAttribute("abstol")) ≈ 1e-5
    MOI.set(model, MOI.RawOptimizerAttribute("abstol"), 2e-5)
    @test MOI.get(model, MOI.RawOptimizerAttribute("abstol")) == 2e-5
    return
end

function test_iteration_limit()
    v = [5.0, 3.0, 1.0]
    w = [2.0, 1.5, 0.3]
    solver = MOI.OptimizerWithAttributes(ECOS.Optimizer, MOI.Silent() => true)
    model = MOI.instantiate(solver, with_bridge_type = Float64)
    maxit = 1
    MOI.set(model, MOI.RawOptimizerAttribute("maxit"), maxit)
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 3)
    MOI.add_constraint.(model, x, MOI.Interval(0.0, 1.0))
    MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(w, x), 0.0),
        MOI.LessThan(3.0),
    )
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(v, x), 0.0),
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.ITERATION_LIMIT
    @test MOI.get(model, MOI.BarrierIterations()) == maxit
    return
end

function test_empty_problem()
    model = MOI.Utilities.Model{Float64}()
    ecos = ECOS.Optimizer()
    MOI.optimize!(ecos, model)
    @test MOI.get(ecos, MOI.TerminationStatus()) == MOI.INVALID_MODEL
    @test MOI.get(ecos, MOI.PrimalStatus()) == MOI.NO_SOLUTION
    @test MOI.get(ecos, MOI.DualStatus()) == MOI.NO_SOLUTION
    return
end

function test_conic_no_variables()
    model = MOI.Utilities.Model{Float64}()
    ecos = ECOS.Optimizer()
    f = MOI.VectorAffineFunction(
        MOI.VectorAffineTerm{Float64}[],
        [1.0, 0.5, 0.5],
    )
    MOI.add_constraint(model, f, MOI.SecondOrderCone(3))
    MOI.optimize!(ecos, model)
    @test MOI.get(ecos, MOI.TerminationStatus()) == MOI.INVALID_MODEL
    @test MOI.get(ecos, MOI.PrimalStatus()) == MOI.NO_SOLUTION
    @test MOI.get(ecos, MOI.DualStatus()) == MOI.NO_SOLUTION
    return
end

end  # module

TestECOS.runtests()
