module TestECOS

using Test
using MathOptInterface
import ECOS

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
    solver = MOI.OptimizerWithAttributes(ECOS.Optimizer, MOI.Silent() => true)
    model = MOI.instantiate(solver, with_bridge_type = Float64)
    MOI.Test.runtests(
        model,
        MOI.Test.Config(
            atol = 1e-3,
            rtol = 1e-3,
            exclude = Any[
                MOI.ConstraintBasisStatus,
                MOI.VariableBasisStatus,
                MOI.ConstraintName,
                MOI.VariableName,
                MOI.ObjectiveBound,
            ],
        ),
        exclude = String[
            # Expected test failures:
            #   Problem is a nonconvex QP
            "test_basic_ScalarQuadraticFunction_EqualTo",
            "test_basic_ScalarQuadraticFunction_GreaterThan",
            "test_basic_ScalarQuadraticFunction_Interval",
            "test_basic_VectorQuadraticFunction_",
            "test_quadratic_SecondOrderCone_basic",
            "test_quadratic_nonconvex_",
            #   MathOptInterface.jl issue #1431
            "test_model_LowerBoundAlreadySet",
            "test_model_UpperBoundAlreadySet",
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
end

function test_iteration_limit()
    # Problem data
    v = [5.0, 3.0, 1.0]
    w = [2.0, 1.5, 0.3]


    solver = MOI.OptimizerWithAttributes(ECOS.Optimizer, MOI.Silent() => true)
    model = MOI.instantiate(solver, with_bridge_type = Float64)

    maxit = 1
    MOI.set(model, MOI.RawOptimizerAttribute("maxit"), maxit)
    MOI.set(model, MOI.Silent(), true)

    x = MOI.add_variables(model, 3)
    for xj in x
        MOI.add_constraint(
            model,
            MOI.SingleVariable(xj),
            MOI.Interval(0.0, 1.0),
        )
    end

    # Constraint
    MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(w, x), 0.0),
        MOI.LessThan(3.0),
    )

    # Set the objective
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(v, x), 0.0),
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)

    MOI.optimize!(model)

    @test MOI.get(model, MOI.TerminationStatus()) == MOI.ITERATION_LIMIT
    @test MOI.get(model, MOI.BarrierIterations()) == maxit
end

end  # module

TestECOS.runtests()
