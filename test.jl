using JuMP
using Ipopt
using MathOptEnzyme
using Test

function _run_solution_benchmark(model, optimizer, backend)
    set_optimizer(model, optimizer)
    optimize!(model; _differentiation_backend = backend)
    println("Timing statistics")
    nlp_block = MOI.get(model, MOI.NLPBlock())
    println(" - ", nlp_block.evaluator.eval_constraint_timer)
    println(" - ", nlp_block.evaluator.eval_constraint_jacobian_timer)
    println(" - ", nlp_block.evaluator.eval_hessian_lagrangian_timer)
    return value.(all_variables(model))
end

function run_solution_benchmark(
    model::JuMP.Model,
    optimizer;
    atol::Float64 = 1e-6,
    rtol::Float64 = 0.0,
)
    @info "Solving with SparseReverseMode"
    moi_solution = _run_solution_benchmark(
        model,
        optimizer,
        MOI.Nonlinear.SparseReverseMode(),
    )
    @info "Solving with serial MathOptSymbolicAD"
    serial_solution = _run_solution_benchmark(
        model,
        optimizer,
        MathOptEnzyme.EnzymeBackend(),
    )
    @info "Validating solutions"
    println("Errors = ", extrema(moi_solution .- serial_solution))
    Test.@test ≈(moi_solution, serial_solution; atol = atol, rtol = rtol)
    return
end

function run_clnlbeam_benchmark(; N::Int)
    h = 1 / N
    model = Model()
    @variable(model, -1 <= t[1:(N+1)] <= 1)
    @variable(model, -0.05 <= x[1:(N+1)] <= 0.05)
    @variable(model, u[1:(N+1)])
    @NLobjective(
        model,
        Min,
        sum(
            h / 2 * (u[i+1]^2 + u[i]^2) +
            350 * h / 2 * (cos(t[i+1]) + cos(t[i])) for i in 1:N
        ),
    )
    for i in 1:N
        @NLconstraint(model, x[i+1] - x[i] == h / 2 * (sin(t[i+1]) + sin(t[i])))
        @constraint(model, t[i+1] - t[i] == h / 2 * (u[i+1] - u[i]))
    end
    # Different rtol because the points we pick can be quite bad
    # run_unit_benchmark(model; rtol = 1e-6)
    run_solution_benchmark(model, Ipopt.Optimizer)
    return
end

run_clnlbeam_benchmark(; N = 10)