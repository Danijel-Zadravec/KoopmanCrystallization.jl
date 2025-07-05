using Roots

struct BackoffOptimization
    determ_problem::DeterministicProblem
    koop_problem::KoopmanProblem
    determ_opt_sol
    koop_solution
    opt_controls::Vector{Float64}
    indicator_values
    offset_values
    objective_values

end


function BackoffOptimization(determ_problem::DeterministicProblem, unc_params;
    target_violation=0.05, tolerance=0.025, max_iterations=10)

    println("Starting BackoffOptimization with Root Finding...")
    println("Target violation rate: $target_violation ± $tolerance")
    println("Max iterations: $max_iterations")

    # Storage for history
    indicator_values = []
    offset_values = []
    objective_values = []

    # Helper function: Evaluate violation rate for a given offset
    function evaluate_violation(offset_val::Float64)
        println("\n  Evaluating offset: $offset_val")

        # Create fresh copies to avoid parameter contamination
        determ_prob_opt = deepcopy(determ_problem)
        determ_prob_unc = deepcopy(determ_problem)

        # Set the backoff value (not modify, to avoid accumulation)
        determ_prob_opt, _ = set_backoff(determ_prob_opt, offset_val)
        determ_prob_unc, _ = set_backoff(determ_prob_unc, offset_val)

        # Deterministic optimization
        determ_opt_prob = DetermOptimization(determ_prob_opt; penalty=10000000.0)
        determ_opt_sol = OptimizationSolution(determ_opt_prob; maxiters=30)

        opt_controls = determ_opt_sol.optimal_controls
        obj_val = determ_opt_sol.objective_value

        # Update controls in uncertainty problem
        determ_prob_unc = modify_controls(determ_prob_unc, opt_controls)

        # Create and solve Koopman problem
        koop_prob = KoopmanProblem(determ_prob_unc, unc_params)
        sol_koop = solve_koop(koop_prob, ireltol=2e-6, quadalg=HCubatureJL(), maxiters=100000) # vrati na 10

        koop_results = sol_koop.u
        if length(koop_results) >= 2
            bound_viol_perc = koop_results[2]
        else
            @warn "Unexpected Koopman results format"
            return NaN
        end

        # Store history
        push!(indicator_values, bound_viol_perc)
        push!(offset_values, offset_val)
        push!(objective_values, obj_val)

        println("    Violation : $(round(bound_viol_perc, digits=4))")
        println("    Objective: $(round(obj_val, digits=2))")

        return bound_viol_perc
    end

    # Objective function for root finding: f(offset) = violation_rate(offset) - target
    function root_objective(offset_val::Float64)
        violation_rate = evaluate_violation(offset_val)
        if isnan(violation_rate)
            return Inf  # Return large value for failed evaluations
        end
        residual = violation_rate - target_violation
        println("    Residual (violation - target): $(round(residual, digits=4))")
        return residual
    end

    println("\n--- Starting Root Finding Process ---")

    # Initial bracket search
    println("Finding initial bracket...")

    # Start with a reasonable range for offset values
    offset_low = 0.0005
    offset_high = 0.0015
    # Initialize variables that will be used outside try/catch blocks
    optimal_offset = NaN
    final_violation_rate = NaN

    println("\n--- Solving with Brent's Method ---")
    try
        # Use Brent's method for robust root finding
        optimal_offset = find_zero(root_objective, (offset_low, offset_high),
            Roots.Brent(),
            atol=tolerance / 10, rtol=tolerance / 10,
            maxevals=max_iterations)

        println("✓ Root finding converged!")
        println("Optimal offset: $(round(optimal_offset, digits=6))")

        # Final evaluation with optimal offset
        final_violation_rate = evaluate_violation(optimal_offset)

    catch e
        @warn "Root finding failed: $e"
        println("Using best point from evaluations...")

        # Find the offset that gave the closest violation rate to target
        if !isempty(indicator_values)
            deviations = abs.(indicator_values .- target_violation)
            best_idx = argmin(deviations)
            optimal_offset = offset_values[best_idx]
            final_violation_rate = indicator_values[best_idx]
        else
            optimal_offset = 0.0
            final_violation_rate = NaN
        end
    end

    println("\n--- BackoffOptimization Complete ---")
    println("Evaluations performed: ", length(indicator_values))
    println("Final violation rate: $(round(final_violation_rate, digits=4))")
    println("Final offset: $(round(optimal_offset, digits=6))")
    println("Target achieved: $(abs(final_violation_rate - target_violation) <= tolerance)")

    # Create final solution with optimal offset
    final_determ_prob_opt = deepcopy(determ_problem)
    final_determ_prob_unc = deepcopy(determ_problem)

    # Set optimal backoff
    final_determ_prob_opt, _ = set_backoff(final_determ_prob_opt, optimal_offset)
    final_determ_prob_unc, _ = set_backoff(final_determ_prob_unc, optimal_offset)

    # Final optimization
    final_determ_opt_prob = DetermOptimization(final_determ_prob_opt; penalty=10000000.0)
    final_determ_opt_sol = OptimizationSolution(final_determ_opt_prob; maxiters=10)
    final_opt_controls = final_determ_opt_sol.optimal_controls

    # Final uncertainty analysis
    final_determ_prob_unc = modify_controls(final_determ_prob_unc, final_opt_controls)
    final_koop_prob = KoopmanProblem(final_determ_prob_unc, unc_params)
    final_sol_koop = solve_koop(final_koop_prob, ireltol=2e-6, quadalg=HCubatureJL(), maxiters=10000) # vrati na 10

    return BackoffOptimization(
        final_determ_prob_opt,
        final_koop_prob,
        final_determ_opt_sol,
        final_sol_koop,
        final_opt_controls,
        indicator_values,
        offset_values,
        objective_values
    )
end