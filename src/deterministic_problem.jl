using Plots
using SciMLStructures
using SymbolicIndexingInterface
using CSV
using DataFrames
using Printf

struct DeterministicProblem
    process_model::ProcessModel
    initial_states
    initial_controls
    param_values
    problem::ODEProblem
end


function create_init_controls(n_steps::Int, start, stop)
    # Create initial control values
    controls = [start + (stop - start) * i / (n_steps) for i in 0:(n_steps)]
    return controls
end


function save_params_states(model::ProcessModel, states, parameters, controls)
    # Helper function to extract symbol name from various types
    function get_symbol_name(var)
        try
            return nameof(var)
        catch
            # Fallback: convert to string and extract name
            var_str = string(var)
            if occursin("(", var_str)
                # Extract name before parentheses: "cp(t)" -> "cp"
                return Symbol(split(var_str, "(")[1])
            else
                return Symbol(var_str)
            end
        end
    end

    # Create initial conditions array
    state_dict = Dict(states)
    base_initial_conditions = []
    for var in model.symbols.variables
        var_symbol = get_symbol_name(var)
        if haskey(state_dict, var_symbol)
            push!(base_initial_conditions, var => state_dict[var_symbol])
        else
            @warn "State variable $var_symbol not found in initial conditions, skipping"
        end
    end

    # Create parameters array
    param_dict = Dict(parameters)
    base_param_values = []
    for param in model.symbols.parameters
        param_symbol = get_symbol_name(param)
        if haskey(param_dict, param_symbol)
            push!(base_param_values, param => param_dict[param_symbol])
        else
            @warn "Parameter $param_symbol not found in parameter values, setting to 0.0"
            push!(base_param_values, param => 0.0)
        end
    end

    # Create control initial values array with individual parameter mapping
    control_initial_values = []

    # Process u0 controls
    for (i, u_symbol) in enumerate(model.symbols.controls)
        if i <= length(controls)
            push!(control_initial_values, u_symbol => controls[i])
        else
            @warn "u control value for step $i not found, setting to 0.000001"
            push!(control_initial_values, u0_symbol => 0.000001)
        end
    end
    return base_initial_conditions, base_param_values, control_initial_values
end


function DeterministicProblem(process_model::ProcessModel, initial_states, parameters, controls)
    base_initial_conditions, base_param_values, control_initial_values = save_params_states(process_model, initial_states, parameters, controls)

    all_parameters = [base_param_values..., control_initial_values...]

    problem = ODEProblem(process_model.model, initial_states, (0, 1800.0), all_parameters)
    return DeterministicProblem(process_model, base_initial_conditions, control_initial_values, base_param_values, problem)
end


function solve_controls(prob::DeterministicProblem, controls; abstol=1e-10, reltol=1e-8)
    """
    Solve the deterministic problem with given tolerances
    """
    odeprob = deepcopy(prob.problem)
    u_params = prob.process_model.symbols.controls
    # Create a setter function for updating parameters
    setter_fn! = setp(odeprob, u_params)
    setter_fn!(odeprob, controls)
    return solve(odeprob, Rodas5P(), abstol=abstol, reltol=reltol), odeprob

end

function modify_controls(prob::DeterministicProblem, new_controls)
    """
    Modify the controls of the deterministic problem
    """
    odeprob = prob.problem
    u_params = prob.process_model.symbols.controls
    setter_fn! = setp(odeprob, u_params)
    setter_fn!(odeprob, new_controls)
    return prob
end

function modify_backoff(prob::DeterministicProblem, offset_corr)
    """
    Modify the backoff parameter in the deterministic problem
    """
    sys = prob.process_model.model
    idx_b = parameter_index(sys, :h_offset)
    offset_prv = getindex(prob.problem.p, idx_b)
    setindex!(prob.problem.p, offset_prv + offset_corr, idx_b)
    return prob, offset_prv + offset_corr
end
function set_backoff(prob::DeterministicProblem, offset_new)
    """
    Modify the backoff parameter in the deterministic problem
    """
    sys = prob.process_model.model
    idx_b = parameter_index(sys, :h_offset)
    setindex!(prob.problem.p, offset_new, idx_b)
    return prob, offset_new
end


function plot_states(sol::ODESolution, determ_problem::DeterministicProblem; plot_title="Crystallization Process States")
    """
    Plot crystallization states on separate axes - each y variable individually
    """

    # Get process model for accessing variable symbols
    model = determ_problem.process_model

    # Time vector for plotting
    tspan = sol.prob.tspan
    t_plot = range(tspan[1], stop=tspan[2], length=1000)
    
    # Convert time to minutes for x-axis
    t_minutes = t_plot ./ 60.0

    # Create individual plots for each y variable
    plots_array = []
    h_offset = sol.ps[:h_offset]
    
    # Plot 1: y1 (Solute Concentration) with saturation limits
    p1 = plot(t_minutes, [sol(t, idxs=:y1) for t in t_plot],
        label="Concentration", xlabel="Time [min]", ylabel="Concentration [g/g solvent]",
        title="Solute Concentration", linewidth=2, color=:blue)
    plot!(p1, t_minutes, [sol(t, idxs=:Cs) for t in t_plot],
        label="Saturation Cs", linestyle=:dash, linewidth=2, color=:green)
    plot!(p1, t_minutes, [sol(t, idxs=:Cm) for t in t_plot],
        label="Metastable Cm", linestyle=:dot, linewidth=2, color=:red)
    plot!(p1, t_minutes, [sol(t, idxs=:Cm) for t in t_plot] .- h_offset,
        label="Metastable Cm corrected", linestyle=:dot, linewidth=2, color=:cyan)
    push!(plots_array, p1)

    # Plot 2: y2 (Zeroth moment of seeded crystals)
    p2 = plot(t_minutes, [sol(t, idxs=:y2) for t in t_plot],
        label="Seeded crystals", xlabel="Time [min]", ylabel="Crystal Count [#/g solvent]",
        title="Seeded Crystal Count", linewidth=2, color=:purple,
        yformatter=:scientific, left_margin=15Plots.mm)
    push!(plots_array, p2)

    # Plot 3: y3 (First moment of seeded crystals)
    p3 = plot(t_minutes, [sol(t, idxs=:y3) for t in t_plot],
        label="Seeded crystals", xlabel="Time [min]", ylabel="1st Moment [μm/g solvent]",
        title="Seeded 1st Moment", linewidth=2, color=:orange,
        yformatter=:scientific, left_margin=15Plots.mm)
    push!(plots_array, p3)

    # Plot 4: y4 (Second moment of seeded crystals)
    p4 = plot(t_minutes, [sol(t, idxs=:y4) for t in t_plot],
        label="Seeded crystals", xlabel="Time [min]", ylabel="2nd Moment [μm²/g solvent]",
        title="Seeded 2nd Moment", linewidth=2, color=:cyan,
        yformatter=:scientific, left_margin=15Plots.mm)
    push!(plots_array, p4)

    # Plot 5: y5 (Third moment of seeded crystals)
    p5 = plot(t_minutes, [sol(t, idxs=:y5) for t in t_plot],
        label="Seeded crystals", xlabel="Time [min]", ylabel="3rd Moment [μm³/g solvent]",
        title="Seeded 3rd Moment", linewidth=2, color=:darkgreen,
        yformatter=:scientific, left_margin=15Plots.mm)
    push!(plots_array, p5)

    # Plot 6: y6 (Zeroth moment of nucleated crystals)
    p6 = plot(t_minutes, [sol(t, idxs=:y6) for t in t_plot],
        label="Nucleated crystals", xlabel="Time [min]", ylabel="Crystal Count [#/g solvent]",
        title="Nucleated Crystal Count", linewidth=2, color=:navy)
    push!(plots_array, p6)

    # Plot 7: y7 (First moment of nucleated crystals)
    p7 = plot(t_minutes, [sol(t, idxs=:y7) for t in t_plot],
        label="Nucleated crystals", xlabel="Time [min]", ylabel="1st Moment [μm/g solvent]",
        title="Nucleated 1st Moment", linewidth=2, color=:brown,
        yformatter=:scientific, left_margin=15Plots.mm)
    push!(plots_array, p7)

    # Plot 8: y8 (Second moment of nucleated crystals)
    p8 = plot(t_minutes, [sol(t, idxs=:y8) for t in t_plot],
        label="Nucleated crystals", xlabel="Time [min]", ylabel="2nd Moment [μm²/g solvent]",
        title="Nucleated 2nd Moment", linewidth=2, color=:magenta,
        yformatter=:scientific, left_margin=15Plots.mm)
    push!(plots_array, p8)

    # Plot 9: y9 (Third moment of nucleated crystals)
    p9 = plot(t_minutes, [sol(t, idxs=:y9) for t in t_plot],
        label="Nucleated crystals", xlabel="Time [min]", ylabel="3rd Moment [μm³/g solvent]",
        title="Nucleated 3rd Moment", linewidth=2, color=:darkred,
        yformatter=:scientific, left_margin=15Plots.mm)
    push!(plots_array, p9)

    # Plot 10: Temperature profile
    p10 = plot(t_minutes, [sol(t, idxs=:T) for t in t_plot],
        label="Temperature", xlabel="Time [min]", ylabel="Temperature [K]",
        title="Temperature Profile", linewidth=2, color=:red)
    push!(plots_array, p10)

    # Plot 11: C_high (Over-concentration violation)
    p11 = plot(t_minutes, [sol(t, idxs=:C_high) for t in t_plot],
        label="Over-concentration", xlabel="Time [min]", ylabel="Over-concentration",
        title="Over-concentration Violation", linewidth=2, color=:orange)
    push!(plots_array, p11)

    # Plot 12: C_low (Under-concentration violation)
    p12 = plot(t_minutes, [sol(t, idxs=:C_low) for t in t_plot],
        label="Under-concentration", xlabel="Time [min]", ylabel="Under-concentration",
        title="Under-concentration Violation", linewidth=2, color=:lightblue)
    push!(plots_array, p12)

    # Create combined plot layout (6x2 grid for 12 plots) with main title
    combined_plot = plot(plots_array..., layout=(6, 2), size=(1200, 1400), 
                        plot_title=plot_title, titlefontsize=16,
                        left_margin=5Plots.mm, right_margin=5Plots.mm)

    # Display combined plot
    display(combined_plot)

    # Calculate and display objective
    t_end = tspan[2]
    y5_final = sol(t_end, idxs=:y5)
    y9_final = sol(t_end, idxs=:y9)
    objective_value = -(y5_final - y9_final)

    println("\n=== Crystallization Results ===")
    println("y1 (final concentration): ", round(sol(t_end, idxs=:y1), digits=6))
    println("y2 (seeded count): ", round(sol(t_end, idxs=:y2), digits=2))
    println("y3 (seeded 1st moment): ", round(sol(t_end, idxs=:y3), digits=2))
    println("y4 (seeded 2nd moment): ", round(sol(t_end, idxs=:y4), digits=2))
    println("y5 (seeded 3rd moment): ", y5_final)
    println("y6 (nucleated count): ", round(sol(t_end, idxs=:y6), digits=2))
    println("y7 (nucleated 1st moment): ", round(sol(t_end, idxs=:y7), digits=2))
    println("y8 (nucleated 2nd moment): ", round(sol(t_end, idxs=:y8), digits=2))
    println("y9 (nucleated 3rd moment): ", y9_final)
    println("Final temperature: ", round(sol(t_end, idxs=:T), digits=2), " K")
    println("C_high (final): ", round(sol(t_end, idxs=:C_high), digits=6))
    println("C_low (final): ", round(sol(t_end, idxs=:C_low), digits=6))

    println("Objective -(y5-y9): ", objective_value)

    return combined_plot
end
