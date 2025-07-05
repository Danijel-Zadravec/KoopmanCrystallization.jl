using KoopmanCrystallization
using SciMLExpectations
using Statistics
n_steps = 30

# Initial conditions
initial_conditions = [
    :y1 => 0.1743,         # Initial solute concentration [g/g solvent]
    :y2 => 66.66,         # Initial zeroth moment of seeded crystals [#/m³]
    :y3 => 1.83e4,         # Initial first moment of seeded crystals [μm/m³]
    :y4 => 5.05e6,          # Initial second moment of seeded crystals [μm²/m³]
    :y5 => 1.93e9,          # Initial third moment of seeded crystals [μm³/m³]
    :y6 => 0.867,          # Initial zeroth moment of nucleated crystals
    :y7 => 0.0,          # Initial first moment of nucleated crystals
    :y8 => 0.0,          # Initial second moment of nucleated crystals
    :y9 => 0.0,           # Initial third moment of nucleated crystals
    :C_high => 0.0,     # Initial over-concentration variable
    :C_low => 0.0,      # Initial under-concentration variable
    :C_high0 => 0.0,   # Initial under-concentration variable, no offset
    :dT_pos => 0.0,     # Initial temperature increase tracking variable
]
param_values = [
    :ρ => 2.66e-12,      # Crystal density [g/μm³]
    :kv => 1.5,          # Volumetric shape factor
    :kb => 285.0,        # Nucleation rate constant [(s·μm³)⁻¹]
    :kg => 1.44e8,       # Growth rate constant [μm/s]
    :Eg_R => 4859.0,     # Growth activation energy / R [K]
    :Eb_R => 7517.0,     # Nucleation activation energy / R [K]
    :b => 1.45,          # Nucleation rate exponent
    :g => 1.5,           # Growth rate exponent
    :h_offset => 0.00,    # Offset for under-concentration variable
]
model = ProcessModel(n_steps);
controls = create_init_controls(n_steps, 323.149, 303.151)
det_prob = DeterministicProblem(model, initial_conditions, param_values, controls);
linear_results, det_ode_prob = solve_controls(det_prob, controls);
plot_states(linear_results, det_prob; plot_title="Linear Temperature Profile Results");
optim_prob = DetermOptimization(det_prob; penalty=10000000.0);
opt_sol_determ = OptimizationSolution(optim_prob; maxiters=30);
opt_results = opt_sol_determ.optimal_solution;
opt_controls = opt_sol_determ.optimal_controls
plot_states(opt_results[1], det_prob; plot_title="Deterministic Optimization Results");

unc_params = [:kb, :kg, :Eg_R, :Eb_R, :b, :g]
det_prob_unc = deepcopy(det_prob);
bo_res = BackoffOptimization(det_prob_unc, unc_params; max_iterations=8, tolerance=0.005);
koop_prob = bo_res.koop_problem;
koop_sol = solve_koop(koop_prob, ireltol=2e-6, quadalg=HCubatureJL(), maxiters=100000)

# Print koop_sol results
println("Koopman Solution Results:")
println("Objective: ", koop_sol.u[1])
println("Cm violation frequency: ", koop_sol.u[2])
println("Cs violation frequency: ", koop_sol.u[3])

plot_states(bo_res.determ_opt_sol.optimal_solution[1], bo_res.determ_problem; plot_title="Optimization with Uncertainty");
objs_mc, l_dev_mc, u_dev_mc = solve_mc_results(koop_prob, 2000000);
avg_objs_mc = mean(objs_mc)
avg_l_dev_mc = mean(l_dev_mc)
avg_u_dev_mc = mean(u_dev_mc)

# Print Monte Carlo average results
println("\nMonte Carlo Results:")
println("Average objective: ", avg_objs_mc)
println("Average Cm violation frequency: ", avg_u_dev_mc)
println("Average Cs violation frequency: ", avg_l_dev_mc)
