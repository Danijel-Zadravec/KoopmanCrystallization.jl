module KoopmanCrystallization

# Import main modules
include("process_model.jl")
include("deterministic_problem.jl")
include("deterministic_optimization.jl")
include("koopman.jl")
include("backoff.jl")

# Export main types
export ProcessModel, ModelSymbols
export DeterministicProblem, DetermOptimization, OptimizationSolution
export KoopmanProblem, BackoffOptimization

# Export main functions
export create_process_model, create_init_controls
export solve_controls, modify_controls, modify_backoff, set_backoff
export solve_koop, solve_mc_results
export plot_states

end
