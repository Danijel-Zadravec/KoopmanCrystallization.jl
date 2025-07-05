using SymbolicIndexingInterface
using SymbolicIndexingInterface: parameter_values, state_values
using SciMLStructures
using SciMLStructures: Tunable, canonicalize, replace, replace!
using PreallocationTools
using ForwardDiff
using Optimization
using Plots

struct DetermOptimization
    deterministic_problem::DeterministicProblem
    loss_params
    objective_function
    initial_control
    optimization_function::OptimizationFunction
    optimization_problem::OptimizationProblem
end


function create_loss_params(determ_problem::DeterministicProblem; penalty=100000.0)
    # Extract the ODEProblem from the DeterministicProblem
    ode_problem = determ_problem.problem
    proc_model = determ_problem.process_model.model
    # Get the parameter values from the ODEProblem
    param_values = parameters(proc_model)
    # Create a DiffCache for caching parameter values
    diff_cache = DiffCache(copy(canonicalize(Tunable(), parameter_values(ode_problem))[1]))

    u_params = determ_problem.process_model.symbols.controls
    # Create a setter function for updating parameters
    setter_fn = setp(ode_problem, u_params)
    return [penalty, ode_problem, diff_cache, setter_fn, proc_model]
end


function create_objective_function(; abstol=1e-10, reltol=1e-8)
    function loss(u, p)

        # Penalty weights for bound violations
        penalty_weight = p[1]

        odeprob = p[2] # ODEProblem stored as parameters to avoid using global variables
        ps = parameter_values(odeprob) # obtain the parameter object from the problem
        diffcache = p[3]
        buffer = get_tmp(diffcache, u)
        copyto!(buffer, canonicalize(Tunable(), ps)[1])
        ps = replace(Tunable(), ps, buffer)
        setter = p[4]
        setter(ps, u)
        newprob = remake(odeprob; p=ps)
        sol_obj = solve(newprob, Rodas5P(), abstol=abstol, reltol=reltol)

        # Extract final values efficiently
        final_state = sol_obj.u[end]
        proc_model = p[5]
        all_vars = unknowns(proc_model)
        all_vars = unknowns(proc_model)

        # Get indices for main objective variables
        y5_idx = findfirst(x -> string(x) == "y5(t)", string.(all_vars))
        y9_idx = findfirst(x -> string(x) == "y9(t)", string.(all_vars))
        Ch_ind = findfirst(x -> string(x) == "C_high(t)", string.(all_vars))
        Cl_ind = findfirst(x -> string(x) == "C_low(t)", string.(all_vars))
        dT_idx = findfirst(x -> string(x) == "dT_pos(t)", string.(all_vars))



        y5_final = final_state[y5_idx]
        y9_final = final_state[y9_idx]
        # Calculate objective: (S1(T) + 0.5*S2(T))/P(T)
        main_objective = -(y5_final - y9_final) / 1000000.0
        total_penalty = 0.0
        total_penalty += final_state[Ch_ind] * 10000
        total_penalty += final_state[Cl_ind] * 10000
        total_penalty += final_state[dT_idx]

        # Total objective with penalties
        total_objective = main_objective + penalty_weight * total_penalty^2
        if isa(total_objective, AbstractFloat)
            println("Objective $(total_objective)")
        end
        return total_objective
    end
    return loss
end

function DetermOptimization(determ_problem::DeterministicProblem; penalty=100000.0)
    loss_params = create_loss_params(determ_problem; penalty=penalty)
    objective_function = create_objective_function()
    optfn = OptimizationFunction(objective_function, Optimization.AutoForwardDiff())
    n_steps = determ_problem.process_model.n_steps
    lower = fill(273.15 + 30, n_steps + 1)
    upper = fill(273.15 + 50, n_steps + 1)
    initial_controls = [pair.second for pair in determ_problem.initial_controls]
    optprob = OptimizationProblem(
        optfn, initial_controls, loss_params,
        lb=lower, ub=upper)
    return DetermOptimization(determ_problem, loss_params, objective_function, initial_controls, optfn, optprob)
end


struct OptimizationSolution
    optimization_problem::DetermOptimization
    opt_result::Any  # The optimization result from solve()
    optimal_controls::Vector{Float64}
    optimal_solution::Any  # The ODE solution with optimal controls
    solve_time::Float64
    success::Bool
    objective_value::Float64
end

#function opt_solve(determ_opt::DetermOptimization; iterations=100, maxiters=10, trace=true, ext_trace=true, linesearch=HagerZhang())
function opt_solve(determ_opt::DetermOptimization; maxiters=20, reltol=1e-6)

    # Unpack the optimization problem
    optprob = determ_opt.optimization_problem
    sol = solve(optprob, Optimization.LBFGS(), maxiters=maxiters, reltol=reltol)

    return sol
end


#function OptimizationSolution(determ_opt::DetermOptimization; iterations=100, maxiters=10, trace=true, ext_trace=true, linesearch=HagerZhang())
function OptimizationSolution(determ_opt::DetermOptimization; maxiters=20)

    """
    Solve the optimization problem and create OptimizationSolution
    
    Args:
        determ_opt: DetermOptimization containing problem setup
    
    Returns:
        OptimizationSolution with results
    """

    # Record start time
    start_time = time()

    # Solve the optimization problem
    #opt_result = opt_solve(determ_opt::DetermOptimization; iterations=iterations, maxiters=maxiters, trace=trace, ext_trace=ext_trace, linesearch=linesearch)
    opt_result = opt_solve(determ_opt::DetermOptimization; maxiters=maxiters)

    # Calculate solve time
    solve_time = time() - start_time

    # Extract optimal controls
    optimal_controls = opt_result.u
    # Check if optimization was successful
    success = opt_result.retcode == ReturnCode.Success || opt_result.retcode == :Success

    # Get objective value
    objective_value = opt_result.objective

    # Solve the ODE
    optimal_solution = solve_controls(determ_opt.deterministic_problem, optimal_controls; abstol=1e-10, reltol=1e-8)

    return OptimizationSolution(
        determ_opt,
        opt_result,
        optimal_controls,
        optimal_solution,
        solve_time,
        success,
        objective_value
    )
end

