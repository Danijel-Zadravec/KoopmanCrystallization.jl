using SciMLExpectations
using Distributions
using SymbolicIndexingInterface
using CSV
using DataFrames
using Statistics
struct KoopmanProblem
    determ_problem::DeterministicProblem
    unc_params
    sm
    gd
    h
    g
    exprob
end

function create_obj()
    function g(sol, p)
        t_end = 1800.0
        y5_final = sol(t_end, idxs=:y5)
        y9_final = sol(t_end, idxs=:y9)
        C_high0_final = sol(t_end, idxs=:C_high0)
        ind_high = C_high0_final > 0.0 ? 1 : 0
        C_low_final = sol(t_end, idxs=:C_low)
        ind_low = C_low_final > 0.0 ? 1 : 0
        observables = [y5_final - y9_final, ind_high, ind_low]

        return observables
    end
    return g
end


function create_obs_all()
    function g(sol, p)
        t_end = 1800.0
        y5_final = sol(t_end, idxs=:y5)
        y9_final = sol(t_end, idxs=:y9)
        C_high0_final = sol(t_end, idxs=:C_high0)
        ind_high = C_high0_final > 0.0 ? 1 : 0
        C_low_final = sol(t_end, idxs=:C_low)
        ind_low = C_low_final > 0.0 ? 1 : 0
        observables = [y5_final - y9_final, ind_high, ind_low]
        return observables
    end
    return g
end

function create_normal_distribution(mean)
    # Create a normal distribution with mean 4859 and standard deviation 4859 * 0.025
    return truncated(Normal(mean, mean * 0.025), mean * 0.95, mean * 1.05)
end

function create_uniform_distribution(mean)
    # Create a uniform distribution around the mean with a small range
    return Uniform(mean * 0.95, mean * 1.05)
end



function create_gd(unc_params, determ_problem)

    prob = determ_problem.problem
    distrs = []
    unc_params =
        for param in unc_params
            push!(distrs, create_normal_distribution(prob.ps[param]))
        end
    return GenericDistribution(distrs...)
end

function create_sm(prob)
    sm = SystemMap(prob, AutoTsit5(Rosenbrock23()), save_everystep=false)
    return sm
end

function create_h(sys, unc_params)
    idxs = [parameter_index(sys, param) for param in unc_params]
    cnt = [0]

    function h(x, u, p)
        for (i, idx) in enumerate(idxs)
            setindex!(p, x[i], idx)
        end
        cnt[1] += 1
        if cnt[1] % 2000 == 0
            println("h called ", cnt[1], " times,")
        end
        return u, p
    end
    return h
end

function KoopmanProblem(problem, unc_params)
    sm = create_sm(problem.problem)
    gd = create_gd(unc_params, problem)
    h = create_h(problem.process_model.model, unc_params)
    g = create_obj()
    exprob = ExpectationProblem(sm, g, h, gd)
    return KoopmanProblem(problem, unc_params, sm, gd, h, g, exprob)
end

function solve_koop(koopman_prob::KoopmanProblem; ireltol=1e-6, quadalg=HCubatureJL(), maxiters=75000)
    # Solve the expectation problem using the Koopman method
    sol = solve(koopman_prob.exprob, Koopman(), ireltol=ireltol, quadalg=quadalg, maxiters=maxiters)
    return sol
end



function modify_ode_prob_params(prob::ODEProblem, sys::ODESystem, unc_params, gd)
    x = rand(gd)
    idxs = [parameter_index(sys, param) for param in unc_params]
    for (i, idx) in enumerate(idxs)
        setindex!(prob.p, x[i], idx)
    end
    return prob, x
end

function solve_mc_ansamble(koopman_prob::KoopmanProblem, repetitions::Int=50000)
    sol = solve(koopman_prob.exprob, MonteCarlo(repetitions))
    return sol
end



function solve_mc_results(koopman_prob::KoopmanProblem, repetitions::Int=1000)
    prob = koopman_prob.determ_problem.problem
    sys = koopman_prob.determ_problem.process_model.model
    unc_params = koopman_prob.unc_params
    gd = koopman_prob.gd
    objs = []
    lbs_dev = []
    ubs_dev = []
    state_names = [:y5, :y9, :C_high0, :C_low]

    # Progress tracking
    progress_interval = max(1, repetitions ÷ 20)  # Show progress every 5%
    for i in 1:repetitions
        # Create a copy of the problem to modify
        prob_copy = remake(prob)

        # Modify parameters using random sampling
        prob_modified, x = modify_ode_prob_params(prob_copy, sys, unc_params, gd)

        # Solve the modified problem - let solver choose timesteps
        sol = solve(prob_modified, Rodas5P(), abstol=1e-8, reltol=1e-6)  # Remove saveat constraint

        obj_val = sol(1800.0, idxs=:y5) - sol(1800.0, idxs=:y9)
        push!(objs, obj_val)
        push!(lbs_dev, sol(1800.0, (idxs = :C_low)) > 0.000001 ? 1 : 0)
        push!(ubs_dev, sol(1800.0, (idxs = :C_high0)) > 0.0000001 ? 1 : 0)
        # Progress reporting
        if i % progress_interval == 0
            println("Completed $i/$repetitions runs")
        end
    end
    return objs, lbs_dev, ubs_dev
end
