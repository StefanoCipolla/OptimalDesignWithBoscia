# Optimal design with FrankWolfe
function solve_opt_frank_wolfe(seed, m, n, time_limit, criterion, corr; write=true, verbose=true)
    o = SCIP.Optimizer()
    if criterion == "A" || criterion == "D"
        A, _, N, ub = build_data(seed, m, n, false, corr)
    elseif criterion == "AF" || criterion == "DF"
        A, C, N, ub = build_data(seed, m, n, true, corr)
    else
        error("Criterion not valid")
    end

    solution = 0.0
    for mode in ["limit", "cont"]
        lmo, x = build_fw_lmo(o, m, N, ub, mode)
        if criterion == "A"
            f, grad! = build_a_criterion(A, false, μ=1e-4)
        elseif criterion == "AF"
            f, grad! = build_a_criterion(A, true, C=C)
        elseif criterion =="D" 
            f, grad! = build_d_criterion(A, false, μ=1e-4)
        elseif criterion == "DF"
            f, grad! = build_d_criterion(A, true, C=C)
        else
            error("Invalid criterion!")
        end

        if criterion == "AF" || criterion == "DF"
            direction = collect(1.0:m)
            x0 = compute_extreme_point(lmo, direction)

            x, _, primal, dual_gap, _, _ = FrankWolfe.blended_pairwise_conditional_gradient(f, grad!, lmo, x0; timeout=10, lazy=true)
            (x, _, primal, dual_gap, traj_data, _), time,_,_,_ = @timed FrankWolfe.blended_pairwise_conditional_gradient(f, grad!, lmo, x0; timeout=time_limit, lazy=true, trajectory = true, max_iteration = 100000)
    
        elseif criterion == "A" || criterion == "D"
            domain_oracle = build_domain_oracle(A, n)
            linesearch = FrankWolfe.Agnostic()
            StepSizeRule = FrankWolfe.MonotonicGenericStepsize(linesearch, domain_oracle)
            _, active_set,_ = start_point(A, m, n, N, ub, "FW_"*mode)
            x, _, primal, dual_gap, _, _ = FrankWolfe.away_frank_wolfe(f, grad!, lmo, active_set; timeout=10, lazy=true, line_search=StepSizeRule)
            _, active_set,_ = start_point(A, m, n, N, ub, "FW_"*mode)
            (x, _, primal, dual_gap, traj_data, _), time,_,_,_ = @timed FrankWolfe.away_frank_wolfe(f, grad!, lmo, active_set; timeout=time_limit, lazy=true, trajectory = true, max_iteration = 100000, line_search=StepSizeRule, verbose = verbose)
    
        end

        iteration = traj_data[end][1]
        if dual_gap ≤ 1e-7
            status = "OPTIMAL"
        elseif time ≥ limit
            status = "Time limit reached"
        elseif iteration ≥ 10000
            status = "Iteration max out"
        end

        # Create integer solution and check feasibility
        #x[x .< 1e-7] .= 0.0
        k = length(findall(x->x!=0, x))
        o = SCIP.Optimizer()
            check_lmo,_ = build_lmo(o,m,N,ub)
        @show Boscia.is_linear_feasible(check_lmo.o, x)
        @show sum(x)
        @show primal
        @show x
        if k >= N
            solution_int = Inf
            feasible = false
        else
            x_int = heuristics(x, N, ub, mode)
            if x_int != Inf
                feasible = Boscia.is_linear_feasible(check_lmo.o, x_int)
                @assert sum(x_int) == N
                solution_int = f(x_int)
            else
                solution_int = Inf
                feasible = false
            end
        end

        type = corr ? "correlated" : "independent"
        
        if write
            # CSV file for the results of all instances.
            df = DataFrame(seed=seed, numberOfExperiments=m, numberOfParameters=n, N=N,time=time, solution_fw=primal, solution_int=solution_int, dual_gap = dual_gap, termination=status, nIteration=iteration, feasibilty=feasible)
            file_name = joinpath(@__DIR__, "../csv/frank_wolfe_" * criterion * "-" * type * "_" *  mode  * "_optimality.csv")
            if !isfile(file_name) 
                CSV.write(file_name, df, append=true, writeheader=true, delim=";")
            else 
                CSV.write(file_name, df, append=true, delim=";")
            end
        end
    end

    return x, solution
end

"""
Build the LMO for FW.
In the limit case:
    min g(A^T diag(x) A)
    s.t. sum x = 1
    x ∈ [0,1]

In the cont case:
    min g(A^T diag(x) A)
    s.t. sum x = s
    x ∈ [0,ub]
"""
function build_fw_lmo(o, m, s, ub, mode)
    MOI.set(o, MOI.Silent(), true)
    MOI.empty!(o)
    x = MOI.add_variables(o, m)
    for i in 1:m
        MOI.add_constraint(o, x[i], MOI.GreaterThan(0.0))
        if mode == "limit"
            MOI.add_constraint(o, x[i], MOI.LessThan(1.0))
        elseif mode == "cont"
            MOI.add_constraint(o, x[i], MOI.LessThan(ub[i]))
        end
    end
    if mode == "limit"
        MOI.add_constraint(
            o,
            MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(ones(m), x), 0.0),
            MOI.EqualTo(1.0)
        )
    elseif mode == "cont"
        MOI.add_constraint(
            o,
            MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(ones(m), x), 0.0),
            MOI.LessThan(s)
        )
        MOI.add_constraint(
            o,
            MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(ones(m), x), 0.0),
            MOI.GreaterThan(1.0)
        )
    end
    lmo = FrankWolfe.MathOptLMO(o)

    return lmo, x
end