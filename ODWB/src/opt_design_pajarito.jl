#### Solving opt design problem with Pajarito

# Pajarito model for the D-optimal problems
function build_D_pajarito_model(seed, m, n, criterion, time_limit, corr, verbose=true)
    if criterion == "DF" 
        A, C, N, ub = build_data(seed, m, n, true, corr)
    else
        A, _, N, ub = build_data(seed, m, n, false, corr)
        @assert N ≥ n
    end
    @show m, n, N, sum(ub) 
    @assert (m > n) && (sum(ub) >= N)

    # setup solvers
    # MIP solver
    oa_solver = optimizer_with_attributes(HiGHS.Optimizer,
        MOI.Silent() => !verbose,
    )
    # SDP solver
    conic_solver = optimizer_with_attributes(Hypatia.Optimizer, 
        MOI.Silent() => !verbose,
    )
    opt = optimizer_with_attributes(Pajarito.Optimizer,
        "time_limit" => time_limit, 
        "iteration_limit" => 100000,
        "oa_solver" => oa_solver, 
        "conic_solver" => conic_solver,
        MOI.Silent() => !verbose,
    )

    model = Model(opt)
    # add variables
    JuMP.@variable(model, x[1:m], Int)
    # we want to do s experiments
    JuMP.@constraint(model, sum(x) == N)

    # Constraints on the total times each experiment can be run
    ub_u = copy(ub)
    unique!(ub_u)
    for u in ub_u
        ind = findall(x->x==u, ub)
        mid = u / 2
        JuMP.@constraint(model, vcat(mid, x[ind] .- mid) in MOI.NormInfinityCone(length(ind) + 1))
    end

    JuMP.@variable(model, t)

    # information matrix lower triangle
    if criterion == "D"
        a1 = [
            JuMP.@expression(model, sum(A[k, i] * x[k] * A[k, j] for k in 1:m)) for
            i in 1:n for j in 1:i
        ]
        JuMP.@constraint(model, vcat(t, 1, a1) in MOI.LogDetConeTriangle(n))
    elseif criterion == "DF"
        a1 = [
            JuMP.@expression(model, C[i,j] + sum(A[k, i] * x[k] * A[k, j] for k in 1:m)) for
            i in 1:n for j in 1:i
        ]
        JuMP.@constraint(model, vcat(t, 1, a1) in MOI.LogDetConeTriangle(n))
        
    end
    JuMP.@objective(model, Max, t)

    return model, x
end

# Pajarito model for the A-optimal problems
# As suggested here: https://github.com/jump-dev/Pajarito.jl/issues/444
function build_A_pajarito_model(seed, m , n, criterion, time_limit, corr, verbose = true)
    if criterion == "AF"
        A, C, N, ub = build_data(seed, m, n, true, corr)
    else
        A, _, N, ub = build_data(seed, m, n, false, corr)
        @assert N ≥ n
    end
    @assert (m > n) && (sum(ub) >= N)

    # setup solvers
    # MIP solver
    oa_solver = optimizer_with_attributes(HiGHS.Optimizer,
        MOI.Silent() => !verbose,
    )
    # SDP solver
    conic_solver = optimizer_with_attributes(Hypatia.Optimizer, 
        MOI.Silent() => !verbose,
    )
    opt = optimizer_with_attributes(Pajarito.Optimizer,
        "time_limit" => time_limit, 
        "iteration_limit" => 100000,
        "oa_solver" => oa_solver, 
        "conic_solver" => conic_solver,
        MOI.Silent() => !verbose,
    )

    model = Model(opt)
    # add variables
    JuMP.@variable(model, x[1:m])
    JuMP.set_integer.(x)
    JuMP.@variable(model, t)
    # we want to do s experiments
    JuMP.@constraint(model, sum(x) == N)
    @objective(model, Min, 4 * t)

    # Constraints on the total times each experiment can be run
    ub_u = copy(ub)
    unique!(ub_u)
    for u in ub_u
        ind = findall(x->x==u, ub)
        mid = u / 2
        JuMP.@constraint(model, vcat(mid, x[ind] .- mid) in MOI.NormInfinityCone(length(ind) + 1))
    end

    if criterion == "A"
        # vectorized information matrix
        X_vec = [
            JuMP.@expression(
                model,
                (i == j ? 1.0 : sqrt(2)) * sum(A[k, i] * x[k] * A[k, j] for k in 1:m)
            ) for i in 1:n for j in 1:i
        ]
        add_homog_spectral(MatNegSqrtConj() , n, vcat(1.0 * t, X_vec), model)
    elseif criterion == "AF"
        X_vec = [
            JuMP.@expression(
                model,
                (i == j ? 1.0 : sqrt(2)) * (C[i,j] + sum(A[k, i] * x[k] * A[k, j] for k in 1:m))
            ) for i in 1:n for j in 1:i
        ]
        add_homog_spectral(MatNegSqrtConj() , n, vcat(1.0 * t, X_vec), model)
    end

    return model, x
end


function solve_opt_pajarito(seed, m, n, time_limit, criterion, corr; write=true, verbose=true)
    if criterion == "DF" || criterion == "D"
        model, x = build_D_pajarito_model(seed, m, n, criterion, 10, corr, false)
        optimize!(model)
        model, x = build_D_pajarito_model(seed, m, n, criterion, time_limit, corr, verbose)
    elseif criterion == "AF" || criterion == "A"
        model, x= build_A_pajarito_model(seed, m, n, criterion, 10, corr, false)
        optimize!(model)
        model, x = build_A_pajarito_model(seed, m, n, criterion, time_limit, corr, verbose)
    end

    # solve 
    optimize!(model)

    # query solution
    status = termination_status(model)
    solution = objective_value(model)
    solution = criterion == "D" || criterion == "DF" ? solution * (-1) : solution
    y = value.(x)
    t = solve_time(model)
    paja_opt = JuMP.unsafe_backend(model)
    numberIter = paja_opt.num_iters
    numberCuts = paja_opt.num_cuts

    # Check feasibility
    if criterion == "A" || criterion == "D"
        A, _, N, ub = build_data(seed, m, n, false, corr)
    elseif criterion == "AF"|| criterion == "DF"
        A, C, N, ub = build_data(seed, m, n, true, corr)
    end
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
    feasible = isfeasible(seed, m, n,criterion, y, corr)
    @show feasible

    # o = JuMP.moi_backend(model)
    type = corr ? "correlated" : "independent"
    @show y
    @show solution
    @show f(y)
    @show feasible

    if write 
        df = DataFrame(seed=seed, numberOfExperiments=m, numberOfParameters=n, time=t, N=N, solution=solution, termination=status, numberIterations=numberIter, numberCuts=numberCuts)
        file_name = joinpath(@__DIR__, "../csv/Pajarito/pajarito_" * criterion * "_" * string(m) * "_" * type * "_optimality.csv")
        if !isfile(file_name)
            CSV.write(file_name, df, append=true, writeheader=true)
        else 
            CSV.write(file_name, df, append=true)
        end
    end
    @assert feasible
    return y
end
