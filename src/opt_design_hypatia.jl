# Solving the optimal design problem using Hypatia
# Adapting example from the Hypatia package: https://github.com/chriscoey/Hypatia.jl/blob/master/examples/doptimaldesign/JuMP.jl

function build_hypatia_model(seed, m, n, criterion, time_limit, mode, corr, verbose)
    if criterion == "A" || criterion == "D"
        A, _, N, ub = build_data(seed, m, n, false, corr)
    elseif criterion == "AF"|| criterion == "DF"
        A, C, N, ub = build_data(seed, m, n, true, corr)
    end

    @assert (m > n) && (maximum(ub) <= N)
    if criterion == "A" || criterion == "D"
        @assert N ≥ n 
    end

    opt = Hypatia.Optimizer(verbose = verbose)
    model = Model(() -> opt)
    set_time_limit_sec(model, time_limit)
    # add variables
    JuMP.@variable(model, x[1:m])
    if mode == "limit"
        # we want to do s experiments in cont case 
        JuMP.@constraint(model, sum(x) == 1)
        mid = 1 / 2
        JuMP.@constraint(model, vcat(mid, x[1:m] .- mid) in MOI.NormInfinityCone(m + 1))

    elseif mode == "cont"
        # we want to do s experiments in cont case 
        JuMP.@constraint(model, sum(x) == N)

        # Constraints on the total times each experiment can be run
        ub_u = copy(ub)
        unique!(ub_u)
        for u in ub_u
            ind = findall(x->x==u, ub)
            mid = u / 2
            JuMP.@constraint(model, vcat(mid, x[ind] .- mid) in MOI.NormInfinityCone(length(ind) + 1))
        end
    end

    JuMP.@variable(model, t)

    # information matrix lower triangle
    if criterion == "D"
        a1 = [
            JuMP.@expression(model, sum(A[k, i] * x[k] * A[k, j] for k in 1:m)) for
            i in 1:n for j in 1:i
        ]
        JuMP.@constraint(model, vcat(t, 1, a1) in MOI.LogDetConeTriangle(n))
        JuMP.@objective(model, Max, t)
    elseif criterion == "DF"
        a1 = [
            JuMP.@expression(model, C[i,j] + sum(A[k, i] * x[k] * A[k, j] for k in 1:m)) for
            i in 1:n for j in 1:i
        ]
        JuMP.@constraint(model, vcat(t, 1, a1) in MOI.LogDetConeTriangle(n))
        JuMP.@objective(model, Max, t)
    elseif criterion == "A"
        # https://discourse.julialang.org/t/how-to-optimize-trace-of-matrix-inverse-with-jump-or-convex/94167/4
        #@constraint(model, vcat(1.0, t1, [S[i, j]  * (i == j ? 1 : sqrt(2)) for i in 1:p for j in 1:i]...) in cone)
        cone = EpiPerSepSpectralCone{Float64}(Hypatia.Cones.NegSqrtSSF(), Hypatia.Cones.MatrixCSqr{Float64, Float64}, n, true)
        @constraint(model, vcat(1.0, t, [ sum(A[k, i] * x[k] * A[k,j] for k in 1:m) * (i == j ? 1 : sqrt(2)) for i in 1:n for j in 1:i]...) in cone)
        @objective(model, Min, 4 * t)
    elseif criterion == "AF"
        cone = EpiPerSepSpectralCone{Float64}(Hypatia.Cones.NegSqrtSSF(), Hypatia.Cones.MatrixCSqr{Float64, Float64}, n, true)
        @constraint(model, vcat(1.0, t, [ (C[i,j] + sum(A[k, i] * x[k] * A[k,j] for k in 1:m)) * (i == j ? 1 : sqrt(2)) for i in 1:n for j in 1:i]...) in cone)
        @objective(model, Min, 4 * t)
    end

    return model, x, N
end

function solve_opt_hypatia(seed, m, n, time_limit, criterion, corr; write=true, verbose = true)
    for mode in ["limit", "cont"]
        # precompilation
        model, x, N = build_hypatia_model(seed, m, n, criterion, 10, mode, corr, false)
        optimize!(model)

        # set up model 
        model, x, N = build_hypatia_model(seed, m, n, criterion, time_limit, mode, corr, verbose)
        # solve and query solution
        optimize!(model)
        status= termination_status(model)
        solution = objective_value(model)
        solution = criterion == "D" || criterion == "DF" ? solution * (-1) : solution
        #primal_status = primal_status(model) # MOI.FEASIBLE_POINT
        y = value.(x)
        t = solve_time(model)

        # Check feasibility
        if criterion == "A" || criterion == "D"
            n, A, _, s, ub = build_data(seed, m, n, false, corr)
        elseif criterion == "AF"|| criterion == "DF"
            n, A, C, s, ub = build_data(seed, m, n, true, corr)
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
        
        # preprossesing to avoid some small values
        o = SCIP.Optimizer()
            check_lmo,_ = build_lmo(o,m,N,ub)
            @show Boscia.is_linear_feasible(check_lmo.o,y)
            @show sum(y)
            @show solution
            @show y
        y[y .< 1e-7] .= 0.0
        k = length(findall(x->x!=0, y))
        if k >= N
            solution_int = Inf
            feasible = false
        else
            @show Boscia.is_linear_feasible(check_lmo.o, y)
            x_int = heuristics(y, N, ub, mode)
            if x_int != Inf
                feasible = Boscia.is_linear_feasible(check_lmo.o, x_int)
                @assert sum(x_int) == N
                solution_int = f(x_int)

                @show solution_int 
            else
                solution_int = Inf
                feasible = false
            end
        end
        
        if write
            type = corr ? "correlated" : "independent"

            df = DataFrame(seed=seed, numberOfExperiments=m, numberOfParameters=n, N=N,time=t, solution=solution, solution_int=solution_int, termination=status, feasible = feasible)
            file_name = joinpath(@__DIR__,"../csv/hypatia_" * criterion * "-" * type * "_" *  mode * "_optimality" * ".csv")
            if !isfile(file_name)
                CSV.write(file_name, df, append=true, writeheader=true)
            else 
                CSV.write(file_name, df, append=true)
            end
        end
    end
    return y
end




