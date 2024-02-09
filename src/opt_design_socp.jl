## The SOCP formulation form "COMPUTING EXACT D-OPTIMAL DESIGNS BY MIXED INTEGER SECOND-ORDER CONE PROGRAMMING" 
## by Guillaume Sagnol and Radoslav Harman
"""
Build SOCP of A-Opt.
"""
function build_A_socp_model(seed, m, n, criterion, time_limit, corr, verbose; lb=nothing, ub = nothing)
    if criterion == "AF" 
        A, C, N, ubounds = build_data(seed, m, n, true, corr)
    else
        A, _, N, ubounds = build_data(seed, m, n, false, corr)
        @assert N ≥ n
    end
    @show m, n, N, sum(ubounds) 
    @assert (m > n) && (sum(ubounds) >= N)

    if lb === nothing
        lb = fill(0.0, m)
    end
    if ub === nothing
        ub = ubounds
    end

    # setup solvers
    # MIP solver (try SCIP as well?)
    oa_solver = optimizer_with_attributes(HiGHS.Optimizer,
        MOI.Silent() => !verbose,
        "mip_feasibility_tolerance" => 1e-10,
        "mip_rel_gap" => 1e-8,
    )
    # SDP solver
    conic_solver = optimizer_with_attributes(Hypatia.Optimizer, 
        MOI.Silent() => !verbose,
    )
    opt = optimizer_with_attributes(Pajarito.Optimizer,
        "time_limit" => time_limit, 
        "iteration_limit" => 1000000,
        "oa_solver" => oa_solver, 
        "conic_solver" => conic_solver,
        MOI.Silent() => !verbose,
    )

    model = Model(opt)
    # add design variables
    JuMP.@variable(model, w[1:m])
    JuMP.@variable(model, x[1:m], Int)
    # Add the probability simplex constraint
    #JuMP.@constraint(model, sum(x) == N)
    JuMP.@constraint(model, sum(w) == 1)
    # Bound constraints
    JuMP.@constraint(model, x ≤ ub)
    JuMP.@constraint(model, x ≥ lb)
    JuMP.@constraint(model, w ≥ fill(0.0, m))

    # auxiliary variables 
    # μ 
    JuMP.@variable(model, μ[1:m])
    # Y_i collected in matrix
    JuMP.@variable(model, Y_mat[1:m, 1:n])

    # additional constraints
    # linear constraints
    JuMP.@constraint(model, x .== N*w)
    JuMP.@constraint(model, sum(transpose(A[i, :]) * Y_mat[i,:] for i in 1:m) .== sum(μ)*I)
    JuMP.@constraint(model, μ ≥ fill(0.0, m))
    # conic constraint
    # second order cone
    for i in 1:m
        JuMP.@constraint(model, vcat(μ[i] + w[i], vcat(μ[i] - w[i], Y_mat[i,:])) in Hypatia.EpiNormEuclCone{Float64}(n+2))
        #JuMP.@constraint(model, vcat(1/2*μ[i], w[i], Y_mat[i,:]) in Hypatia.EpiPerSquareCone{Float64}(n+2))
    end

    # objective function
    JuMP.@objective(model, Max, sum(μ))

    return model, x
end

"""
Build the SOCP of the D-criterion
"""
function build_D_socp_model(seed, m, n, criterion, time_limit, corr, verbose; lb=nothing)
    if criterion == "DF" 
        A, C, N, ub = build_data(seed, m, n, true, corr)
    else
        A, _, N, ub = build_data(seed, m, n, false, corr)
        @assert N ≥ n
    end
    @show m, n, N, sum(ub) 
    @assert (m > n) && (sum(ub) >= N)

    if lb === nothing
        lb = fill(0.0, m)
    end

    # setup solvers
    # MIP solver (try SCIP as well?)
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
    # add design variables
    JuMP.@variable(model, w[1:m])
    JuMP.@variable(model, x[1:m], Int)
    # Add the probability simplex constraint
    #JuMP.@constraint(model, sum(x) == N)
    JuMP.@constraint(model, sum(w) == 1)
    # Bound constraints
    JuMP.@constraint(model, x ≤ ub)
    JuMP.@constraint(model, x ≥ lb)
    JuMP.@constraint(model, w ≥ fill(0.0, m))

    # auxiliary variables 
    # J 
    JuMP.@variable(model, J_mat[1:n, 1:n])
    for i in 1:n 
        for j in i:n
            JuMP.fix(J_mat[i,j], 0; force=true)
        end
    end
    # Z_i collect in a matrix
    JuMP.@variable(model, Z_mat[1:m, 1:n])
    # t_ij collected in a matrix
    JuMP.@variable(model, T_mat[1:m, 1:n])
    # Epigraph variable
    JuMP.@variable(model, t)

    # additional constraints
    # linear constraints
    JuMP.@constraint(model, x .== N*w)
    JuMP.@constraint(model, sum(transpose(A[i, :]) * Z_mat[i,:] for i in 1:m) .== J_mat)
    for j in 1:n 
        JuMP.@constraint(model, sum(T_mat[i,j] for i in 1:m) ≤ J_mat[j,j])
    end
    # conic constraint
    # second order cone
    JuMP.@constraint(model, vcat(t, LinearAlgebra.diag(J_mat)) in Hypatia.HypoGeoMeanCone{Float64}(n+1))
    for i in 1:m
        for j in 1:n 
          #  JuMP.@constraint(model, vcat(T_mat[i,j] + x[i], Z_mat[i,j]) in Hypatia.EpiNormEuclCone{Float64}(2))
            JuMP.@constraint(model, vcat(1/2*T_mat[i,j], w[i], Z_mat[i,j]) in Hypatia.EpiPerSquareCone{Float64}(3))
        end
    end

    # objective function
    JuMP.@objective(model, Max, t)

    return model, x
end

function solve_opt_socp(seed, m, n, time_limit, criterion, corr; write=true, verbose=true)
    if criterion == "DF" || criterion == "D"
        model, x = build_D_socp_model(seed, m, n, criterion, 10, corr, false)
        optimize!(model)
        model, x = build_D_socp_model(seed, m, n, criterion, time_limit, corr, verbose)
    elseif criterion == "AF" || criterion == "A"
        model, x= build_A_socp_model(seed, m, n, criterion, 10, corr, false)
        optimize!(model)
        model, x = build_A_socp_model(seed, m, n, criterion, time_limit, corr, verbose)
    end

    # solve 
    optimize!(model)

    # query solution
    status = termination_status(model)
    solution = objective_value(model)
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
    scaled_solution = if feasible
        f(y)
    else
        Inf
    end

    # o = JuMP.moi_backend(model)
    type = corr ? "correlated" : "independent"
    @show y
    @show solution
    @show scaled_solution
    @show feasible
    @show status

    if write 
        df = DataFrame(seed=seed, numberOfExperiments=m, numberOfParameters=n, time=t, N=N, solution=solution, scaled_solution=scaled_solution, termination=status, numberIterations=numberIter, numberCuts=numberCuts, feasible = feasible)
        file_name = "/scratch/opt/dhendryc/MasterThesis/optDesign/csv/SOCP/socp_" * criterion * "_" * string(m) * "_" * type * "_optimality.csv"
        if !isfile(file_name)
            CSV.write(file_name, df, append=true, writeheader=true)
        else 
            CSV.write(file_name, df, append=true)
        end
    end
    return y
end