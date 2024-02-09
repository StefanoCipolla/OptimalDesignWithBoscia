## Test the SOCP formulation
using optDesign
using Boscia
using Pajarito
using SCIP
using Hypatia
using FrankWolfe
using Random
using Test
using LinearAlgebra
using MathOptInterface
const MOI = MathOptInterface
using JuMP


seed = 2
Random.seed!(seed)

dimensions = [20,30,40] # 30,40
facs = [10,4]#[4,6]
time_limit = 300 #1200
verbose = false

@testset "A Optimal" begin
    for m in dimensions 
        for k in facs 
            n = Int(floor(m/k))
            # independent data
            x_b_ind, _ = optDesign.solve_opt(seed, m ,n, time_limit, "A", false; write = false, verbose=verbose)
            x_socp = optDesign.solve_opt_socp(seed, m, n, time_limit, "A", false; write=false, verbose=verbose)

            # Testing
            A, C, N, ub = optDesign.build_data(seed, m , n, false, false)
            check_lmo = optDesign.build_blmo(m, N, ub)
            f, _ = optDesign.build_a_criterion(A, false)

            domain_oracle = optDesign.build_domain_oracle(A,n)
            function g(x)
               # if domain_oracle(x)
                    X = LinearAlgebra.Symmetric(transpose(A)*LinearAlgebra.diagm(x)*A)
                    return (LinearAlgebra.tr(X^(-1)))^(-1)
                #end
                #return 0
            end

            @show g(x_b_ind)
            @show g(x_socp)

            @test g(x_b_ind) >= g(x_socp)

            @test Boscia.is_linear_feasible(check_lmo, x_socp)
            if optDesign.isfeasible(seed, m, n, "A", x_socp, false)
                @show f(x_b_ind)
                @show f(x_socp)
                @test f(x_b_ind) ≤ f(x_socp)
            end

            model, y = optDesign.build_A_socp_model(seed, m, n, "A", time_limit, false, false, lb=x_b_ind, ub=x_b_ind)
            optimize!(model)
            status = termination_status(model)
            solution = objective_value(model)
            y = round.(value.(y))
            @show solution
            @show y
        end
    end
end 

@testset "D Optimal" begin
    for m in dimensions 
        for k in facs 
            n = Int(floor(m/k))
            # independent data
            x_b_ind, _ = optDesign.solve_opt(seed, m ,n, time_limit, "D", false; write = false, verbose=verbose)
            x_socp = optDesign.solve_opt_socp(seed, m, n, time_limit, "D", false; write=false, verbose=verbose)

            # Testing
            A, C, N, ub = optDesign.build_data(seed, m , n, false, false)
            check_lmo = optDesign.build_blmo(m, N, ub)
            f, _ = optDesign.build_d_criterion(A, false)

            domain_oracle = optDesign.build_domain_oracle(A,n)
            function g(x)
                if domain_oracle(x)
                    X = LinearAlgebra.Symmetric(transpose(A)*LinearAlgebra.diagm(x)*A)
                    return (det(X))^(1/n)
                end
                return 0
            end

            @show g(x_b_ind)
            @show g(x_socp)

            @test Boscia.is_linear_feasible(check_lmo, x_socp)
            if optDesign.isfeasible(seed, m, n, "D", x_socp, false)
                @show f(x_b_ind)
                @show f(x_socp)
                @test f(x_b_ind) ≤ f(x_socp)
            end
            model, y = optDesign.build_A_socp_model(seed, m, n, "D", time_limit, false, false, lb=x_b_ind)
            optimize!(model)
            status = termination_status(model)
            solution = objective_value(model)
            y = value.(y)
            @show solution
            @show y
        end
    end 
end 

@testset "A Fusion" begin
 #=   for m in dimensions 
        for k in facs 
            n = Int(floor(m/k))
        
            x_b_corr, _ = optDesign.solve_opt(seed, m, n, time_limit, "AF", false; write = false, verbose=verbose)
            x_socp = optDesign.solve_opt_socp(seed, m, n, time_limit, "AF", false; write=false, verbose=verbose)

            # Testing
            A, C, N, ub = optDesign.build_data(seed, m , n, true, false)
            check_lmo = optDesign.build_blmo(m, N, ub)
            f, _ = optDesign.build_a_criterion(A, true, C=C)


            @test Boscia.is_linear_feasible(check_lmo, x_socp)
           if optDesign.isfeasible(seed, m, n, "AF", x_socp, false)
                @show f(x_b_corr)
                @show f(x_socp)
                @test f(x_b_ind) ≤ f(x_socp)
            end
        end
    end =#
end

@testset "D Fusion" begin
  #=  for m in dimensions 
        for k in facs 
            n = Int(floor(m/k))
            
            x_b_corr, _ = optDesign.solve_opt(seed, m, n, time_limit, "DF", false; write = false, verbose=verbose)
            x_socp = optDesign.solve_opt_socp(seed, m, n, time_limit, "DF", false; write=false, verbose=verbose)

            # Testing
            A, C, N, ub = optDesign.build_data(seed, m , n, true, false)
            check_lmo = optDesign.build_blmo(m, N, ub)
            f, _ = optDesign.build_d_criterion(A, true, C=C)


            @test Boscia.is_linear_feasible(check_lmo, x_socp)
            if optDesign.isfeasible(seed, m, n, "DF", x_socp, false)
                @show f(x_b_corr)
                @show f(x_socp)
                @test f(x_b_ind) ≤ f(x_socp)
            end
        end
    end =#
end 
