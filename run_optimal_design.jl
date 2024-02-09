# Script for running the experiments
using optDesign
using Printf
using Test
using DataFrames
using CSV

"""
DON'T FORGET TO ADD e AGAIN ONCE YOU ARE DONE DEBUGGING!!
"""
mode = ENV["MODE"]
criterion = ENV["CRITERION"]
type = ENV["TYPE"]
m = parse(Int, ENV["DIMENSION"])
corr = if type == "IND"
    false 
elseif type == "CORR"
    true
else 
    error("Type not found")
end
#run(`hostname`)
#=mode = "Pajarito"
criterion = "A"
corr = true
m = 30 #number of experiments
=#

#num_experiments = [20,30,50,60,80,100,120,150,180,200]
#num_experiments = [250,300,350,400] # Only for Boscia!
ratio_para = [4,10]#[4,10]#[4, 6, 8, 10]
#ratio_para = [4]
time_limit = 3600 # one hour time limit

@show criterion, mode, corr

if !(criterion in ["A", "D", "DF", "AF"])
    error("Invalid criterion!")
end
#for m in num_experiments
    for k in ratio_para
        n = Int(floor(m/k))
        for seed in 1:5
           # if (n==25 && seed in [1]) #|| (n==22 && seed in [1,3,4,5]) || (n==18 && seed in [1,2,4,5]) #|| (n==20 && seed in [1,3])
           #     continue
           # end
            @show m, n, seed
            try
                if mode == "Boscia"
                    optDesign.solve_opt(seed, m, n, time_limit, criterion,corr)
                elseif mode == "SCIP"
                   # if criterion in ["A", "D"]
                    #    error("SCIP OA does not work with the optimal problems!")
                    #end
                    optDesign.solve_opt_scip(seed, m, n, time_limit, criterion, corr)
                elseif mode == "Pajarito"
                    optDesign.solve_opt_pajarito(seed, m, n, time_limit, criterion, corr)
                elseif mode == "Custom"
                   # if criterion in ["AF", "DF"]
                    #    error("The custom BB does not work with the fusion problems!")
                    #end
                    optDesign.solve_opt_custom(seed, m, n, time_limit, criterion, corr)
                elseif mode == "FrankWolfe"
                    optDesign.solve_opt_frank_wolfe(seed, m, n, time_limit, criterion, corr)
                elseif mode == "Hypatia"
                    optDesign.solve_opt_hypatia(seed, m , n, time_limit, criterion, corr)
                else 
                    error("Invalid mode!")
                end
            catch e
                println(e)
                error_file = criterion * "_opt_" * mode * "_" * type * ".txt" 
                open(error_file,"a") do io
                    println(io, seed, " ", m, " ", mode, " : ", e)
                end
            end
        end
    end
#end
#=
if criterion == "A"
    for m in dimensions
        for k in ratio_para
            n = Int(floor(m/k))
            for seed in 1:5
                @show m, n, seed
                try
                    if mode == "Boscia"
                        optDesign.solve_opt(seed, m, criterion,corr=corr)
                    elseif mode == "SCIP"
                        #error("SCIP OA does not work with this setting!")
                    elseif mode == "FrankWolfe"
                        optDesign.solve_opt_frank_wolfe(seed, m, time_limit, criterion, corr)
                    elseif mode == "Hypatia"
                        optDesign.solve_opt_hypatia(seed, m , n, time_limit, criterion, corr)
                    elseif mode == "Pajarito"
                        optDesign.solve_opt_pajarito(seed, m, n, time_limit, criterion, corr)
                    elseif mode == "Custom"
                        optDesign.solve_opt_custom(seed, m, n, time_limit, criterion, corr)
                    else 
                        error("Invalid mode!")
                    end
                catch 
                    println(e)
                    error_file = criterion * "_opt_" * mode * "_" * type * ".txt" #[20,30,50,60,80,100,
                    open(error_file,"a") do io
                        println(io, seed, " ", m, " ", mode, " : ", e)
                    end
                end
            end
        end
    end
elseif criterion == "D"
    for m in [20,30]#[20,30] #[20,30,50,60,80,100,150,200,250,300]
        for seed in 1:5
            @show m, seed
            try
                if mode == "Boscia"
                    optDesign.solve_opt(seed, m, criterion, corr=corr)
                elseif mode == "SCIP"
                    #error("SCIP OA does not work with this setting!")
                elseif mode == "FrankWolfe"
                    optDesign.solve_opt_frank_wolfe(seed, m, criterion, corr)
                elseif mode == "Hypatia"
                    optDesign.solve_opt_hypatia(seed, m ,criterion, corr)
                elseif mode == "Pajarito"
                    optDesign.solve_opt_pajarito(seed, m, criterion)
                else
                    error("Invalid mode!")
                end
            catch 
                println(e)
                error_file = criterion * "_opt_" * mode * "_" * type * ".txt"
                open(error_file,"a") do io
                    println(io, seed, " ", m, " ", mode, " : ", e)
                end
            end
        end
    end
elseif criterion == "DF"
    for m in [20,30]#[20,30,50,60,80,100,150,200,250,300]
        for seed in 1:5 #1:5
            @show m, seed
            try
                if mode == "Boscia"
                    optDesign.solve_opt(seed, m, criterion,corr=corr)
                elseif mode == "SCIP"
                    optDesign.solve_opt_scip(seed, m, criterion, corr)
                elseif mode == "FrankWolfe"
                    optDesign.solve_opt_frank_wolfe(seed, m, criterion, corr)
                elseif mode == "Hypatia"
                    optDesign.solve_opt_hypatia(seed, m ,criterion, corr)
                elseif mode == "Pajarito"
                    optDesign.solve_opt_pajarito(seed, m, criterion)
                else
                    error("Invalid mode!")
                end
            catch 
                println(e)
                error_file = criterion * "_opt_" * mode * "_" * type * ".txt"
                open(error_file,"a") do io
                    println(io, seed, " ", m, " ", mode, " : ", e)
                end
            end
        end
    end
elseif criterion == "AF"
    for m in [20,30]#[20,30,50,60,80,100,150,200,250,300] 
        for seed in 1:5
            @show m, seed
            try
                if mode == "Boscia"
                    optDesign.solve_opt(seed, m, criterion, corr=corr)
                elseif mode == "SCIP"
                    optDesign.solve_opt_scip(seed, m, criterion, corr)
                elseif mode == "FrankWolfe"
                    optDesign.solve_opt_frank_wolfe(seed, m, criterion, corr)
                elseif mode == "Hypatia"
                    optDesign.solve_opt_hypatia(seed, m ,criterion, corr)
                elseif mode == "Pajarito"
                    optDesign.solve_opt_pajarito(seed, m, criterion, corr)
                else
                    error("Invalid mode!")
                end
            catch 
                println(e)
                error_file = criterion * "_opt_" * mode * "_" * type * ".txt"
                open(error_file,"a") do io
                    println(io, seed, " ", m, " ", mode, " : ", e)
                end
            end
        end
    end
else
    error("Invalid criterion!")
end =#







