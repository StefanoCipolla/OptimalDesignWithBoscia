
using optDesign
using Printf
using LinearAlgebra
using DataFrames
using CSV

# for the prelimenary runs
#dimensions = [20,30,50,60,80,100,120,150,180,200]
#fracs = [4,6,8,10]

# for exclusive 
dimensions = [50,60,80,100,120]
fracs = [4,10]

for m in dimensions #[20,30,50,60,80,100,
    for k in fracs
        n = Int(floor(m/k))
        for seed in 1:5
            A, _, N, _ = optDesign.build_data(seed, m, n, false, false)
            eigs = eigvals(A'*A)
            maxeig= maximum(eigs)
            mineig = minimum(eigs)
            df = DataFrame(seed=seed, numberOfExperiments=m, numberOfParameters=n, frac=k, numberOfAllowedEx=N, eigmax = maxeig, eigmin=mineig, ratio = maxeig/mineig)
            file_name = joinpath(@__DIR__, "csv/opt_independent_data.csv")
            if !isfile(file_name) 
                CSV.write(file_name, df, append=true, writeheader=true, delim=";")
            else 
                CSV.write(file_name, df, append=true, delim=";")
            end

            A, _, N, _ = optDesign.build_data(seed, m, n, false, true)
            eigs = eigvals(A'*A)
            maxeig= maximum(eigs)
            mineig = minimum(eigs)
            df = DataFrame(seed=seed, numberOfExperiments=m, numberOfParameters=n, frac=k, numberOfAllowedEx=N, eigmax = maxeig, eigmin=mineig, ratio = maxeig/mineig)
            file_name = joinpath(@__DIR__, "csv/opt_correlated_data.csv")
            if !isfile(file_name) 
                CSV.write(file_name, df, append=true, writeheader=true, delim=";")
            else 
                CSV.write(file_name, df, append=true, delim=";")
            end

            A, _, N, _ = optDesign.build_data(seed, m, n, true, false)
            eigs = eigvals(A'*A)
            maxeig= maximum(eigs)
            mineig = minimum(eigs)
            df = DataFrame(seed=seed, numberOfExperiments=m, numberOfParameters=n, frac=k, numberOfAllowedEx=N, eigmax = maxeig, eigmin=mineig, ratio = maxeig/mineig)
            file_name = joinpath(@__DIR__, "csv/fusion_independent_data.csv")
            if !isfile(file_name) 
                CSV.write(file_name, df, append=true, writeheader=true, delim=";")
            else 
                CSV.write(file_name, df, append=true, delim=";")
            end

            A, _, N, _ = optDesign.build_data(seed, m, n, true, true)
            eigs = eigvals(A'*A)
            maxeig= maximum(eigs)
            mineig = minimum(eigs)
            df = DataFrame(seed=seed, numberOfExperiments=m, numberOfParameters=n, frac=k, numberOfAllowedEx=N, eigmax = maxeig, eigmin=mineig, ratio = maxeig/mineig)
            file_name = joinpath(@__DIR__, "csv/fusion_correlated_data.csv")
            if !isfile(file_name) 
                CSV.write(file_name, df, append=true, writeheader=true, delim=";")
            else 
                CSV.write(file_name, df, append=true, delim=";")
            end

        end
    end
end