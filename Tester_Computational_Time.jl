using ODWB
using StatsPlots
using Statistics
using CSV
using DataFrames


    seed = 1
    n = 10
    criterion = "D"
    time_limit = 36000
    corr = true
    long_run = false

    ms = 100:100:200
    Ns = n:5:20
    num_trials = 10
    nrows = ceil(Int, length(Ns) / 3)

    # Prepare to collect raw timing data for CSV and plotting
    all_data = DataFrame(N=Int[], m=Int[], time=Float64[])

    for N_val in Ns
        for m in ms
            for trial in 1:num_trials
                try
                    A, C, _, ub, C_hat = ODWB.build_data_ext(seed + trial, m, n, criterion in ["AF", "DF", "GTIF"], corr, scaling_C=long_run)
                    elapsed = @elapsed ODWB.solve_opt_ext(seed + trial, m, n, time_limit, criterion, corr;
                        A_external=A,
                        C_external=C,
                        N_external=float(N_val),
                        ub_external=ub,
                        C_hat_external=C_hat)
                    push!(all_data, (N=N_val, m=m, time=elapsed))

                     file_name = "/Users/stefanocipolla/Library/CloudStorage/OneDrive-UniversityofSouthampton/Southampton_Work/Data_Boscia/Boscia_generated" *string(n)* "_" * string(m) * "_" * string(N_val) * "_" * string(trial) * "_.m"
                     file = ODWB.matopen(file_name, "w")
                     ODWB.write(file, "X", A)
                     ODWB.write(file, "N", N)
                     ODWB.close(file)
                catch e
                    println("Error for m=$m, N=$N_val, trial=$trial: $e")
                    push!(all_data, (N=N_val, m=m, time=NaN))
                end
            end
        end
    end

    # Save all timing data
    file_name = "/Users/stefanocipolla/Library/CloudStorage/OneDrive-UniversityofSouthampton/Southampton_Work/OptimalDesignWithBoscia/Results/Boscia" * string(n) * "_.csv"
    CSV.write(file_name, all_data)


