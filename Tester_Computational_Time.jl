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

    ms = 100:100:500
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
                catch e
                    println("Error for m=$m, N=$N_val, trial=$trial: $e")
                    push!(all_data, (N=N_val, m=m, time=NaN))
                end
            end
        end
    end

    # Save all timing data
    CSV.write("raw_timing_data_n10.csv", all_data)


