# plot_results.jl
using StatsPlots
using CSV
using DataFrames
using Measures

    df = CSV.read("raw_timing_data_n10.csv", DataFrame)

    # Ensure m values are treated as categorical strings
    df.m_str = string.(df.m)

    Ns = unique(df.N)
    sort!(Ns)
    nrows = ceil(Int, length(Ns) / 3)

    plots = []

    for N_val in Ns
        df_N = filter(row -> row.N == N_val && !isnan(row.time), df)
        pal = palette(:auto, length(unique(df_N.m_str)))
        #p = @df df_N boxplot(:m_str, :time, group=:m_str,
            #xlabel="m", ylabel="Time (s)", title="N = $N_val", legend=false)

        p = @df df_N violin(:m_str, :time, group=:m_str, 
                            palette   = pal,
                            fillalpha = 0.9,    # slightly transparent violins
                            legend    = false)
        @df df_N boxplot!(p,:m_str, :time, group=:m_str, 
                            palette   = pal,
                            fillalpha = 0.2,    # slightly transparent violins
                            legend    = false)
        @df df_N dotplot!(p,:m_str, :time, group=:m_str, 
                          legend=false, marker=(:black, stroke(0)))
        xlabel!(p, "m")
        ylabel!(p, "Time (s)")
        title!(p, "N = $N_val",)
            
        push!(plots, p)
    end

    plt = plot(plots..., layout=(nrows, 3), 
               size=(1200, 300 * nrows), 
               bottom_margin = 7mm,
               left_margin  = 7mm)
    
    savefig(plt, "boxplot_computation_time_fixed_n10.png")
    display(plt)