using PyPlot
using DataFrames
using CSV
#using Plots

# function plot_boscia_vs_scip(example; boscia=true, scip_oa=false, ipopt=false, afw=true, ss=true, as=true, as_ss=true, boscia_methods=true)
function plot_termination(corr)
    criteria = ["AF","A","DF","D"]
    type = corr ? "correlated" : "independent"
    for criterion in criteria
    if  criterion == "AF" || criterion == "DF"
        df_full = DataFrame(CSV.File(joinpath(@__DIR__, "csv/" * criterion * "_optimality_"* type * "_non_grouped.csv")))
        df = select(df_full, [:ratio, :timeBoscia, :terminationBoscia, :timeScip, :terminationScip, :timePajarito, :terminationPajarito])

        df[df.timeBoscia.>1800, :timeBoscia] .= 1800
        df[df.timeScip.>1800, :timeScip] .= 1800
        df[df.timePajarito.>1800, :timePajarito] .= 1800
    elseif criterion == "D" || criterion == "A"
        df_full = DataFrame(CSV.File(joinpath(@__DIR__, "csv/" * criterion * "_optimality_"* type * "_non_grouped.csv")))
        df = select(df_full, [:ratio, :timeBoscia, :terminationBoscia, :timePajarito, :terminationPajarito])

        df[df.timeBoscia.>1800, :timeBoscia] .= 1800
        df[df.timePajarito.>1800, :timePajarito] .= 1800
    end

    time_limit = 1800

    colors = ["b", "m", "c", "r", "g", "y", "k", "peru"]
    markers = ["o", "s", "^", "P", "X", "H", "D"]
    linestyle = ["-", ":", "-.", "--"]

    PyPlot.matplotlib[:rc]("text", usetex=true)
    PyPlot.matplotlib[:rc]("font", size=12, family="cursive")
    PyPlot.matplotlib[:rc]("axes", labelsize=14)
    PyPlot.matplotlib[:rc]("text.latex", preamble=raw"""
    \usepackage{libertine}
    \usepackage{libertinust1math}
    """)

    #ax.scatter(df_temp[!,"dimension"], df_temp[!,"time_boscia"], label="BO (ours)", color=colors[1], marker=markers[1])
    #ax.scatter(df_temp[!,"dimension"], df_temp[!,"time_scip"], label="SCIP", color=colors[3], marker=markers[2])

    
    fig = plt.figure(figsize=(6.5,3.5))
    ax = fig.add_subplot(111)
    df_boscia = deepcopy(df)
    filter!(row -> !(row.terminationBoscia == 0),  df_boscia)
    time_boscia = sort(df_boscia[!,"timeBoscia"])
    push!(time_boscia, 1.1 * time_limit)
    ax.plot(time_boscia, [1:nrow(df_boscia); nrow(df_boscia)], label="Boscia", color=colors[1], marker=markers[1], markevery=0.1, linestyle=linestyle[1])
   
    if criterion == "AF" || criterion == "DF"
        df_scip = deepcopy(df)
        filter!(row -> !(row.terminationScip == 0),  df_scip)     
        time_scip = sort(df_scip[!,"timeScip"])
        push!(time_scip, 1.1 * time_limit)
        ax.plot(time_scip, [1:nrow(df_scip); nrow(df_scip)], label="SCIP+OA", color=colors[end], marker=markers[2], markevery=0.1, linestyle=linestyle[2])
    end

    df_pajarito = deepcopy(df)
    filter!(row -> !(row.terminationPajarito == 0),  df_pajarito)
    time_pajarito = sort(df_pajarito[!,"timePajarito"])
    push!(time_pajarito, 1.1 * time_limit)
    ax.plot(time_pajarito, [1:nrow(df_pajarito); nrow(df_pajarito)], label="Pajarito", color=colors[2], marker=markers[3], markevery=0.1, linestyle=linestyle[3])
    
    ylabel("Solved instances")
    #locator_params(axis="y", nbins=4)
    xlabel("Time (s)")
    ax.set_xscale("log")
    ax.grid()
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), fontsize=12,
    fancybox=true, shadow=false, ncol=2)
   
  #=  if criterion == "A"
        title("A-Opt: Solved instances over time" , loc="center")
    elseif criterion == "AF"
        title("A-Fusion: Solved instances over time" , loc="center")
    elseif criterion == "D"
        title("D-Opt: Solved instances over time" , loc="center")
    elseif criterion == "DF"
        title("D-Fusion: Solved instances over time" , loc="center")
    end=#


    fig.tight_layout()

    over_what = "time" 

    file = "csv/plots/" * criterion * "_" * type * "_" * over_what *  ".pdf"

    PyPlot.savefig(file)
    end
end

function plot_termination_ratio(corr, time=true, both=true)

    criteria = ["AF","A","DF","D"]
    for criterion in criteria
    if  criterion == "AF" || criterion == "DF"
        if both
            df_full = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/" * criterion * "_optimality_correlated_non_grouped.csv")))
            df = select(df_full, [:ratio, :timeBoscia, :terminationBoscia, :timeScip, :terminationScip, :timePajarito, :terminationPajarito, :timeCustomBB, :terminationCustomBB])

            df[df.timeBoscia.>3600, :timeBoscia] .= 3600
            df[df.timeScip.>3600, :timeScip] .= 3600
            df[df.timePajarito.>3600, :timePajarito] .= 3600
            df[df.timeCustomBB.>3600, :timeCustomBB] .= 3600

            df_full_ind = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/" * criterion * "_optimality_independent_non_grouped.csv")))
            df_ind = select(df_full_ind, [:ratio, :timeBoscia, :terminationBoscia, :timeScip, :terminationScip, :timePajarito, :terminationPajarito, :timeCustomBB, :terminationCustomBB])
            df_ind[df_ind.timeBoscia.>3600, :timeBoscia] .= 3600
            df_ind[df_ind.timeScip.>3600, :timeScip] .= 3600
            df_ind[df_ind.timePajarito.>3600, :timePajarito] .= 3600
            df_ind[df_ind.timeCustomBB.>3600, :timeCustomBB] .= 3600
        elseif corr 
            df_full = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/" * criterion * "_optimality_correlated_non_grouped.csv")))
            df = select(df_full, [:ratio, :timeBoscia, :terminationBoscia, :timeScip, :terminationScip, :timePajarito, :terminationPajarito])

            df[df.timeBoscia.>3600, :timeBoscia] .= 3600
            df[df.timeScip.>3600, :timeScip] .= 3600
            df[df.timePajarito.>3600, :timePajarito] .= 3600
        else
            df_full = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/" * criterion * "_optimality_independent_non_grouped.csv")))
            df = select(df_full, [:ratio, :timeBoscia, :terminationBoscia, :timeScip, :terminationScip, :timePajarito, :terminationPajarito])
            df[df.timeBoscia.>3600, :timeBoscia] .= 3600
            df[df.timeScip.>3600, :timeScip] .= 3600
            df[df.timePajarito.>3600, :timePajarito] .= 3600
        end
    elseif criterion == "D" || criterion == "A"
        if both
            df_full = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/" * criterion * "_optimality_correlated_non_grouped.csv")))
            df = select(df_full, [:ratio, :timeBoscia, :terminationBoscia, :timePajarito, :terminationPajarito, :timeCustomBB, :terminationCustomBB])

            df[df.timeBoscia.>3600, :timeBoscia] .= 3600
            df[df.timePajarito.>3600, :timePajarito] .= 3600
            df[df.timeCustomBB.>3600,:timeCustomBB] .= 3600

            df_full_ind = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/" * criterion * "_optimality_independent_non_grouped.csv")))
            df_ind = select(df_full_ind, [:ratio, :timeBoscia, :terminationBoscia, :timePajarito, :terminationPajarito, :timeCustomBB, :terminationCustomBB])

            df_ind[df_ind.timeBoscia.>3600, :timeBoscia] .= 3600
            df_ind[df_ind.timePajarito.>3600, :timePajarito] .= 3600
            df_ind[df.timeCustomBB.>3600,:timeCustomBB] .= 3600
        elseif corr 
            df_full = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/" * criterion * "_optimality_correlated_non_grouped.csv")))
            df = select(df_full, [:ratio, :timeBoscia, :terminationBoscia, :timePajarito, :terminationPajarito])

            df[df.timeBoscia.>3600, :timeBoscia] .= 3600
            df[df.timePajarito.>3600, :timePajarito] .= 3600
            df[df.timeCustomBB.>3600,:timeCustomBB] .= 3600
        else
            df_full = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/" * criterion * "_optimality_independent_non_grouped.csv")))
            df = select(df_full, [:ratio, :timeBoscia, :terminationBoscia, :timePajarito, :terminationPajarito])

            df[df.timeBoscia.>3600, :timeBoscia] .= 3600
            df[df.timePajarito.>3600, :timePajarito] .= 3600
            df[df.timeCustomBB.>3600,:timeCustomBB] .= 3600
        end
    end

    colors = ["b", "m", "c", "r", "g", "y", "k", "peru"]
    markers = ["o", "s", "^", "P", "X", "H", "D"]
    linestyle = ["-", ":", "-.", "--"]

    PyPlot.matplotlib[:rc]("text", usetex=true)
    PyPlot.matplotlib[:rc]("font", size=9, family="cursive")
    PyPlot.matplotlib[:rc]("axes", labelsize=12)
    PyPlot.matplotlib[:rc]("text.latex", preamble=raw"""
    \usepackage{libertine}
    \usepackage{libertinust1math}
    """)

    if criterion in ["A","D"]
        fig, axs = plt.subplots(3, sharex=true, sharey=true, figsize=(6.5,4.5))
    else
        fig, axs = plt.subplots(4, sharex=true, sharey=true, figsize=(6.5,6.5))
    end

    linewidth = 2

    df_boscia_ind = deepcopy(df_ind)
    filter!(row -> !(row.terminationBoscia == 0),  df_boscia_ind)
    x_boscia_ind = time ? sort(df_boscia_ind[!,"timeBoscia"]) : sort(df_boscia_ind[!,"ratio"])
    Boscia_plot = axs[1].plot(x_boscia_ind, 1:nrow(df_boscia_ind), label="Boscia independent", color=colors[1], marker=markers[1], markevery=0.1, linestyle=linestyle[2], linewidth=linewidth)

    df_boscia = deepcopy(df)
    filter!(row -> !(row.terminationBoscia == 0),  df_boscia)
    x_boscia = time ? sort(df_boscia[!,"timeBoscia"]) : sort(df_boscia[!,"ratio"])
    Boscia_plot = axs[1].plot(x_boscia, 1:nrow(df_boscia), label="Boscia correlated", color=colors[6], marker=markers[4], markevery=0.1, linestyle=linestyle[1], linewidth=linewidth)
    axs[1].grid()

   df_pajarito_ind = deepcopy(df_ind)
    filter!(row -> !(row.terminationPajarito == 0),  df_pajarito_ind)
    x_pajarito_ind = time ? sort(df_pajarito_ind[!,"timePajarito"]) : sort(df_pajarito_ind[!,"ratio"])
    Pajarito_plot= axs[2].plot(x_pajarito_ind, 1:nrow(df_pajarito_ind), label="Pajarito independent", color=colors[2], marker=markers[3], markevery=0.1, linestyle=linestyle[2], linewidth=linewidth) 

    df_pajarito = deepcopy(df)
    filter!(row -> !(row.terminationPajarito == 0),  df_pajarito)
    x_pajarito = time ? sort(df_pajarito[!,"timePajarito"]) : sort(df_pajarito[!,"ratio"])
    Pajarito_plot= axs[2].plot(x_pajarito, 1:nrow(df_pajarito), label="Pajarito correlated", color=colors[7], marker=markers[6], markevery=0.1, linestyle=linestyle[2], linewidth=linewidth)
    axs[2].grid()

    
    if criterion == "AF" || criterion == "DF"
       df_scip_ind = deepcopy(df_ind)
        filter!(row -> !(row.terminationScip == 0),  df_scip_ind)   
        x_scip_ind = time ? sort(df_scip_ind[!, "timeScip"]) : sort(df_scip_ind[!,"ratio"])  
        SCIP_plot =axs[3].plot(x_scip_ind, 1:nrow(df_scip_ind), label="SCIP+OA independent", color=colors[end], marker=markers[2], markevery=0.1, linestyle=linestyle[2], linewidth=linewidth) 

        df_scip = deepcopy(df)
        filter!(row -> !(row.terminationScip == 0),  df_scip)   
        x_scip = time ? sort(df_scip[!,"timeScip"]) : sort(df_scip[!,"ratio"])  
        SCIP_plot =axs[3].plot(x_scip, 1:nrow(df_scip), label="SCIP+OA correlated", color=colors[5], marker=markers[7], markevery=0.1, linestyle=linestyle[4],linewidth=linewidth)

        axs[3].grid()
        axs[3].legend(loc="lower right")#, bbox_to_anchor=(0.5, -0.3), fontsize=12,fancybox=true, shadow=false, ncol=2) 
    end   
   # else
   k = criterion in ["AF","DF"] ? 4 : 3
        df_bb_ind = deepcopy(df_ind)
        filter!(row -> !(row.terminationCustomBB == 0),  df_bb_ind)   
        x_bb_ind = time ? sort(df_bb_ind[!, "timeCustomBB"]) : sort(df_bb_ind[!,"ratio"])  
        BB_plot =axs[k].plot(x_bb_ind, 1:nrow(df_bb_ind), label="Custom BnB independent", color=colors[4], marker=markers[4], markevery=0.1, linestyle=linestyle[2], linewidth=linewidth) 

        df_bb = deepcopy(df)
        filter!(row -> !(row.terminationCustomBB == 0),  df_bb)   
        x_bb = time ? sort(df_bb[!,"timeCustomBB"]) : sort(df_bb[!,"ratio"])  
        SCIP_plot =axs[k].plot(x_bb, 1:nrow(df_bb), label="Custom BnB correlated", color=colors[3], marker=markers[5], markevery=0.1, linestyle=linestyle[3], linewidth=linewidth)

        axs[k].grid()
        axs[k].legend(loc="lower right")#, bbox_to_anchor=(0.5, -0.3), fontsize=12,fancybox=true, shadow=false, ncol=2) 
   # end

    #plot(Boscia_plot, Pajarito_plot, title="Overall title of the plot")
   # ylabel("Solved instances", loc="top")
    axs[2].set_ylabel("Solved instances", loc="center")
    #locator_params(axis="y", nbins=4)
    if time
        xlabel("Time")
    else
        xlabel("Ratio of eigenvalues")
    end

    axs[1].legend(loc="lower right")#, bbox_to_anchor=(0.5, -0.3), fontsize=12, fancybox=true, shadow=false, ncol=2)
    axs[2].legend(loc="lower right")#, bbox_to_anchor=(0.5, -0.3), fontsize=12,fancybox=true, shadow=false, ncol=2)
    #=
    if criterion == "A"
        PyPlot.suptitle("A-Opt: Solved instances over the ratio of eigenvalues")
    elseif criterion == "AF"
        PyPlot.suptitle("A-Fusion: Solved instances over the ratio of eigenvalues")
    elseif criterion == "D"
        PyPlot.suptitle("D-Opt: Solved instances over the ratio of eigenvalues")
    elseif criterion == "DF"
        PyPlot.suptitle("D-Fusion: Solved instances over the ratio of eigenvalues")
    end=#

    fig.tight_layout()

    type  = corr ? "correlated" : "independent"
    by_what = time ? "time" : "ratio"

    file = ""
    if both
        file = "csv/plots/" * criterion * "_" * by_what * ".pdf"
    else
        file = "csv/plots/" * criterion * "_"  * type * "_"  * by_what * ".pdf"
    end

    PyPlot.savefig(file)
    end
end

"""
Plot for the curavture constant
"""
function plot_curvature()

    colors = ["b", "m", "c", "r", "g", "y", "k", "peru"]
    markers = ["o", "s", "^", "P", "X", "H", "D"]
    linestyle = ["-", ":", "-.", "--"]

    PyPlot.matplotlib[:rc]("text", usetex=true)
    PyPlot.matplotlib[:rc]("font", size=12, family="cursive")
    PyPlot.matplotlib[:rc]("axes", labelsize=14)
    PyPlot.matplotlib[:rc]("text.latex", preamble=raw"""
    \usepackage{libertine}
    \usepackage{libertinust1math}
    """)

    fig = plt.figure(figsize=(6.5,6.5)) #figsize=(6.5,6.5)
    ax = fig.add_subplot(111)

    function quadratic(a,x)
        return a*x.^2
    end

    # lb, time 
    x = collect(-5:0.5:5)
    ax.plot(x, quadratic(0.1,x), label="a=0.1", color=colors[1], marker=markers[1], markevery=0.1, alpha=.5, linestyle=linestyle[2])
    ax.plot(x, quadratic(0.7,x), label="a=0.7", color=colors[end], marker=markers[2], markevery=0.1, alpha=.5, linestyle=linestyle[3])
    ax.plot(x, quadratic(1.0,x), label="a=1.0", color=colors[2], marker=markers[3], markevery=0.1, alpha=.5, linestyle=linestyle[1])
    ax.plot(x, quadratic(1.5,x), label="a=1.5", color=colors[3], marker=markers[4], markevery=0.1, alpha=.5, linestyle=linestyle[4])
    ax.plot(x, quadratic(4.0,x), label="a=4.0", color=colors[4], marker=markers[5], markevery=0.1, alpha=.5, linestyle=linestyle[2])
    ax.plot(x, quadratic(10.0,x), label="a=10.0", color=colors[5], marker=markers[6], markevery=0.1, alpha=.5, linestyle=linestyle[3])

    ylabel("y")
    #locator_params(axis="y", nbins=4)
    xlabel("x")
    ax.grid()

   ax.legend(loc="upper center", fontsize=12, #, bbox_to_anchor=(0.5, -0.2)
   fancybox=true, shadow=false, ncol=3)
    fig.tight_layout()
    file = joinpath(@__DIR__, "csv/plots/example_curvature.pdf")
    PyPlot.savefig(file)
end

function plot_progress(corr, criterion)

    type = corr ? "correlated" : "independent"
    direc = corr ? "correlated/" : "linearly_independent/"
    #criteria = ["A","AF","D","DF"]


    colors = ["b", "m", "c", "r", "g", "y", "k", "peru"]
    markers = ["o", "s", "^", "P", "X", "H", "D"]
    linestyle = ["-", ":", "-.", "--"]

    PyPlot.matplotlib[:rc]("text", usetex=true)
    PyPlot.matplotlib[:rc]("font", size=12, family="cursive")
    PyPlot.matplotlib[:rc]("axes", labelsize=14)
    PyPlot.matplotlib[:rc]("text.latex", preamble=raw"""
    \usepackage{libertine}
    \usepackage{libertinust1math}
    """)

    #for criterion in criteria
        for dim in [50,60,80,100,120]
            for k in [4,10]
                for seed in 1:5
                    n = Int(floor(dim/k))
                    df = DataFrame(CSV.File(joinpath(@__DIR__, "csv/full_runs_boscia/boscia_" *  criterion * "_optimality_" * type * "_" * string(dim) * "-" * string(n) * "_" * string(seed) * ".csv")))

                    #  df[df.time.>1800, :time] .= 1800

                    df[!, :time] = df[!,:time]./1000
                
                    fig = plt.figure(figsize=(6.5,3.5))
                    ax = fig.add_subplot(111)
            
                    # lb, time 
                    ax.plot(df[!,"time"], df[!,"lowerBound"], label="Lower bound", color=colors[1], marker=markers[1], markevery=0.1, alpha=.5, linestyle=linestyle[3])
                    ax.plot(df[!,"time"], df[!,"upperBound"], label="Incumbent", color=colors[end], marker=markers[2], markevery=0.1, alpha=.5, linestyle=linestyle[1])
            
                    ylabel("Lower bound")
                    #locator_params(axis="y", nbins=4)
                    xlabel("Time (s)")
                    # ax.set_xscale("log")
                    ax.grid()

                    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), fontsize=12,
                    fancybox=true, shadow=false, ncol=2)
                    fig.tight_layout()
                    save_direc = corr ? "correlated/" : "independent/"
                    file = joinpath(@__DIR__, "csv/plots/Progress_Plots/" * criterion * "/"  * criterion * "-" * type * "_" * string(dim) * "_" * string(n) * "_" * string(seed) * ".pdf")
                    PyPlot.savefig(file)
                end
            end
        end
   # end

    return true
end

function plot_dual_gap(corr, criterion, seed, m, n)

    type = corr ? "correlated" : "independent"
    direc = corr ? "correlated/" : "linearly_independent/"
    #criteria = ["A","AF","D","DF"]


    colors = ["b", "m", "c", "r", "g", "y", "k", "peru"]
    markers = ["o", "s", "^", "P", "X", "H", "D"]
    linestyle = ["-", ":", "-.", "--"]

    PyPlot.matplotlib[:rc]("text", usetex=true)
    PyPlot.matplotlib[:rc]("font", size=12, family="cursive")
    PyPlot.matplotlib[:rc]("axes", labelsize=14)
    PyPlot.matplotlib[:rc]("text.latex", preamble=raw"""
    \usepackage{libertine}
    \usepackage{libertinust1math}
    """)

    #for criterion in criteria
      #  for dim in [50,60,80,100,120]
      #      for k in [4,10]
      #          for seed in 1:5
                 #   n = Int(floor(dim/k))
                    df = DataFrame(CSV.File(joinpath(@__DIR__, "csv/full_runs_boscia/boscia_" *  criterion * "_optimality_" * type * "_" * string(m) * "-" * string(n) * "_" * string(seed) * ".csv")))

                    #  df[df.time.>1800, :time] .= 1800

                    df[!, :time] = df[!,:time]./1000
                
                    fig = plt.figure(figsize=(6.5,3.5))
                    ax = fig.add_subplot(111)
            
                    # lb, time 
                    ax.plot(df[!,"time"], df[!,"upperBound"] - df[!,"lowerBound"], label="Dual Gap", color=colors[2], marker=markers[3], markevery=0.1, alpha=.5, linestyle=linestyle[1])
                   # ax.plot(df[!,"time"], df[!,"upperBound"], label="Incumbent", color=colors[end], marker=markers[2], markevery=0.1, alpha=.5, linestyle=linestyle[1])
            
                    ylabel("Dual Gap")
                    #locator_params(axis="y", nbins=4)
                    xlabel("Time (s)")
                    # ax.set_xscale("log")
                    ax.grid()

                    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), fontsize=12, fancybox=true, shadow=false, ncol=2)
                    fig.tight_layout()
                    file = joinpath(@__DIR__, "csv/plots/Dual_Gap_Plots/"  * criterion * "-" * type * "_" * string(m) * "_" * string(n) * "_" * string(seed) * ".pdf")
                    PyPlot.savefig(file)
       #         end
       #     end
       # end
   # end

    return true
end

function plot_progress_dual_gap(corr, criterion, seed, m, n)

    type = corr ? "correlated" : "independent"
    direc = corr ? "correlated/" : "linearly_independent/"
    #criteria = ["A","AF","D","DF"]


    colors = ["b", "m", "c", "r", "g", "y", "k", "peru"]
    markers = ["o", "s", "^", "P", "X", "H", "D"]
    linestyle = ["-", ":", "-.", "--"]

    PyPlot.matplotlib[:rc]("text", usetex=true)
    PyPlot.matplotlib[:rc]("font", size=12, family="cursive")
    PyPlot.matplotlib[:rc]("axes", labelsize=14)
    PyPlot.matplotlib[:rc]("text.latex", preamble=raw"""
    \usepackage{libertine}
    \usepackage{libertinust1math}
    """)


    df = DataFrame(CSV.File(joinpath(@__DIR__, "csv/full_runs_boscia/boscia_" *  criterion * "_optimality_" * type * "_" * string(m) * "-" * string(n) * "_" * string(seed) * ".csv")))

    #  df[df.time.>1800, :time] .= 1800

    df[!, :time] = df[!,:time]./1000

    fig, axs = plt.subplots(1,2, sharex=false, sharey=false, figsize=(13.0,4.5))
                
    #fig = plt.figure(figsize=(13.0,3.5))
    #axs = fig.add_subplot(112)

    axs[1].plot(df[!,"time"], df[!,"lowerBound"], label="Lower bound", color=colors[1], marker=markers[1], markevery=0.1, alpha=.5, linestyle=linestyle[3])
    axs[1].plot(df[!,"time"], df[!,"upperBound"], label="Incumbent", color=colors[end], marker=markers[2], markevery=0.1, alpha=.5, linestyle=linestyle[1])
            
    axs[1].set_ylabel("Lower bound")
    #locator_params(axis="y", nbins=4)
    axs[1].set_xlabel("Time (s)")
   # axs[1].legend("upper center")
    axs[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), fontsize=12, fancybox=true, shadow=false, ncol=2)
    axs[1].grid()
            
    # lb, time 
    axs[2].plot(df[!,"time"], df[!,"upperBound"] - df[!,"lowerBound"], label="Dual Gap", color=colors[2], marker=markers[3], markevery=0.1, alpha=.5, linestyle=linestyle[1])
    # ax.plot(df[!,"time"], df[!,"upperBound"], label="Incumbent", color=colors[end], marker=markers[2], markevery=0.1, alpha=.5, linestyle=linestyle[1])
            
    axs[2].set_ylabel("Dual Gap")
    #locator_params(axis="y", nbins=4)
    axs[2].set_xlabel("Time (s)")
    #axs[2].legend("upper center")
    axs[2].legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), fontsize=12, fancybox=true, shadow=false, ncol=2)
    #ax.set_xscale("log")
    axs[2].grid()

   # ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), fontsize=12, fancybox=true, shadow=false, ncol=2)
    fig.tight_layout()
    file = joinpath(@__DIR__, "csv/plots/Dual_Gap_Plots/progress_dual_"  * criterion * "-" * type * "_" * string(m) * "_" * string(n) * "_" * string(seed) * ".pdf")
    PyPlot.savefig(file)

#=    n = Int(floor(dim/k))
                    df = DataFrame(CSV.File(joinpath(@__DIR__, "csv/full_runs_boscia/boscia_" *  criterion * "_optimality_" * type * "_" * string(dim) * "-" * string(n) * "_" * string(seed) * ".csv")))

                    #  df[df.time.>1800, :time] .= 1800

                    df[!, :time] = df[!,:time]./1000
                
                    fig = plt.figure(figsize=(6.5,3.5))
                    ax = fig.add_subplot(111)
            
                    # lb, time 
                    ax.plot(df[!,"time"], df[!,"lowerBound"], label="Lower bound", color=colors[1], marker=markers[1], markevery=0.1, alpha=.5, linestyle=linestyle[3])
                    ax.plot(df[!,"time"], df[!,"upperBound"], label="Incumbent", color=colors[end], marker=markers[2], markevery=0.1, alpha=.5, linestyle=linestyle[1])
            
                    ylabel("Lower bound")
                    #locator_params(axis="y", nbins=4)
                    xlabel("Time (s)")
                    # ax.set_xscale("log")
                    ax.grid()

                    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), fontsize=12,
                    fancybox=true, shadow=false, ncol=2)
                    fig.tight_layout()
                    save_direc = corr ? "correlated/" : "independent/"
                    file = joinpath(@__DIR__, "csv/plots/Progress_Plots/" * criterion * "/"  * criterion * "-" * type * "_" * string(dim) * "_" * string(n) * "_" * string(seed) * ".pdf")
                    PyPlot.savefig(file)=#

    return true
end
