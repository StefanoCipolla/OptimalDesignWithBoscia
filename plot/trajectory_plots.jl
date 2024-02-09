using LinearAlgebra
using Plots
#using PyPlot
using Random
using optDesign
seed = 30
pyplot()

function plot_contour(x, y, z, type, problem, func)
    PyPlot.matplotlib[:rc]("text", usetex=true)
    PyPlot.matplotlib[:rc]("font", size=12, family="cursive")
    PyPlot.matplotlib[:rc]("axes", labelsize=14)
    PyPlot.matplotlib[:rc]("text.latex", preamble=raw"""
    \usepackage{libertine}
    \usepackage{libertinust1math}
    """)

   # fig = plt.figure(figsize=(6.5,3.5))
   # ax = fig.add_subplot(111)
   Plots.text("abcde", "Arial")

    contour(x, y, z, xlabel="x-axis", ylabel="y-axis") #title="Contour Plot of f(x, y)"
    #xlabel("x-axis")
    #ylabel("y-axis")

    file = joinpath(@__DIR__, "../csv/plots/contour_plot_" * func * "_" * type * "_" * problem * ".pdf")
    Plots.savefig(file)
end

for fusion in [true, false]
    for corr in [true, false]
        Random.seed!(seed)
        #A = randn(10,2)
        A, _, _, _ = optDesign.build_data(seed, 10, 2, fusion, corr)
        type = corr ? "CORR" : "IND"
        problem = fusion ? "F" : "O"

        C = rand(10, 10)
        C = C * C'

        function inf_matrix(x,y,fusion)
            if fusion
                return Symmetric(C + A * Diagonal([x,y]) * A') 
            else
                return Symmetric(A * Diagonal([x,y]) * A') + 20I
            end
        end
        f(x, y) = -1 * logdet(inf_matrix(x, y, fusion))
        g(x, y) = tr(inv(inf_matrix(x, y, fusion)))
        h(x, y) = tr(inf_matrix(x, y, fusion))

        # Generate data points
        x = range(0.001, stop=10, length=100)
        y = range(0.001, stop=10, length=100)

        # Create a grid of (x, y) values
        grid_ = Iterators.product(x, y)

        #fig = plt.figure(figsize=(6.5,3.5))
        #ax = fig.add_subplot(111)

        # Evaluate the function at each grid point
        zf = [f(xi, yi) for (xi, yi) in grid_]
        zf = reshape(zf, length(x), length(y))
       # fig = contour(x, y, zf) #title="Contour Plot of f(x, y)"
       # xlabel("x-axis")
       # ylabel("y-axis")
       plot_contour(x, y, zf, type, problem, "log_det")

        #ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), fontsize=12, fancybox=true, shadow=false, ncol=2)
        #fig.tight_layout()
        #save_direc = corr ? "correlated/" : "independent/"
        #file = joinpath(@__DIR__, "../csv/plots/contour_plot_log_det_" * type * "_" * problem * ".pdf")
        #Plots.savefig(file)

        zg = [g(xi, yi) for (xi, yi) in grid_]
        zg = reshape(zg, length(x), length(y))
        #fig = contour(x, y, zg) #title="Contour Plot of g(x, y)"
        #xlabel("x-axis")
        #ylabel("y-axis")

        #file = joinpath(@__DIR__, "../csv/plots/contour_plot_trace_inverse_" * type * "_" * problem * ".pdf")
        #Plots.savefig(file)
        plot_contour(x, y, zg, type, problem, "trace_inverse")

        zh = [h(xi, yi) for (xi, yi) in grid_]
        zh = reshape(zh, length(x), length(y))
        #fig = contour(x, y, zh) #title="Contour Plot of h(x, y)"
        #xlabel("x-axis")
        #ylabel("y-axis")

        #file = joinpath(@__DIR__, "../csv/plots/contour_plot_trace_" * type * "_" * problem * ".pdf")
        #Plots.savefig(file)
        plot_contour(x, y, zh, type, problem, "trace")
    end
end
