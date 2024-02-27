# Optimal Experiment Design with Boscia

The repository for the "Solving the Optimal Experiment Design Problem with Mixed-Integer Convex Methods" [paper](https://arxiv.org/abs/2312.11200).

The experiment scripts are in the ODWB folder. The models and set-ups for the different solvers can be found in the source folder in the ODWB folder. 

The models and solvers tested are [Boscia.jl](https://github.com/ZIB-IOL/Boscia.jl) with the classical problem formulation, [SCIP.jl](https://github.com/scipopt/SCIP.jl) using the epigraph formulation, [Pajarito.jl](https://github.com/jump-dev/Pajarito.jl) implementing the direct conic model and a custom solver based on this [paper](https://link.springer.com/article/10.1007/s11222-021-10043-5)

The experiments are set up to be run on the cluster of the Zuse Institute Berlin. If you want to test individual instances, you can use Julia's environmental variables and call the ```run_optimal_design.jl``` file directly.
