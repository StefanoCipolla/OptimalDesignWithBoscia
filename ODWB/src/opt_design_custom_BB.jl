"""
    Custom Branch-and-Bound for the Exact Optimal Design Problem
"""

"""
Node struct for the tree. Knows its bounds, the solution vector, number of iterations and time taken.
"""
mutable struct CustomBBNode <: Bonobo.AbstractNode
    std::Bonobo.BnBNodeInfo
    lower_bounds::Vector{Float64}
    upper_bounds::Vector{Float64}
    solution::Vector{Float64}
    time::Float64
    iteration_count::Int64
end

mutable struct ConstrainedBoxMIProblem{F,G,L}
    f::F
    grad::G
    linesearch::L
    nvars::Int64
    integer_vars::Vector{Int64}
    N::Float64
    time_limit::Union{Int64, Float64}
    rel_gap::Float64
    abs_gap::Float64
    Stage::Boscia.Solve_Stage
end

ConstrainedBoxMIProblem(f, g, linesearch, n, int_vars, N, time_limit) =
    ConstrainedBoxMIProblem(f, g, linesearch, n, int_vars, N, time_limit, 1e-2, 1e-6, Boscia.SOLVING)

    """
Returns the solution vector of the relaxed node problem.
"""
function Bonobo.get_relaxed_values(tree::Bonobo.BnBTree{<:CustomBBNode}, node)
    return node.solution
end

"""
Create the information for the childern nodes.
"""
function Bonobo.get_branching_nodes_info(tree::Bonobo.BnBTree{<:CustomBBNode}, node, vidx)
    x = Bonobo.get_relaxed_values(tree, node)
    @assert node.lower_bounds[vidx] < node.upper_bounds[vidx]
    frac_val = x[vidx]

    # Left child keeps the lower bounds and gets a new upper bound at vidx.
    left_bounds = copy(node.upper_bounds)
    left_bounds[vidx] = floor(frac_val)
    left_solution = copy(node.solution)
    left_solution[vidx] = floor(frac_val)


    node_info_left = (lower_bounds=node.lower_bounds, 
                    upper_bounds = left_bounds,
                    solution = left_solution,
                    time = 0.0,
                    iteration_count = 0
    )

    # Right child keeps the upper bounds and gets a new lower bound at vidx.
    right_bounds = copy(node.lower_bounds)
    right_bounds[vidx] = ceil(frac_val)
    right_solution = copy(node.solution)
    right_solution[vidx] = ceil(frac_val)

    node_info_right = (lower_bounds=right_bounds, 
                    upper_bounds = node.upper_bounds,
                    solution = right_solution,
                    time = 0.0,
                    iteration_count = 0
    )

    return node_info_left, node_info_right
end

"""
Solve the box constrained approximate optimal design problem.
"""
function Bonobo.evaluate_node!(tree, node)
    # count time and iterations
    time_ref = Dates.now()
    iter_count = 0

    # get feasible start point
    x = copy(node.solution)
    x = project_onto_feasible_region(x, node, tree.root.N)

    # If the node is not reachable, i.e. we can never satisfy the knapsack add_constraint
    if x === nothing
        return NaN, NaN
    end

    @assert check_feasibilty(x, node, tree)

    gradient = similar(x)
    tree.root.grad(gradient, x)
    linesearch_workspace = FrankWolfe.build_linesearch_workspace(FrankWolfe.Adaptive(), x, gradient)

    jdx = find_max_gradient(x, -gradient, node.upper_bounds)
    kdx = find_min_gradient(x, -gradient, node.lower_bounds)

    while gradient[jdx]/gradient[kdx] - 1 > 1e-5 && iter_count < 10000
        # find theta
        d = zeros(length(x))
        d[jdx] = 1.0
        d[kdx] = -1.0
        theta_max = min(node.upper_bounds[jdx] - x[jdx], x[kdx] - node.lower_bounds[kdx])
        theta =  tree.root.linesearch(
            node,
            tree.root.f,
            tree.root.grad,
            gradient, 
            x,
            d,
            theta_max,
            linesearch_workspace,
            iter_count,
            jdx,
            kdx
        )

        # weight shift
        x[jdx] += theta
        x[kdx] -= theta

        # compute the new jdx and kdx
        tree.root.grad(gradient, x)
        jdx = find_max_gradient(x, -gradient, node.upper_bounds)
        kdx = find_min_gradient(x, -gradient, node.lower_bounds)
        iter_count += 1

        if !check_feasibilty(x, node, tree)
            @show node.upper_bounds - x
            @show x - node.lower_bounds
            @show sum(x)
        end
        @assert check_feasibilty(x, node, tree)
    end

    time = float(Dates.value(Dates.now() - time_ref))
    node.time = time
    node.iteration_count = iter_count

    obj_value = tree.root.f(x)
    node.solution = x


    if isapprox(sum(isapprox.(x, round.(x); atol=1e-6, rtol=5e-2)), tree.root.nvars)  && check_feasibilty(round.(x), node, tree)
        #println("Integer solution found.")
        node.solution = round.(x)
        primal = tree.root.f(node.solution)
        return obj_value, primal
    end

    return obj_value, NaN
end

"""
Project point on the feasible region which is the intersection of the hyperbox
and the probability simplex scaled by integer N.
"""
function project_onto_feasible_region(x, node, N)
    # the difference in weight you need to distribute
    s_prod = sum(x) - N
    # we are already on the probability simplex
    if s_prod == 0
        return x
    end
    # If the difference is negative, we might crash the upper bounds. 
    # Else, we could crash the lower bounds.
    bounds = s_prod > 0 ? node.lower_bounds : node.upper_bounds
    b = (x - s_prod*ones(length(x))-bounds)*sign(s_prod)
    set = findall(x -> x > 0, b)

    # If the set is empty, we reached an infeasible node!
    if isempty(set)
        return nothing
    end

    dir = zeros(length(x))
    dir[set] .= 1
    x = x - s_prod/sum(dir) * dir

    @assert isapprox(sum(x), N; atol=1e-4, rtol=1e-2)
    @assert sum(node.upper_bounds - x .>= 0) == length(x)
    @assert sum(x- node.lower_bounds .>= 0) == length(x)

    return x
end

"""
Check feasibility of a given point.
"""
function check_feasibilty(x, lb, ub, N)
    m = length(x)
    return sum(x.>=lb) == m && sum(x.<=ub) == m && isapprox(sum(x), N) 
end
check_feasibilty(x, node, tree) = check_feasibilty(x, node.lower_bounds, node.upper_bounds, tree.root.N)

"""
Find maximum gradient entry where the upper bound has yet been met.
"""
function find_max_gradient(x, gradient, ub)
    bool_vec = x .< ub
    idx_set = findall(x -> x == 1, bool_vec) 
    return findfirst(x -> x == maximum(gradient[idx_set]), gradient)
end

"""
Find minimum gradient entry where the lower bound has yet been met.
"""
function find_min_gradient(x, gradient, lb)
    bool_vec = x .> lb
    idx_set = findall(x -> x == 1, bool_vec) 
    return findfirst(x -> x == minimum(gradient[idx_set]), gradient)
end

"""
Returns the list of indices corresponding to the integral variables.
"""
function Bonobo.get_branching_indices(root::ConstrainedBoxMIProblem)
    return root.integer_vars
end


"""
Set up the problem and solve it with the custom BB algorithm.
"""
function solve_opt_custom(seed, m, n, time_limit, criterion, corr; write = true, verbose= true)
    # build data
    A, C, N, ub = build_data(seed, m, n, criterion in ["AF","DF"], corr)
    # build function
    p = 1
    if criterion == "A" || criterion == "AF"
        p = -1
    elseif criterion == "D" || criterion == "DF"
        p = 0
    end
    @assert p <= 0
    if criterion in ["A","D"]
        C = nothing
    end
    f, grad!, linesearch = build_matrix_means_objective(A, p,C=C)

    # create problem and tree
    if criterion in ["A","D"]
        x0 = greedy_incumbent(A, m, n, N, ub)
    else
        x0= greedy_incumbent_fusion(A,m,n,N,ub)
    end
    root = ConstrainedBoxMIProblem(f, grad!, linesearch, m, collect(1:m), N, time_limit)
    nodeExample = CustomBBNode(
        Bonobo.BnBNodeInfo(1, 0.0, 0.0),
        zeros(m),
        ub, 
        fill(-1.0, length(x0)),
        0.0,
        0
    )

    Node = typeof(nodeExample)
    Value = Vector{Float64}
    tree = Bonobo.initialize(;
        traverse_strategy = Bonobo.BFS(),
        Node = Node,
        Solution = Bonobo.DefaultSolution{Node, Value},
        root = root,
        sense = :Min
    )

    Bonobo.set_root!(tree, (
        lower_bounds = zeros(length(x0)),
        upper_bounds = ub,
        solution = x0,
        time=0.0,
        iteration_count=0
    ))

    time_ref = Dates.now()
    callback = build_callback(tree, time_ref, verbose)

    # dummy solution - In case the process stops because of the time limit and no solution as been found yet.
    dummy_solution = Bonobo.DefaultSolution(Inf, fill(-1.0, length(x0)), nodeExample)
    if isapprox(sum(isapprox.(x0, round.(x0))), tree.root.nvars)
        # Minus f(x0) because the tree internally treats this as a min problem
        first_solution = Bonobo.DefaultSolution(tree.root.f(x0), x0, nodeExample)
        push!(tree.solutions, first_solution)
        tree.incumbent_solution = first_solution
        tree.incumbent = tree.root.f(x0)
    else
        push!(tree.solutions, dummy_solution)
        tree.incumbent_solution = dummy_solution
    end

    # Optimize
    sol_data = @timed Bonobo.optimize!(tree,callback=callback)

    time = sol_data.time
    # post processing 
    sol_x = convert.(Int, Bonobo.get_solution(tree))
    if tree.root.Stage == Boscia.OPT_TREE_EMPTY
        @assert sum(sol_x .!= dummy_solution.solution) == length(sol_x)
    end
    @show sol_x
    solution = f(sol_x)
    @show solution

    # check feasibility
    feasible = isfeasible(seed, m, n, criterion, sol_x, corr, N=N)
    if !feasible
        @show m, sum(ub - sol_x .>= 0)
        @show m, sum(sol_x .>= 0)
        @show N, sum(sol_x)
    end
    @assert feasible

    # Calculate objective with respect to the criteria used in Boscia
    if criterion in ["A","AF"]
        f_check, _ = build_a_criterion(A, criterion == "AF", C=C)
    else
        f_check, _ = build_d_criterion(A, criterion == "DF", C=C)
    end
    solution_scaled = f_check(sol_x)
    #solution_scaled = criterion == "A" ? exp(solution) : solution - n*log(n)
    @show solution_scaled

    # Check the solving stage
    @assert tree.root.Stage != Boscia.SOLVING
    status = if tree.root.Stage in [Boscia.OPT_GAP_REACHED, Boscia.OPT_TREE_EMPTY]
        "optimal"
    else
        "Time limit reached"
    end

    if write 
        # write data into file
        time_per_nodes = time/tree.num_nodes
        type = corr ? "correlated" : "independent"
        df = DataFrame(seed=seed, numberOfExperiments=m, numberOfParameters=n, N=N, time=time, time_per_nodes=time_per_nodes, solution=solution, solution_scaled=solution_scaled, number_nodes=tree.num_nodes, termination=status)
        file_name = "/home/htc/dhendryc/research_projects/MasterThesis/optDesign/csv/CustomBB/customBB_" * criterion * "_" * string(m) * "_" * type * "_optimality.csv"
        if !isfile(file_name)
            CSV.write(file_name, df, append=true, writeheader=true)
        else 
            CSV.write(file_name, df, append=true)
        end
    end
    return sol_x
end

function build_callback(tree, time_ref, verbose)
    return function callback(tree, node; node_infeasible = false, worse_than_incumbent = false)

        time = float(Dates.value(Dates.now() - time_ref))
        if (node.id == 1 || mod(tree.num_nodes, 100) == 0 ) && verbose
            println("ID: $(node.id) num_nodes: $(tree.num_nodes) LB: $(tree.lb) Incumbent: $(tree.incumbent) Total Time in s: $(time/1000.0) Time node in s: $(node.time/1000.0) Iterations: $(node.iteration_count)")
            if node_infeasible
                println("Node not feasible!")
            elseif worse_than_incumbent
                println("Node cut because it is worse than the incumbent")
            end
        end

        # break if time is met
        if tree.root.time_limit < Inf
            if time / 1000.0 ≥ tree.root.time_limit
                @assert tree.root.Stage == Boscia.SOLVING
                tree.root.Stage = Boscia.TIME_LIMIT_REACHED
            end
        end
    end
end

"""
Checks if the branch and bound can be stopped.
By default (in Bonobo) stops then the priority queue is empty. 
"""
function Bonobo.terminated(tree::Bonobo.BnBTree{<:CustomBBNode})
    # time limit reached
    if tree.root.Stage == Boscia.TIME_LIMIT_REACHED
        return true
    end

    # absolute gap reached
    absgap = tree.incumbent - tree.lb
    if absgap ≤ tree.root.abs_gap
        tree.root.Stage = Boscia.OPT_GAP_REACHED
        return true
    end

    # relative gap reached
    dual_gap = if signbit(tree.incumbent) != signbit(tree.lb)
        Inf
    elseif tree.incumbent == tree.lb
        0.0
    else
        absgap / min(abs(tree.incumbent), abs(tree.lb))
    end
    if dual_gap ≤ tree.root.rel_gap
        tree.root.Stage = Boscia.OPT_GAP_REACHED
        return true
    end

    # tree empty
    if isempty(tree.nodes)
        tree.root.Stage = Boscia.OPT_TREE_EMPTY
        return true
    end

    return false
end

"""
Own optimize function.
"""
function Bonobo.optimize!(tree::Bonobo.BnBTree{<:CustomBBNode}; callback=(args...; kwargs...)->())
    while !Bonobo.terminated(tree)
        node = Bonobo.get_next_node(tree, tree.options.traverse_strategy)
        lb, ub = Bonobo.evaluate_node!(tree, node)
        # if the problem was infeasible we simply close the node and continue
        if isnan(lb) && isnan(ub)
            Bonobo.close_node!(tree, node)
            callback(tree, node; node_infeasible=true)
            continue
        end

        Bonobo.set_node_bound!(tree.sense, node, lb, ub)
       
        # if the evaluated lower bound is worse than the best incumbent -> close and continue
        if node.lb >= tree.incumbent
            Bonobo.close_node!(tree, node)
            callback(tree, node; worse_than_incumbent=true)
            continue
        end

        tree.node_queue[node.id] = (node.lb, node.id)
        _ , prio = peek(tree.node_queue)
        @assert tree.lb <= prio[1]
        tree.lb = prio[1]

        updated = Bonobo.update_best_solution!(tree, node)
        if updated
            Bonobo.bound!(tree, node.id)
            if isapprox(tree.incumbent, tree.lb; atol=tree.options.atol, rtol=tree.options.rtol)
                tree.root.Stage = Boscia.OPT_GAP_REACHED
                callback(tree, node)
                break
            end
        end

        Bonobo.close_node!(tree, node)
        Bonobo.branch!(tree, node)
        callback(tree, node)
    end
    Bonobo.sort_solutions!(tree.solutions, tree.sense)
end
