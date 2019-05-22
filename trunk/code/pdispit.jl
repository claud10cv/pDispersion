using Clustering
using Distances
using Dates

function pdispersion_binary_search(p; with_lb = 1)
    clean_data_points()
    init_solver_status()
#    memGB = 1e-9 * (sizeof(Int64) * data.nnodes * data.nnodes)
    if data.nnodes <= 2500
        E = build_full_matrix()
        ub = maximum(E)
        lb, ub, opt = binarysearch(E, p, ub; with_lb = with_lb)
        solver_status.endTime = Dates.now()
        solver_status.endStatus = solver_status.ok ? :optimal : :tilim
        return lb, ub, opt
    else
        solver_status.ok = false
        solver_status.endStatus = :memlim
        return 0, 0, []
    end
end

function pdispersion_decremental_clustering(p)
    clean_data_points()
    init_solver_status()
    println("computing lower bound")
    lb, nothing = compute_lower_bound(p)
    println("lower bound = $lb")
    dim, nnodes = size(data.D)
    E, groups = build_initial_groups(p)
    while maximum([E[i, i] for i in eachindex(groups)]) >= lb
        groups, E = split_groups_and_build_matrix(groups, E; use_diagonal = true)
    end
    ub = maximum(E)
    opt = []
    while true
        val = 0
        if !isempty(opt)
            opt, val = compute_easy_solution(E, p, opt)
        end
        if val < ub
            lb, ub, opt = binarysearch(E, p, ub)
        end
        if !solver_status.ok
            break
        end
        maxsize = maximum([E[g, g] for g in opt])
        if maxsize > 0
            groups, E = split_groups_and_build_matrix(groups, E; restrict_to = opt)
        else
            break
        end
        # println("val = $val")
    end
    solver_status.endTime = Dates.now()
    opt = [groups[u][1] for u in opt]
    solver_status.endStatus = solver_status.ok ? :optimal : :tilim
    solver_status.endGroupsNb = length(groups)
    println("optimal solution of $ub")
    ub, opt, groups
end
