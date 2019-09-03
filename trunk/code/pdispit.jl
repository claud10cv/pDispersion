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
    lb, bks = compute_lower_bound(p)
	println("bks = $bks")
    println("lower bound = $lb")
    dim, nnodes = size(data.D)
    E, groups = build_initial_groups(p)
    while maximum([E[i, i] for i in eachindex(groups)]) >= lb
        groups, E = split_groups_and_build_matrix(groups, E; use_diagonal = true)
    end
    ub = maximum(E)
	if !isempty(data.qdata.Q)
		ub = data.qdata.dQ
	end
    opt = []
   if lb >= ub
	nothing, nbks = size(bks)
	for i in 1 : nbks
		for u in 1 : data.nnodes
			if data.D[:, u] == bks[:, i]
				push!(opt, u)
				break
			end
		end
	end
	solver_status.endTime = Dates.now()
	solver_status.endStatus = :optimal
	solver_status.endGroupsNb = data.nnodes
	return lb, opt, [[u] for u in 1 : data.nnodes]
   end
   while lb < ub
        val = 0
        if !isempty(opt)
            opt, val = compute_easy_solution(E, p, opt)
        end
        if val < ub
            rlb, ub, opt = binarysearch(E, p, ub)
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
    sol, lb = compute_heuristic_lower_bound([groups[u] for u in opt])
    # opt = [groups[u][1] for u in opt]
    avgSize = round(Int64, sum(length(groups[u]) for u in opt) * 100 / p) / 100
    solver_status.endStatus = solver_status.ok ? :optimal : :tilim
    solver_status.endGroupsNb = length(groups)
    println("upper bound of $ub")
    lb, ub, sol, groups, avgSize
end
