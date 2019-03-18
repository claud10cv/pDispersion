using Clustering
using Distances
using Dates

function pdispersion_binary_search(p)
    init_solver_status()
    E = build_full_matrix()
    ub = maximum(E)
    lb, ub, opt = binarysearch(E, p, ub)
    solver_status.endTime = Dates.now()
    lb, ub, opt
end

function pdispersion_decremental_clustering(p)
    init_solver_status()
    lb = compute_lower_bound(p)
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
        else
            println("found a solution by simple inspection!")
        end
        println("solution with $(length(groups)) groups of value $val")
        if !solver_status.ok
            break
        end
        maxsize = maximum([E[g, g] for g in opt])
        if maxsize > 0
            groups, E = split_groups_and_build_matrix(groups, E; restrict_to = opt)
        else
            break
        end
    end
    solver_status.endTime = Dates.now()
    opt = [groups[u][1] for u in opt]
    ub, opt
end

function split_groups_and_build_matrix(groups, E; restrict_to = [], use_diagonal = false)
    ngroups = length(groups)
    if isempty(restrict_to)
        restrict_to = collect(1 : ngroups)
    end
    if use_diagonal
        Ediag = [E[i, i] for i in restrict_to]
        emax, imax = findmax(Ediag)
        imax = restrict_to[imax]
    else
        cands = []
        for u in restrict_to, v in restrict_to
            if v > u && E[u, u] + E[v, v] > 0
                dist = E[u, v]
                push!(cands, (u, v, E[u, v]))
            end
        end
        if isempty(cands) return copy(groups), copy(E)
        else
            sort!(lt = (x, y) -> x[3] < y[3], cands)
            u, v = cands[1][1], cands[1][2]
            if E[u, u] > E[v, v]
                imax = u
                emax = E[u, u]
            else
                imax = v
                emax = E[v, v]
            end
        end
    end

    if emax <= 0
        return copy(groups), copy(E)
    else
        group = groups[imax]
        nnodes = length(group)
        if nnodes > 2
            F = data.D[:, group]
            res = kmeans(convert.(Float64, F), 2)
            g1 = [group[v] for v in 1 : nnodes if res.assignments[v] == 1]
            g2 = [group[v] for v in 1 : nnodes if res.assignments[v] == 2]
        elseif nnodes == 2
            g1 = [group[1]]
            g2 = [group[2]]
        end
        newgroups = copy(groups)
        newgroups[imax] = g1
        push!(newgroups, g2)
        newE = vcat(E, zeros(Int64, ngroups)')
        newE = hcat(newE, zeros(Int64, ngroups + 1))
        for j in 1 : ngroups + 1
            d1j = d2j = 0
            for v in newgroups[j]
                d1j = max(d1j, maximum([distance(u, v) for u in g1]))
                d2j = max(d2j, maximum([distance(u, v) for u in g2]))
            end
            newE[imax, j] = newE[j, imax] = d1j
            newE[ngroups + 1, j] = newE[j, ngroups + 1] = d2j
        end
        return newgroups, newE
    end
end
