using Distances
using Gurobi
using JuMP
using Dates

function distance(u, v)
    if params.wtype == :ceil
        ceil(Int64, euclidean(data.D[:, u], data.D[:, v]))
    else
        round(Int64, euclidean(data.D[:, u], data.D[:, v]))
    end
end

function build_full_matrix()
    E = zeros(Int64, data.nnodes, data.nnodes)
    for i in 1 : data.nnodes - 1, j in i + 1 : data.nnodes
        E[i, j] = E[j, i] = distance(i, j)
    end
    E
end

function clean_data_points()
    init_data_points = data.nnodes
    E = [(data.D[1, u], data.D[2, u]) for u in 1 : data.nnodes]
    # println(E)
    sort!(E)
    unique!(E)
    data.nnodes = length(E)
    data.D = zeros(Float64, 2, data.nnodes)
    for u in 1 : data.nnodes
        data.D[1, u] = E[u][1]
        data.D[2, u] = E[u][2]
    end
    end_data_points = data.nnodes
    println("removed $(init_data_points - end_data_points) points from the data")
end

function maximum_weighted_clique_exact(nnodes, adj, weights)
    m = Model(solver = GurobiSolver(OutputFlag = 0, Threads = 1))
    @variable(m, x[1 : nnodes], Bin)
    for i in 1 : nnodes, j in i + 1 : nnodes
        if adj[i, j]
            @constraint(m, x[i] + x[j] <= 1)
        end
    end
    @objective(m, Max, JuMP.dot(weights, x))
    status = solve(m)
    xvals = round.(Int64, getvalue(x))
    opt = [i for i in 1 : nnodes if xvals[i] > 0]
end

function maximum_weighted_clique_heuristic(nnodes, adj, weights)
    cands = collect(1 : nnodes)
    incands = trues(nnodes)
    sol = []
    while length(cands) > 1
        deg = zeros(Int64, nnodes)
        for u in cands
            deg[u] = sum(!adj[u, v] for v in 1 : nnodes if incands[v] && v != u)
        end
        sort!(lt = (u, v) -> deg[u] * weights[u] > deg[v] * weights[v], cands)
        u = popfirst!(cands)
        push!(sol, u)
        torem = [u]
        for v in cands
            if v != u && adj[u, v]
                push!(torem, v)
            end
        end
        for v in torem
            incands[v] = false
        end
        filter!(x -> incands[x], cands)
    end
    if !isempty(cands) push!(sol, cands[1])
    end
    sort!(sol)
    sol
end

function maximum_weighted_clique(nnodes, adj, weights)
    maximum_weighted_clique_heuristic(nnodes, adj, weights)
end

function build_initial_groups(p)
    E = zeros(Int64, p, p)
    res = kmeans(convert.(Float64, data.D), p)
    groups = [[] for i in 1 : p]
    for u in 1 : data.nnodes
        k = res.assignments[u]
        push!(groups[k], u)
    end
    for k in 1 : p
        for u in groups[k]
            du = maximum([distance(u, v) for v in groups[k] if v >= u])
            E[k, k] = max(E[k, k], du)
        end
    end
    for k in 1 : p - 1, l in k + 1 : p
        for u in groups[k]
            du = maximum([distance(u, v) for v in groups[l]])
            E[k, l] = max(E[k, l], du)
        end
        E[l, k] = E[k, l]
    end
    E, groups
end

function compute_lower_bound2(p)
    function d0(u)
        if params.wtype == :ceil
            ceil(Int64, euclidean(data.D[:, u], [0, 0]))
        else
            round(Int64, euclidean(data.D[:, u], [0, 0]))
        end
    end
    opt_coords = zeros(2, 0)
    let
        dists = [d0(u) for u in 1 : data.nnodes]
        val, u = findmax(dists)
        opt_coords = hcat(opt_coords, data.D[:, u])
    end
    function d(u, v)
        if params.wtype == :ceil
            ceil(Int64, euclidean(data.D[:, u], opt_coords[:, v]))
        else
            round(Int64, euclidean(data.D[:, u], opt_coords[:, v]))
        end
    end
    function mind(u)
        minimum([d(u, v) for v in 1 : size(opt_coords, 2)])
    end
    lb = typemax(Int64)
    while size(opt_coords, 2) < p
        dists = [mind(u) for u in 1 : data.nnodes]
        val, u = findmax(dists)
        opt_coords = hcat(opt_coords, data.D[:, u])
        lb = min(lb, val)
    end
    lb
end
function compute_lower_bound(p)
    sorted = collect(1 : data.nnodes)
    lb = 0
    nonsuccess = 0
    while nonsuccess < 10
        Random.shuffle!(sorted)
        initcenters = data.D[:, sorted[1 : p]]
        res = kmeans!(convert.(Float64, data.D), convert.(Float64, initcenters))
        groups = [[u for u in 1 : data.nnodes if res.assignments[u] == k] for k in 1 : p]
        centers = []
        for k in 1 : p
            minimax = 1e+20
            argminmax = 0
            for u in groups[k]
                roundfn = params.wtype == :ceil ? ceil : round
                dmax = roundfn(Int64, euclidean(data.D[:, u], initcenters[:, k]))#maximum([ceil(Int64, euclidean(D[:, u], D[:, v])) for v in groups[k]])
                if dmax < minimax
                    minimax = dmax
                    argminmax = u
                end
            end
            push!(centers, argminmax)
        end
        newlb = minimum([distance(u, v) for u in centers, v in centers if v > u])
        if newlb > lb
            lb = newlb
            nonsuccess = 0
        else nonsuccess += 1
        end
    end
    lb
end

function compute_easy_solution(E, p, oldopt)
    ngroups = size(E, 1)
    opt = copy(oldopt)
    push!(opt, ngroups)
    if length(opt) != p + 1
        return [], 0
    else
        val = 0
        newopt = []
        for u in opt
            sol = [v for v in opt if v != u]
            dmin = minimum(E[u, v] for u in sol, v in sol if v > u)
            if dmin > val
                val = dmin
                newopt = sol
            end
        end
        return newopt, val
    end
end

function split_groups_and_build_matrix(groups, E; restrict_to = [], use_diagonal = false)
    ngroups = length(groups)
    if isempty(restrict_to)
        restrict_to = collect(1 : ngroups)
    end
    if maximum([length(g) for g in groups[restrict_to]]) <= 1
        return copy(groups), copy(E)
    end
    if use_diagonal
        Ediag = [(i, E[i, i], length(groups[i])) for i in restrict_to]
        sort!(lt = (u, v) -> u[2] > v[2] || (u[2] == v[2] && u[3] > v[3]), Ediag)
        imax = Ediag[1][1]
    else
        cands = []
        for u in restrict_to, v in restrict_to
            if v > u && length(groups[u]) + length(groups[v]) > 2
                push!(cands, (u, v, E[u, v]))
            end
        end
        sort!(lt = (x, y) -> x[3] < y[3], cands)
        u, v = cands[1][1], cands[1][2]
        if E[u, u] > E[v, v] || (E[u, u] == E[v, v] && length(groups[u]) > length(groups[v]))
            imax = u
            emax = E[u, u]
        else
            imax = v
            emax = E[v, v]
        end
    end

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

function set_maximum_time(s)
    params.max_time = s
end

function elapsed()
    b = Dates.now()
    ms = Dates.value(b - solver_status.initTime)
    s = ms / 1000
end

function total_elapsed()
    ms = Dates.value(solver_status.endTime - solver_status.initTime)
    s = ms / 1000
end

function init_solver_status()
    solver_status.initTime = Dates.now()
    solver_status.endTime = Dates.now()
    solver_status.ok = true
    solver_status.endStatus = :none
    solver_status.endGroupsNb = 0
end

function optimal()
    return solver_status.endStatus == :optimal
end

function get_status()
    return solver_status.endStatus
end

function get_number_groups()
	return solver_status.endGroupsNb
end
