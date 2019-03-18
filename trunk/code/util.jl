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

function maximum_weighted_clique_exact(nnodes, adj, weights)
    m = Model(solver = GurobiSolver(OutputFlag = 0))
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
end

function optimal()
    return solver_status.endStatus == :optimal
end

function getStatus()
    return solver_status.endStatus
end
