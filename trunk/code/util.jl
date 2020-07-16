using Distances
using Gurobi
using JuMP
using Dates

function rounding(T, w)
    if params.wtype == :ceil
        return ceil(T, w)
    else return round(T, w)
    end
end

function distance(u::T where T<:Integer, v::T where T<:Integer)
    q = size(data.qdata.Q, 2)
    if params.wtype == :orlib
        m1 = orlibdata.dmat[u, v]
    else
        m1 = rounding(Int64, new_euclidean(data.D[:, u], data.D[:, v]))
    end
    if q > 0
        m2 = data.qdata.dQ
        m3 = data.qdata.dToQ[u]#minimum([rounding(Int64, euclidean(data.D[:, u], data.Q[:, w])) for w in 1 : q])
        m4 = data.qdata.dToQ[v]#minimum([rounding(Int64, euclidean(data.D[:, u], data.Q[:, w])) for w in 1 : q])
        return minimum([m1, m2, m3, m4])
    else
        return m1
    end
end

function distance(u::T where T<:Integer)
    q = size(data.qdata.Q, 2)
    if q > 0
        m1 = data.qdata.dQ
        m2 = data.qdata.dToQ[u]
        return min(m1, m2)
    else return typemax(Int64)
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
    println("Cleaning repeated data points")
    init_data_points = data.nnodes
    if size(data.qdata.Q, 2) > 0
        E = [(data.D[1, u], data.D[2, u], data.qdata.dToQ[u]) for u in 1 : data.nnodes]
    else
        E = [(data.D[1, u], data.D[2, u]) for u in 1 : data.nnodes]
    end
    sort!(E)
    unique!(E)
    data.nnodes = length(E)
    data.D = zeros(Float64, 2, data.nnodes)
    if size(data.qdata.Q, 2) > 0
        data.qdata.dToQ = zeros(Int64, data.nnodes)
    end
    for u in 1 : data.nnodes
        data.D[1, u] = E[u][1]
        data.D[2, u] = E[u][2]
        if size(data.qdata.Q, 2) > 0
            data.qdata.dToQ[u] = E[u][3]
        end
    end
    if size(data.qpdata.Q, 2) > 0
        q = size(data.qpdata.Q, 2)
        data.qpdata.inQ = falses(data.nnodes)
        for u in 1 : data.nnodes
            for i in 1 : q
                if data.qpdata.Q[:, i] == data.D[:, u]
                    data.qpdata.inQ[u] = true
                    break
                end
            end
        end
    end
    end_data_points = data.nnodes
    println("removed $(init_data_points - end_data_points) points from the data containing $init_data_points points")
end

function compute_heuristic_lower_bound(g)
    println("computing heuristic lb")
    opt, lb = [], 0
    unsuccessful = 0
    while unsuccessful < 10
        sol = [rand(q) for q in g]
        val = minimum([distance(u, v) for u in sol, v in sol if u != v])
        if val > lb
            opt, lb = sol, val
            unsuccessful = 0
        else unsuccessful += 1
        end
    end
    opt, lb
end

function reduce_data_using_Q(lb)
    println("reducing data using Q")
    init_data_points = data.nnodes
    if size(data.qdata.Q, 2) > 0
        E = [(data.D[1, u], data.D[2, u], data.qdata.dToQ[u]) for u in 1 : data.nnodes if distance(u) >= lb]
    else
        E = [(data.D[1, u], data.D[2, u]) for u in 1 : data.nnodes if distance(u) >= lb]
    end
    sort!(E)
    unique!(E)
    data.nnodes = length(E)
    data.D = zeros(Float64, 2, data.nnodes)
    if size(data.qdata.Q, 2) > 0
        data.qdata.dToQ = zeros(Int64, data.nnodes)
    end
    for u in 1 : data.nnodes
        data.D[1, u] = E[u][1]
        data.D[2, u] = E[u][2]
        if size(data.qdata.Q, 2) > 0
            data.qdata.dToQ[u] = E[u][3]
        end
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
    nq = size(data.qpdata.Q, 2)
    E = zeros(Int64, p, p)
    res = kmeans(convert.(Float64, data.D), p - nq)
    groups = [[] for i in 1 : p]
    currq = 0
    for u in 1 : data.nnodes
        if nq > 0 && data.qpdata.inQ[u]
            currq += 1
            push!(groups[currq], u)
        else
            k = res.assignments[u] + nq
            push!(groups[k], u)
        end
    end
    println("nq = $nq")
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
    fixed = size(data.qpdata.Q, 2) > 0
    q = fixed ? size(data.qpdata.Q, 2) : size(data.qdata.Q, 2)
    println("computing lower bound with p = $p and q = $q")
    opt_coords = fixed ? copy(data.qpdata.Q) : copy(data.qdata.Q)
    nopt = size(opt_coords, 2)
    if nopt <= 0
        dists = params.wtype == :orlib ? [maximum(orlibdata.dmat[u, :]) for u in 1 : data.nnodes] : [rounding(Int64, new_euclidean(data.D[:, u], [0, 0])) for u in 1 : data.nnodes]
        val, u = findmax(dists)
        opt_coords = hcat(opt_coords, data.D[:, u])
        nopt += 1
    end
    function mind(u)
        minimum([rounding(Int64, new_euclidean(data.D[:, u], opt_coords[:, v])) for v in 1 : nopt])
    end
    if q > 0
        lb = data.qdata.dQ
    else lb = typemax(Int64)
    end
    maxsize = fixed ? p : p + q
    while nopt < maxsize
        dists = [mind(u) for u in 1 : data.nnodes]
        # println("dists = $dists")
        val, u = findmax(dists)
        # println("val = $val")
        opt_coords = hcat(opt_coords, data.D[:, u])
        nopt += 1
        lb = min(lb, val)
    end
    println("lower bound = $lb")
    lb, opt_coords
end

function compute_easy_solution(E, p, oldopt)
    nqp = size(data.qpdata.Q, 2)
    ngroups = size(E, 1)
    opt = copy(oldopt)
    push!(opt, ngroups)
    if length(opt) != p + 1
        return [], 0
    else
        val = 0
        newopt = []
        for u in opt
            if u > nqp
                sol = [v for v in opt if v != u]
                dmin = minimum(E[u, v] for u in sol, v in sol if v > u)
                if dmin > val
                    val = dmin
                    newopt = sol
                end
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

function get_nnodes()
    data.nnodes
end

function set_data(name, nnodes, D)
    data = Data(name, nnodes, D)
end

function set_params(wtype, max_time)
    params = Parameters(wtype, max_time)
end

function get_coordinates(opt)
    data.D[:, opt]
end

function reset()
    params = Parameters(:round, 21600)
    data = Data("", 0, [])
    solver_status = SolverStatus(Dates.now(), Dates.now(), true, :none, 0)
end

function set_initial_Q(coords)
    data.qdata.Q = coords
    q = size(coords, 2)
    dQ = -1
    for i in 1 : q - 1, j in i + 1 : q
        euc = new_euclidean(coords[:, i], coords[:, j])
        if dQ < 0
            if params.wtype == :ceil
                dQ = ceil(Int64, euc)
            else dQ = round(Int64, euc)
            end
        elseif params.wtype == :ceil
            dQ = min(dQ, ceil(Int64, euc))
        else dQ = min(dQ, round(Int64, euc))
        end
    end
    data.qdata.dQ = dQ
    data.qdata.dToQ = zeros(Int64, data.nnodes)
    for u in 1 : data.nnodes
        data.qdata.dToQ[u] = minimum([rounding(Int64, new_euclidean(data.D[:, u], coords[:, v])) for v in 1 : q])
    end
end

function set_initial_Qp(coords)
    q = size(coords, 2)
    data.qpdata.Q = coords
    data.qpdata.inQ = falses(data.nnodes)
    nfix = 0
    for i in 1 : q
        for u in 1 : data.nnodes
            if coords[:, i] == data.D[:, u]
                data.qpdata.inQ[u] = true
                nfix += 1
                break
            end
        end
    end
    println("successfully fixed $nfix nodes")
end

function get_random_coordinates(q)
    data.D[:, rand(1 : data.nnodes, q)]
end

function new_euclidean(U, V)
    if params.wtype == :orlib
        return orlibdata.dmat[round(Int64, U[1]), round(Int64, V[1])]
    else
        return euclidean(U, V)
    end
end
