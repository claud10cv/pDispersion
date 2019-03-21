using Gurobi
using JuMP
using Random
using Distances
using MathProgBase

function pdisp(D, p, lb, ub)
    nnodes = size(D, 1)
    dmax = maximum(D)
    Dlist = [[] for i in 1 : dmax + 1]
    unfeas = []
    adj_unfeas = trues(nnodes, nnodes)
    for i in 1 : nnodes - 1, j in i + 1 : nnodes
        if D[i, j] >= lb && D[i, j] <= ub
            push!(Dlist[D[i, j] + 1], (i, j))
        elseif D[i, j] < lb
            push!(unfeas, (i, j))
            adj_unfeas[i, j] = adj_unfeas[j, i] = false
        end
    end
    dmap = [d for d in 0 : dmax if !isempty(Dlist[d + 1])]
    if isempty(dmap) return [], 0
    end
    ndists = length(dmap)
    adj_k = [trues(nnodes, nnodes) for k in 1 : ndists]
    for k in 2 : ndists
        adj_k[k] = copy(adj_k[k - 1])
        for (i, j) in Dlist[dmap[k - 1]]
            adj_k[k][i, j] = adj_k[k][j, i] = false
        end
    end
    nunfeas = length(unfeas)
    maxtime = max(1, params.max_time - elapsed())
    m = Model(solver = GurobiSolver(OutputFlag = 0,
                                    Threads = 1,
                                    MIPFocus = 1,
                                    TimeLimit = maxtime + 1))
    @variable(m, x[1 : nnodes], Bin)
    @variable(m, z[1 : ndists], Bin)
    @objective(m, Max, dmap[1] * z[1] + sum((dmap[k] - dmap[k - 1]) * z[k] for k in 2 : ndists))
    @constraint(m, sum(x) == p)
    @constraint(m, monot[k in 2 : ndists], z[k] - z[k - 1] <= 0)
    @constraint(m, [k in 2 : ndists, (i, j) in Dlist[dmap[k - 1] + 1]], x[i] + x[j] + z[k] <= 2)
    @constraint(m, [i in 1 : nunfeas], x[unfeas[i][1]] + x[unfeas[i][2]] <= 1)

    function lazycb(cb)
        cbnodes = MathProgBase.cbgetexplorednodes(cb)
        if cbnodes >= 100 return
        end
        xvals = getvalue(x)
        zvals = getvalue(z)
        if sum(abs.(xvals - round.(Int64, xvals))) > 1e-5 && sum(abs.(zvals[2 : ndists] - round.(Int64, zvals[2 : ndists]))) > 1e-5
            return
        end
        let
            sol = maximum_weighted_clique(nnodes, adj_unfeas, xvals)
            if sum(xvals[u] for u in sol) > 1.1
                # println("adding clique inequality of type I")
                @lazyconstraint(cb, sum(x[u] for u in sol) <= 1)
                return
            end
        end
        for k in 2 : ndists
            # sum(x_i i in S) + (|S| - 1)z[k] <= |S|
            weights = [xvals[u] + zvals[k] - 1 for u in 1 : nnodes]
            sol = maximum_weighted_clique(nnodes, adj_k[k], weights)
            nsol = length(sol)
            if nsol > 0 && sum(xvals[u] for u in sol) + (nsol - 1) * zvals[k] > nsol + 1e-1
                # println("adding clique inequality of type II")
                @lazyconstraint(cb, sum(x[u] for u in sol) + (nsol - 1) * z[k] <= nsol)
            end
        end
    end

    feas = []
    function infocb(cb)
        xvals = round.(Int64, getvalue(x))
        sol = [i for i in 1 : nnodes if xvals[i] > 0]
        push!(feas, sol)
        throw(CallbackAbort())
    end

    addlazycallback(m, lazycb; fractional = true)
    addinfocallback(m, infocb, when = :MIPSol)

    status = solve(m)
    if !isempty(feas)
        sol = feas[1]
        optval = minimum([D[i, j] for i in sol, j in sol if j > i])
        return sol, optval
    else return [], 0
    end
end

function binarysearch(D, p, ub; with_lb = 1)
    newlb = ub
    newub = ub
    mult = 2
    sol = []
    while true
        newsol, newoptval = pdisp(D, p, newlb, newub)
        if !isempty(newsol) #problem feasible
            newlb = newoptval
            sol = newsol
            break
        elseif elapsed() < params.max_time
            newub = newlb - 1
            newlb = max(with_lb, newlb - mult)
            mult *= 2
        else
            solver_status.ok = false
            return 0, newub, sol
        end
    end
    while newlb < newub
        r = ceil(Int64, (newlb + newub) / 2)
        newsol, newoptval = pdisp(D, p, r, newub)
        if !isempty(newsol) #problem feasible
            newlb = newoptval
            sol = newsol
        elseif elapsed() < params.max_time
            newub = r - 1
        else
            solver_status.ok = false
            break
        end
    end
    return newlb, newub, sol
end
