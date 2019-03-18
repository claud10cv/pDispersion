using Gurobi
using JuMP
using Random
using Distances
using MathProgBase

function pdisp(D, p, lb, ub)
    # println(D)
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
    ndists = length(dmap)
    adj_k = [trues(nnodes, nnodes) for k in 1 : ndists]
    for k in 2 : ndists
        adj_k[k] = copy(adj_k[k - 1])
        for (i, j) in Dlist[dmap[k - 1]]
            adj_k[k][i, j] = adj_k[k][j, i] = false
        end
    end
    nunfeas = length(unfeas)
    m = Model(solver = GurobiSolver(OutputFlag = 0))
    @variable(m, x[1 : nnodes], Bin)
    @variable(m, z[2 : ndists], Bin)
    @objective(m, Max, sum((dmap[k] - dmap[k - 1]) * z[k] for k in 2 : ndists))
    @constraint(m, sum(x) == p)
    @constraint(m, monot[k in 3 : ndists], z[k] - z[k - 1] <= 0)
    @constraint(m, [k in 2 : ndists, (i, j) in Dlist[dmap[k - 1] + 1]], x[i] + x[j] + z[k] <= 2)
    @constraint(m, [i in 1 : nunfeas], x[unfeas[i][1]] + x[unfeas[i][2]] <= 1)

    function lazycb(cb)
        cbnodes = MathProgBase.cbgetexplorednodes(cb)
        if cbnodes >= 1 return
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
    addlazycallback(m, lazycb; fractional = true)
    status = solve(m)
    if status == :Infeasible
        return [], 0
    else
        xvals = round.(Int64, getvalue(x))
        optval = round(Int64, dmap[1] + getobjectivevalue(m))
        sol = [i for i in 1 : nnodes if xvals[i] > 0]
        return sol, optval
    end
end

function pdisp_binarysearch(D, p, ub)
    newlb = ub
    newub = ub
    mult = 2
    sol, optval = [], -1
    while true
        newsol, newoptval = pdisp(D, p, newlb, newub)
        if !isempty(newsol) #problem feasible
            newlb = newoptval
            sol = newsol
            optval = newoptval
            break
        else
            newub = newlb - 1
            newlb = max(1, newlb - mult)
            mult *= 2
        end
    end
    println("finished powsearch with lb = $newlb and ub = $newub")
    while newlb < newub
        r = ceil(Int64, ( newlb + newub ) / 2)
        newsol, newval = pdisp(D, p, r, newub)
        if !isempty(newsol)
            newlb = optval = newval
            sol = newsol
        else
            newub = r - 1
        end
    end
    println("finished binsearch with lb = $newlb and ub = $newub")
    sol, optval
end
