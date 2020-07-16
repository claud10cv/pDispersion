cd("../../code")
include("$(pwd())/PDispersion.jl")
p = parse(Int64, ARGS[2])
q = parse(Int64, ARGS[3])
policy = ARGS[4]
reduce = ARGS[5]
outfile = "../runs/cc2019/$(ARGS[1])_$(p)_$(q)_$(policy)_$(reduce).log"
out = open(outfile, "w")

PDispersion.read_instance_tsplib("../instances/$(ARGS[1]).tsp")
if policy == "greedy"
    nothing, Q = PDispersion.compute_lower_bound(q)
elseif policy == "optimal"
	PDispersion.set_maximum_time(86400)
    nothing, opt, nothing = PDispersion.pdispersion_decremental_clustering(q)
	if !PDispersion.optimal()
		println(out, "NA")
		close(out)
		exit(1)
	end
    Q = PDispersion.get_coordinates(opt)
elseif policy == "random"
    Q = PDispersion.get_random_coordinates(q)
end
PDispersion.set_initial_Q(Q)
lb0, nothing = PDispersion.compute_lower_bound(p)
init_nnodes = PDispersion.get_nnodes()
if reduce == "y"
    PDispersion.reduce_data_using_Q(lb0)
end
final_nnodes = PDispersion.get_nnodes()
PDispersion.set_maximum_time(86400)
init_time = time_ns()
lb, ub, opt, groups, avgSize = PDispersion.pdispersion_decremental_clustering(p)
elapsed = round(Int64, (time_ns() - init_time) * 1e-8) * 1e-1
println(out, "BestBound $lb")
println(out, "Elapsed $elapsed")
println(out, "InitialLB $lb0")
println(out, "Eliminated $(init_nnodes - final_nnodes)")
print(out, "Q")
for i in 1 : q
    print(out, " $(Q[1, i]) $(Q[2, i])")
end
println(out)
if PDispersion.optimal()
	println(out, "Optimal y")
	print(out, "BestSOL")
	optcoords = PDispersion.get_coordinates(opt)
	for i in 1 : p
	    print(out, " $(optcoords[1, i]) $(optcoords[2, i])")
	end
	println(out)
else println(out, "Optimal n")
end
close(out)
