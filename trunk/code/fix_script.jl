cd("../../code")
include("$(pwd())/PDispersion.jl")
include("$(pwd())/extractQ.jl")
p = parse(Int64, ARGS[2])
q = parse(Int64, ARGS[3])
outfile = "../runs/cc2019/$(ARGS[1])_$(p)_$(q)_fix.log"
out = open(outfile, "w")
fileToExtract = "../runs/cc2019/greedy-optimal/$(ARGS[1])_$(p)_$(q)_greedy.log"
Q = extractQ(fileToExtract)
PDispersion.read_instance_tsplib("../instances/$(ARGS[1]).tsp")
PDispersion.set_initial_Qp(Q)
PDispersion.set_maximum_time(86400)
init_time = time_ns()
lb, ub, opt, groups, avgs = PDispersion.pdispersion_decremental_clustering(p + q)
elapsed = round(Int64, (time_ns() - init_time) * 1e-8) * 1e-1
println(out, "BestBound $lb")
println(out, "Elapsed $elapsed")
println(out, "InitialLB 0")
println(out, "Eliminated 0")
print(out, "Q")
for i in 1 : q
    print(out, " $(Q[1, i]) $(Q[2, i])")
end
println(out)
if PDispersion.optimal()
	println(out, "Optimal y")
	print(out, "BestSOL")
	optcoords = PDispersion.get_coordinates(opt)
	for i in 1 : p + q
	    print(out, " $(optcoords[1, i]) $(optcoords[2, i])")
	end
	println(out)
else println(out, "Optimal n")
end
close(out)
