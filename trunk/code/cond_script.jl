#cd("../../code")
include("$(pwd())/PDispersion.jl")
filename = ARGS[1]
p = parse(Int64, ARGS[2])
q = parse(Int64, ARGS[3])
policy = ARGS[4]
tok = split(filename, ['.', '/']; keepempty = false)
ntok = length(tok)
outfile = "$(tok[ntok - 1])_$(p)_$(q)_$policy.log"
out = open(outfile, "w")

PDispersion.read_instance_tsplib(filename)
if policy == "greedy"
    nothing, Q = PDispersion.compute_lower_bound(q)
elseif policy == "optimal"
    nothing, opt, nothing = PDispersion.pdispersion_decremental_clustering(q)
    Q = PDispersion.get_coordinates(opt)
elseif policy == "random"
    Q = PDispersion.get_random_coordinates(q)
end
PDispersion.set_initial_Q(Q)
lb0, nothing = PDispersion.compute_lower_bound(p)
init_nnodes = PDispersion.get_nnodes()
PDispersion.reduce_data_using_Q(lb0)
final_nnodes = PDispersion.get_nnodes()
lb, opt, groups = PDispersion.pdispersion_decremental_clustering(p)
println(out, "Best $lb")
println(out, "Initial $lb0")
println(out, "Eliminated $(init_nnodes - final_nnodes)")
print(out, "Q")
for i in 1 : q
    print(out, " $(Q[1, i]) $(Q[2, i])")
end
println(out)
print(out, "Opt")
optcoords = PDispersion.get_coordinates(opt)
for i in 1 : p
    print(out, " $(optcoords[1, i]) $(optcoords[2, i])")
end
println(out)
