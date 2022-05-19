# pDispersion
Exact solver for the p-dispersion problem

## Basic usage
```julia
PDispersion.read_instance_tsplib(filename) # read instance
PDispersion.set_maximum_time(3600) # set time limit of one hour
Q = PDispersion.get_random_coordinates(2) # select two points at random
PDispersion.set_initial_Q(Q) # fix these two points
lb, ub, opt, groups, avgSize = PDispersion.pdispersion_decremental_clustering(5) # solve to select additional five points 
```
