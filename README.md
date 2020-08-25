# pDispersion
Exact solver for the p-dispersion problem

## Basic usage
```julia
PDispersion.read_instance_tsplib(filename)
PDispersion.set_maximum_time(3600)
lb, ub, opt, groups, avgSize = PDispersion.pdispersion_decremental_clustering(5)
```
