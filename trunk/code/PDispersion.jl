module PDispersion
    include("typedefs.jl")
    include("util.jl")
    include("pdisp.jl")
    include("pdispit.jl")
    include("filereader.jl")
    include("plotter.jl")
    export pdispersion_decremental_clustering
    export pdispersion_binary_search
    export read_instance_tsplib
    export total_elapsed
    export optimal
end
