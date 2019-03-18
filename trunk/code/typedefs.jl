mutable struct Parameters
    wtype
end

mutable struct Data
    name
    nnodes
    D
end

params = Parameters(:round)
data = Data("", 0, [])
