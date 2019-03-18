using Dates

mutable struct Parameters
    wtype
    max_time
end

mutable struct Data
    name
    nnodes
    D
end

mutable struct SolverStatus
    initTime::DateTime
    endTime::DateTime
    ok::Bool # if false, optimization has been aborted due to time limits
end

params = Parameters(:round, 21600)
data = Data("", 0, [])
solver_status = SolverStatus(Dates.now(), Dates.now(), true)
