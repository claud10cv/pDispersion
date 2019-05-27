using Dates

mutable struct Parameters
    wtype
    max_time
end

mutable struct QData
    Q
    dQ
    dToQ
end

mutable struct Data
    name
    nnodes
    D
    qdata::QData
end

mutable struct SolverStatus
    initTime::DateTime
    endTime::DateTime
    ok::Bool # if false, optimization has been aborted due to time limits
    endStatus
    endGroupsNb::Int64
end

params = Parameters(:round, 21600)
data = Data("", 0, zeros(2, 0), QData(zeros(2, 0), typemax(Int64), []))
solver_status = SolverStatus(Dates.now(), Dates.now(), true, :none, 0)
