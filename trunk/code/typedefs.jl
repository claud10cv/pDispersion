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

mutable struct QPData
    Q
    inQ
end

mutable struct Data
    name
    nnodes
    D
    qdata::QData
    qpdata::QPData
end

mutable struct OrlibData
    dmat
end

mutable struct SolverStatus
    initTime::DateTime
    endTime::DateTime
    ok::Bool # if false, optimization has been aborted due to time limits
    endStatus
    endGroupsNb::Int64
end

params = Parameters(:round, 21600)
data = Data("", 0, zeros(2, 0), QData(zeros(2, 0), typemax(Int64), []), QPData(zeros(2, 0), []))
solver_status = SolverStatus(Dates.now(), Dates.now(), true, :none, 0)
orlibdata = OrlibData(zeros(0, 0))
