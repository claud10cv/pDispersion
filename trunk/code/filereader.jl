function read_instance_tsplib(filename)
    start_coords = false
    open(filename) do f
        lines = readlines(f)
        for line in lines
            sline = split(line)
            if sline[1] == "EOF"
                break
            elseif start_coords
                p = parse(Int64, sline[1])
                x = parse(Float64, sline[2])
                y = parse(Float64, sline[3])
                data.D[1, p] = x
                data.D[2, p] = y
            elseif sline[1] == "NAME"
                data.name = sline[3]
            elseif sline[1] == "DIMENSION"
                data.nnodes = parse(Int64, sline[3])
                data.D = zeros(Float64, 2, data.nnodes)
            elseif sline[1] == "EDGE_WEIGHT_TYPE"
                if sline[3] == "EUC_2D"
                    params.wtype = :round
                elseif sline[3] == "CEIL_2D"
                    params.wtype = :ceil
                elseif sline[3] == "GEOM"
                    params.wtype = :geom
                end
            elseif sline[1] == "NODE_COORD_SECTION"
                start_coords = true
            end
        end
    end
    println("instance file $filename parsed successfully")
end

function read_instance_orlib(filename)
    params.wtype = :orlib
    data.name = filename
    data.nnodes = 0
    f = open(filename, "r")
    nedges = p = 0
    let
        sline = split(readline(f))
        data.nnodes = parse(Int64, sline[1])
        data.D = zeros(Int64, 2, data.nnodes)
        orlibdata.dmat = 10000000 * ones(Int64, data.nnodes, data.nnodes)
        nedges = parse(Int64, sline[2])
        p = parse(Int64, sline[3])
        for u in 1 : data.nnodes
            data.D[1, u] = data.D[2, u] = u
            orlibdata.dmat[u, u] = 0
        end
    end
    for e in 1 : nedges
        sline = split(readline(f))
        u = parse(Int64, sline[1])
        v = parse(Int64, sline[2])
        c = parse(Int64, sline[3])
        orlibdata.dmat[u, v] = orlibdata.dmat[v, u] = c
    end
    close(f)
    for k in 1 : data.nnodes, i in 1 : data.nnodes, j in 1 : data.nnodes
        orlibdata.dmat[i, j] = min(orlibdata.dmat[i, j], orlibdata.dmat[i, k] + orlibdata.dmat[k, j])
    end
    println("instance file $filename parsed successfully")
    p
end
