function readinstance_tsplib(filename)
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
