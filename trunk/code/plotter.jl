using Plots
using BoundingSphere

function plotter(groups, opt, val; radius_factor = 0.25)
    ngroups = get_number_groups()
    Centers = zeros(ngroups, 2)
    Radius = zeros(ngroups)
    for ng in 1 : ngroups
        pts = [data.D[:, u] for u in groups[ng]]
        center, radius = boundingsphere(pts)
        Centers[ng, 1] = center[1]
        Centers[ng, 2] = center[2]
        Radius[ng] = radius
    end
    gr()
    Radius = max.(2, radius_factor * Radius)
    newCenters = vcat(Centers, data.D')
    newRadius = vcat(Radius, 2 * ones(size(data.D, 2)))
    p1 = build_left()
    p2 = build_right(p1, Centers, Radius, opt, val)
    p = plot(p1, p2, layout = 2)
end

function build_left()
    p1 = plot(t = :scatter,
                data.D[1, :]',
                data.D[2, :]',
                color=:orange,
                markersize = 2,
                legend = false,
                axis = false,
                grid = false,
                aspect_ratio = :equal,
                size = (1920, 1080))
end

function build_right(p1, Centers, Radius, opt, val)
    p2 = plot(t = :scatter,
                data.D[1, :]',
                data.D[2, :]',
                color=:orange,
                markersize = 2,
                legend = false,
                axis = false,
                grid = false,
                aspect_ratio = :equal,
                size = (1920, 1080))
    plot!(t = :scatter,
                Centers[:, 1],
                Centers[:, 2],
                markercolor = nothing,
                markeralpha = 0.0,
                markerstrokecolor = :blue,
                markerstrokealpha = 0.1,
                legend = false,
                markersize = Radius,
                axis = false,
                grid = false,
                aspect_ratio = :equal,
                size = (1920, 1080))
    for u in opt, v in opt
        if v > u
            dist = distance(u, v)
            s = Shape([(data.D[1, u], data.D[2, u]), (data.D[1, v], data.D[2, v])])
            if dist <= val
                plot!(s, linestyle = :solid, linewidth = 4, linealpha = 0.25, color = :red)
            else
                plot!(s, linestyle = :dot, linewidth = 3, linealpha = 0.25, color = :blue)
            end
        end
    end
    p2
end
