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
    p1 = scatter(data.D[1, :]',
                data.D[2, :]',
                color=:orange,
                markersize = 2,
                legend = false,
                axis = false,
                grid = false,
                aspect_ratio = :equal,
                size = (1920, 1080))
    p2 = scatter(Centers[:, 1],
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
                plot!(s, linestyle = :solid, linewidth = 2, linealpha = 0.25, color = :red)
            else
                plot!(s, linestyle = :dot, linewidth = 0.5, linealpha = 0.1, color = :cyan)
            end
        end
    end
    # layout = @layout [ a{0.5w, 0.5h} b{0.5w, 0.5h} ]
    p = plot(p1, p2, layout = 2)
end
