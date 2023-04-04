"""
    triplot(!)(points, triangles, boundary_nodes, convex_hull, constrained_edges, args...; kwargs...)
    triplot(!)(tri::Triangulation, args...; kwargs...)

Plots a triangulation. The available attributes are:

- `markersize=11`

Markersize for the points. 
- `show_ghost_edges=false`

Whether to plot interpretations of the ghost edges. If the ghost edges 
are on the outermost boundary, then a ghost edge is a line starting at 
the solid vertex and parallel with the vertex of the domain, pointing 
outwards. If the ghost edge is an interior ghost edge, then an edge 
is shown that connects with the centroid of the corresponding interior.
- `recompute_centers=true`

Only relevant if `show_ghost_edges`. This will recompute the representative centers for 
each region if needed. 
- `show_all_points=false`

Whether to plot all points. If `false`, only plots those that correspond 
to a solid triangle. 
- `point_color=:red`

Only relevant if `show_all_points`. Colour to use for the points.
- `strokecolor=:black`

Colour to use for the edges of the triangles.
- `triangle_color=(:white, 0.0)`

Colour to use for the triangles. 
- `ghost_edge_color=:blue`

Only relevant if `show_ghost_edges`. Colour to use for the ghost edges. 
- `ghost_edge_linewidth=1`

Only relevant if `show_ghost_edges`. Linewidth to use for the ghost edges.
- `strokewidth=1`

The width of the triangle edges.
- `ghost_edge_extension_factor=10.0`

Only relevant if `show_ghost_edges`. For the outermost ghost edges, we extrapolate 
a line to decide its extent away from the boundary. This is the factor 
used.
- `plot_convex_hull=true`

Whether to also plot the convex hull of the points.
- `convex_hull_color=:grey'

Only relevant if `plot_convex_hull`. The linecolor to use for the convex hull.
- `convex_hull_linestyle=:dash`

Only relevant if `plot_convex_hull`. The linestyle to use for the convex hull.
- `convex_hull_linewidth=2`

Only relevant if `plot_convex_hull`. The linewidth to use for the convex hull.
- `show_constrained_edges=true`

Only relevant if the triangulation is constrained. This colours all constrained edges.
- `constrained_edge_color=:magenta`

Only relevant if `show_constrained_edges`. The colour to use for the constrained edges.
"""
MakieCore.@recipe(Triplot, points, triangles, boundary_nodes, convex_hull, constrained_edges) do scene
    return MakieCore.Attributes(; markersize=11,
        show_ghost_edges=false,
        recompute_centers=false,
        show_all_points=false,
        point_color=:red,
        strokecolor=:black,
        triangle_color=(:white, 0.0),
        ghost_edge_color=:blue,
        ghost_edge_linewidth=1,
        strokewidth=1,
        ghost_edge_extension_factor=10.0,
        plot_convex_hull=true,
        convex_hull_color=:red,
        convex_hull_linestyle=:dash,
        convex_hull_linewidth=2,
        show_constrained_edges=true,
        constrained_edge_color=:magenta,
        constrained_edge_linewidth=2)
end
function MakieCore.convert_arguments(plot::Type{<:Triplot}, tri::Triangulation)
    return (get_points(tri), get_triangles(tri), get_boundary_nodes(tri),
        get_convex_hull(tri), get_all_constrained_edges(tri))
end

function MakieCore.plot!(p::Triplot)
    ## Extract
    points = p[:points]
    triangles = p[:triangles]
    boundary_nodes = p[:boundary_nodes]
    convex_hull = p[:convex_hull]
    constrained_edges = p[:constrained_edges]
    markersize = p[:markersize]
    show_ghost_edges = p[:show_ghost_edges]
    recompute_centers = p[:recompute_centers]
    show_all_points = p[:show_all_points]
    point_color = p[:point_color]
    strokecolor = p[:strokecolor]
    triangle_color = p[:triangle_color]
    ghost_edge_color = p[:ghost_edge_color]
    ghost_edge_linewidth = p[:ghost_edge_linewidth]
    strokewidth = p[:strokewidth]
    ghost_edge_extension_factor = p[:ghost_edge_extension_factor]
    plot_convex_hull = p[:plot_convex_hull]
    convex_hull_color = p[:convex_hull_color]
    convex_hull_linestyle = p[:convex_hull_linestyle]
    convex_hull_linewidth = p[:convex_hull_linewidth]
    show_constrained_edges = p[:show_constrained_edges]
    constrained_edge_color = p[:constrained_edge_color]
    constrained_edge_linewidth = p[:constrained_edge_linewidth]

    ## Define all the necessary observables
    points_2f = MakieCore.Observable(NTuple{2,Float64}[])
    triangle_mat = MakieCore.Observable(NTuple{3,Int64}[])
    triangle_mat_2 = MakieCore.Observable{Matrix{Int64}}()
    ghost_edges = MakieCore.Observable(NTuple{2,Float64}[])
    convex_hull_points = MakieCore.Observable(NTuple{2,Float64}[])
    constrained_edge_points = MakieCore.Observable(NTuple{2,Float64}[])

    ## Define the function for updating the observables 
    function update_plot(points, triangles, boundary_nodes, convex_hull, constrained_edges)
        boundary_map = construct_boundary_map(boundary_nodes)
        if !has_multiple_segments(boundary_nodes) && num_boundary_edges(boundary_nodes) == 0
            recompute_centers[] &&
                compute_representative_points!(points, get_indices(convex_hull))
        else
            recompute_centers[] && compute_representative_points!(points, boundary_nodes)
        end

        ## Clear out the previous observables
        empty!(points_2f[])
        empty!(triangle_mat[])
        empty!(ghost_edges[])
        empty!(convex_hull_points[])
        empty!(constrained_edge_points[])

        ## Fill out the points 
        for pt in each_point(points)
            x, y = getxy(pt)
            push!(points_2f[], (x, y))
        end

        ## Fill out the triangles 
        for T in each_triangle(triangles)
            if !is_ghost_triangle(T)
                i, j, k = indices(T)
                push!(triangle_mat[], (i, j, k))
            end
        end
        triangle_mat_2[] = Matrix(reinterpret(reshape, Int64, triangle_mat[])')

        ## Now get the ghost edges if needed 
        if show_ghost_edges[]
            for index in values(boundary_map)
                curve_index = get_curve_index(index)
                representative_coordinates = get_representative_point_coordinates(curve_index,
                    number_type(points))
                cx, cy = getxy(representative_coordinates)
                bn = !isempty(index) ? get_boundary_nodes(boundary_nodes, index) :
                     get_indices(convex_hull) # If index is empty, there are no constrained boundaries, meaning we want to get the convex hull
                ## TODO: This needs to get fixed when we use delete_point!, since that will give a convex hull that could end up not 
                ## being right for the ghost edges in the case where we delete a boundary point.
                n_edge = num_boundary_edges(bn)
                for i in 1:n_edge
                    u = get_boundary_nodes(bn, i)
                    q = get_point(points, u)
                    qx, qy = getxy(q)
                    if !is_interior_curve(curve_index)
                        end_edge_x = cx * (1 - ghost_edge_extension_factor[]) +
                                     ghost_edge_extension_factor[] * qx
                        end_edge_y = cy * (1 - ghost_edge_extension_factor[]) +
                                     ghost_edge_extension_factor[] * qy
                    else
                        end_edge_x = cx
                        end_edge_y = cy
                    end
                    push!(ghost_edges[], (qx, qy))
                    push!(ghost_edges[], (end_edge_x, end_edge_y))
                end
            end
        end

        ## Get the convex hull edges if needed 
        if plot_convex_hull[]
            idx = get_indices(convex_hull)
            for i in idx
                pt = get_point(points, i)
                x, y = getxy(pt)
                push!(convex_hull_points[], (x, y))
            end
        end

        ## Get the constrained edges if needed 
        if show_constrained_edges[]
            for e in each_edge(constrained_edges)
                p1 = get_point(points, initial(e))
                p2 = get_point(points, terminal(e))
                push!(constrained_edge_points[], getxy(p1), getxy(p2))
            end
        end
    end

    ## Connect the plot so that it updates whenever we change a value 
    MakieCore.Observables.onany(update_plot, points, triangles, boundary_nodes, convex_hull, constrained_edges)

    ## Call it once to prepopulate with current values 
    update_plot(points[], triangles[], boundary_nodes[], convex_hull[], constrained_edges[])

    ## Now plot 
    poly!(p, points_2f, triangle_mat_2;
        strokewidth=strokewidth[],
        strokecolor=strokecolor[],
        color=triangle_color[])
    if show_all_points[]
        scatter!(p, points_2f; markersize=markersize[], color=point_color[])
    end
    if show_ghost_edges[]
        linesegments!(p, ghost_edges; color=ghost_edge_color[],
            linewidth=ghost_edge_linewidth[])
    end
    if plot_convex_hull[]
        lines!(p, convex_hull_points; color=convex_hull_color[],
            linewidth=convex_hull_linewidth[], linestyle=convex_hull_linestyle[])
    end
    if show_constrained_edges[]
        linesegments!(p, constrained_edge_points; color=constrained_edge_color[],
            linewidth=constrained_edge_linewidth[])
    end
    return p
end
