module DelaunayTriangulationMakieCoreExt

using DelaunayTriangulation
import MakieCore: MakieCore, @recipe

"""
    triplot(!)(points, triangles, boundary_nodes, convex_hull, constrained_edges, representative_point_list; kwargs...)
    triplot(!)(tri::Triangulation; kwargs...)

Plots a triangulation. 

# Attributes 
- `markersize=11`

Size of the points.
- `show_ghost_edges=false`

Whether to show ghost edges.
- `recompute_centers=false`

Whether to recompute `tri.representative_point_list`.
- `show_all_points=false`

Whether to show all points, regardless of whether they appear in the triangulation.
- `point_color=:red`

Colour of the points.
- `strokecolor=:black`

Colours of the triangle edges.
- `triangle_color=(:white, 0.0)`

Colours of the triangles.
- `ghost_edge_color=:blue`

Colours of the ghost edges.
- `ghost_edge_linewidth=1`

Width of the ghost edges.
- `strokewidth=1`

Width of the triangle edges.
- `ghost_edge_extension_factor=10.0`

Factor that scales the length of the ghost edges.
- `show_convex_hull=true`

Whether to show the convex hull.
- `convex_hull_color=:red`

Colour of the convex hull.
- `convex_hull_linestyle=:dash`

Linestyle for the convex hull edges.
- `convex_hull_linewidth=2`

Width of the convex hull.
- `show_constrained_edges=true`

Whether to show the constrained edges.
- `constrained_edge_color=:magenta`

Colour of the constrained edges.
- `constrained_edge_linewidth=2`

Width of the constrained edges.
"""
MakieCore.@recipe(_Triplot, points, triangles, boundary_nodes, convex_hull, constrained_edges, representative_point_list) do scene
    return MakieCore.Attributes(;
        markersize=11,
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
        show_convex_hull=true,
        convex_hull_color=:red,
        convex_hull_linestyle=:dash,
        convex_hull_linewidth=2,
        show_constrained_edges=true,
        constrained_edge_color=:magenta,
        constrained_edge_linewidth=2
    )
end
function MakieCore.convert_arguments(plot::Type{<:_Triplot}, tri::DelaunayTriangulation.Triangulation)
    return (
        DelaunayTriangulation.get_points(tri),
        DelaunayTriangulation.get_triangles(tri),
        DelaunayTriangulation.get_boundary_nodes(tri),
        DelaunayTriangulation.get_convex_hull(tri),
        DelaunayTriangulation.get_all_constrained_edges(tri),
        DelaunayTriangulation.get_representative_point_list(tri)
    )
end

function MakieCore.plot!(p::_Triplot)
    ## Extract
    points = p[:points]
    triangles = p[:triangles]
    boundary_nodes = p[:boundary_nodes]
    convex_hull = p[:convex_hull]
    constrained_edges = p[:constrained_edges]
    representative_point_list = p[:representative_point_list]
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
    show_convex_hull = p[:show_convex_hull]
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
    function update_plot(points, triangles, boundary_nodes, convex_hull, constrained_edges, representative_point_list)
        boundary_map = DelaunayTriangulation.construct_boundary_map(boundary_nodes)
        if recompute_centers[]
            if !DelaunayTriangulation.has_multiple_segments(boundary_nodes) && DelaunayTriangulation.num_boundary_edges(boundary_nodes) == 0
                DelaunayTriangulation.compute_representative_points!(representative_point_list, points, DelaunayTriangulation.get_indices(convex_hull))
            else
                DelaunayTriangulation.compute_representative_points!(representative_point_list, points, boundary_nodes)
            end
        end

        ## Clear out the previous observables
        empty!(points_2f[])
        empty!(triangle_mat[])
        empty!(ghost_edges[])
        empty!(convex_hull_points[])
        empty!(constrained_edge_points[])

        ## Fill out the points 
        for pt in DelaunayTriangulation.each_point(points)
            x, y = DelaunayTriangulation.getxy(pt)
            push!(points_2f[], (x, y))
        end

        ## Fill out the triangles 
        for T in DelaunayTriangulation.each_triangle(triangles)
            if !DelaunayTriangulation.is_ghost_triangle(T)
                i, j, k = DelaunayTriangulation.indices(T)
                push!(triangle_mat[], (i, j, k))
            end
        end
        triangle_mat_2[] = Matrix(reinterpret(reshape, Int64, triangle_mat[])')

        ## Now get the ghost edges if needed 
        if show_ghost_edges[]
            for index in values(boundary_map)
                curve_index = DelaunayTriangulation.get_curve_index(index)
                representative_coordinates = DelaunayTriangulation.get_representative_point_coordinates(representative_point_list, curve_index)
                cx, cy = DelaunayTriangulation.getxy(representative_coordinates)
                bn = !isempty(index) ? DelaunayTriangulation.get_boundary_nodes(boundary_nodes, index) :
                     DelaunayTriangulation.get_indices(convex_hull) # If index is empty, there are no constrained boundaries, meaning we want to get the convex hull
                ## TODO: This needs to get fixed when we use delete_point!, since that will give a convex hull that could end up not 
                ## being right for the ghost edges in the case where we delete a boundary point.
                n_edge = DelaunayTriangulation.num_boundary_edges(bn)
                for i in 1:n_edge
                    u = DelaunayTriangulation.get_boundary_nodes(bn, i)
                    q = DelaunayTriangulation.get_point(points, u)
                    qx, qy = DelaunayTriangulation.getxy(q)
                    if !DelaunayTriangulation.is_interior_curve(curve_index)
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
        if show_convex_hull[]
            idx = DelaunayTriangulation.get_indices(convex_hull)
            for i in idx
                pt = DelaunayTriangulation.get_point(points, i)
                x, y = DelaunayTriangulation.getxy(pt)
                push!(convex_hull_points[], (x, y))
            end
        end

        ## Get the constrained edges if needed 
        if show_constrained_edges[]
            for e in DelaunayTriangulation.each_edge(constrained_edges)
                p1 = DelaunayTriangulation.get_point(points, DelaunayTriangulation.initial(e))
                p2 = DelaunayTriangulation.get_point(points, DelaunayTriangulation.terminal(e))
                push!(constrained_edge_points[], DelaunayTriangulation.getxy(p1), DelaunayTriangulation.getxy(p2))
            end
        end
    end

    ## Connect the plot so that it updates whenever we change a value 
    MakieCore.Observables.onany(update_plot, points, triangles, boundary_nodes, convex_hull, constrained_edges, representative_point_list)

    ## Call it once to prepopulate with current values 
    update_plot(points[], triangles[], boundary_nodes[], convex_hull[], constrained_edges[], representative_point_list[])

    ## Now plot 
    MakieCore.poly!(p, points_2f, triangle_mat_2; strokewidth=strokewidth[], strokecolor=strokecolor[], color=triangle_color[])
    show_all_points[] && MakieCore.scatter!(p, points_2f; markersize=markersize[], color=point_color[])
    show_ghost_edges[] && MakieCore.linesegments!(p, ghost_edges; color=ghost_edge_color[], linewidth=ghost_edge_linewidth[])
    show_convex_hull[] && MakieCore.lines!(p, convex_hull_points; color=convex_hull_color[], linewidth=convex_hull_linewidth[], linestyle=convex_hull_linestyle[])
    show_constrained_edges[] && MakieCore.linesegments!(p, constrained_edge_points; color=constrained_edge_color[], linewidth=constrained_edge_linewidth[])
    return p
end

"""
    _Voronoiplot(!)(vor::VoronoiTessellation)

Plot a Voronoi tessellation.

# Attributes 
- `markersize=11`

Size of the generators.
- `show_generators=true`

Whether to show the generators.
- `generator_color=:black`

Colour of the generators.
- `strokecolor=:black`

Colour of the edges.
- `polygon_color=(:white, 0)`

Colour of the polygons.
- `strokewidth=1`

Width of the edges.
- `unbounded_edge_extension_factor=2.0`

Factor that scales the bounding box for clipping unbounded edges.
- `colormap=th.colormap`

Colormap for the polygons.
- `colorrange=get(th.attributes, :colorrange, MakieCore.automatic)`

Range of the colormap.
- `cycle=[:color]`

Cycle of the colormap.
"""
@recipe(_Voronoiplot, vorn) do scene
    th = MakieCore.default_theme(scene, MakieCore.Mesh)
    return MakieCore.Attributes(;
        markersize=11,
        show_generators=true,
        generator_color=:black,
        strokecolor=:black,
        polygon_color=(:white, 0),
        strokewidth=1,
        unbounded_edge_extension_factor=2.0,
        colormap=th.colormap,
        colorrange=get(th.attributes, :colorrange, MakieCore.automatic),
        cycle=[:color]
    )
end

function MakieCore.plot!(p::_Voronoiplot)
    ## Extract 
    vorn = p[:vorn]
    markersize = p[:markersize]
    show_generators = p[:show_generators]
    generator_color = p[:generator_color]
    strokecolor = p[:strokecolor]
    polygon_color = p[:polygon_color]
    strokewidth = p[:strokewidth]
    unbounded_edge_extension_factor = p[:unbounded_edge_extension_factor]
    colormap = p[:colormap]
    colorrange = p[:colorrange]
    cycle = p[:cycle]
    bbox = DelaunayTriangulation.polygon_bounds(vorn[], unbounded_edge_extension_factor[])
    xmin, xmax, ymin, ymax = bbox
    bbox = [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)]
    bbox_order = [1, 2, 3, 4, 1]

    if !(typeof(polygon_color[]) <: AbstractVector)
        n = DelaunayTriangulation.num_polygons(vorn[])
        polygon_color[] = [polygon_color[] for _ in 1:n]
    end

    ## Define all the necessary observables 
    generators_2f = MakieCore.Observable(NTuple{2,Float64}[])
    polygons = MakieCore.Observable([NTuple{2,Float64}[] for _ in DelaunayTriangulation.each_polygon(vorn[])])

    ## Define the plotting function 
    function update_plot(vorn)
        empty!(generators_2f[])
        empty!(polygons[])
        resize!(generators_2f[], DelaunayTriangulation.num_generators(vorn))
        resize!(polygons[], DelaunayTriangulation.num_generators(vorn))
        for i in DelaunayTriangulation.each_generator(vorn)
            generators_2f[][i] = DelaunayTriangulation.get_generator(vorn, i)
            polygons[][i] = DelaunayTriangulation.get_polygon_coordinates(vorn, i, bbox, bbox_order)
        end
    end

    ## Connect the plot so that it updates whenever we change a value 
    MakieCore.Observables.onany(vorn)

    ## Call it once to prepopulate with current values 
    update_plot(vorn[])

    ## Now plot 
    for i in eachindex(polygons[])
        MakieCore.poly!(p, polygons[][i], color=polygon_color[][i], strokecolor=strokecolor[],
            strokewidth=strokewidth[], colormap=colormap[],
            colorrange=colorrange[], cycle=cycle[])
    end
    if show_generators[]
        MakieCore.scatter!(p, generators_2f[], markersize=markersize[], color=generator_color[])
    end
    return p
end

function DelaunayTriangulation.get_polygon_colors(vorn::DelaunayTriangulation.VoronoiTessellation, cmap)
    F = DelaunayTriangulation.number_type(vorn)
    gtr = [DelaunayTriangulation.get_generator(vorn, i) for i in DelaunayTriangulation.each_generator(vorn)]
    gtr_mat = reinterpret(reshape, F, gtr)
    colors = get(cmap, gtr_mat, :extrema)
    return [(a + b) / 2 for (a, b) in eachcol(colors)]
end

@doc (@doc _triplot) DelaunayTriangulation.triplot(args...; kwargs...) = _triplot(args...; kwargs...)
@doc (@doc _triplot!) DelaunayTriangulation.triplot!(args...; kwargs...) = _triplot!(args...; kwargs...)
@doc (@doc _voronoiplot) DelaunayTriangulation.voronoiplot(args...; kwargs...) = _voronoiplot(args...; kwargs...)
@doc (@doc _voronoiplot!) DelaunayTriangulation.voronoiplot!(args...; kwargs...) = _voronoiplot!(args...; kwargs...)

end # DelaunayTriangulationMakieCoreExt