# While we transition out of the Makie dep, need to put the Makie recipes here for tests 
using Makie
using ..DelaunayTriangulation
const DelTri = DelaunayTriangulation
NEEDS_PLOT_DEFS = !isdefined(Makie, :triplot) && !isdefined(DelaunayTriangulation, :triplot)
if NEEDS_PLOT_DEFS
    @recipe(Triplot, triangles) do scene
        sc = default_theme(scene, Scatter)
        return Attributes(;
            show_points=false,
            show_convex_hull=false,
            show_ghost_edges=false,
            show_constrained_edges=false,
            recompute_centers=false,
            markersize=theme(scene, :markersize),
            marker=theme(scene, :marker),
            strokecolor=theme(scene, :patchstrokecolor),
            strokewidth=1,
            linestyle=:solid,
            triangle_color=(:white, 0.0),
            point_color=sc.color,
            convex_hull_color=:red,
            convex_hull_linestyle=:dash,
            convex_hull_linewidth=theme(scene, :linewidth),
            ghost_edge_color=:blue,
            ghost_edge_linestyle=theme(scene, :linestyle),
            ghost_edge_linewidth=theme(scene, :linewidth),
            ghost_edge_extension_factor=0.1,
            bounding_box=Makie.automatic,
            constrained_edge_color=:magenta,
            constrained_edge_linestyle=theme(scene, :linestyle),
            constrained_edge_linewidth=theme(scene, :linewidth))
    end
    function get_all_triangulation_points!(points, tri)
        empty!(points)
        sizehint!(points, DelTri.num_points(tri))
        for p in DelTri.each_point(tri)
            x, y = DelTri.getxy(p)
            push!(points, Point2f(x, y))
        end
        return points
    end
    function get_present_triangulation_points!(points, tri)
        empty!(points)
        sizehint!(points, DelTri.num_solid_vertices(tri))
        for i in DelTri.each_solid_vertex(tri)
            p = DelTri.get_point(tri, i)
            x, y = DelTri.getxy(p)
            push!(points, Point2f(x, y))
        end
        return points
    end
    function get_triangulation_triangles!(triangles, tri)
        empty!(triangles)
        sizehint!(triangles, DelTri.num_solid_triangles(tri))
        for T in DelTri.each_solid_triangle(tri)
            i, j, k = DelTri.indices(T)
            push!(triangles, Makie.TriangleFace(i, j, k))
        end
        return triangles
    end
    function get_triangulation_ghost_edges!(ghost_edges, extent, tri, bounding_box)
        @assert extent > 0.0 "The ghost_edge_extension_factor must be positive."
        empty!(ghost_edges)
        sizehint!(ghost_edges, 2DelTri.num_ghost_edges(tri))
        if bounding_box === Makie.automatic
            if DelTri.has_boundary_nodes(tri)
                xmin, xmax, ymin, ymax = DelTri.polygon_bounds(DelTri.get_points(tri), DelTri.get_boundary_nodes(tri))
            else
                xmin, xmax, ymin, ymax = DelTri.polygon_bounds(DelTri.get_points(tri),
                    DelTri.get_convex_hull_indices(tri))
            end
            Δx = xmax - xmin
            Δy = ymax - ymin
            a, b, c, d = (xmin - extent * Δx, xmax + extent * Δx, ymin - extent * Δy, ymax + extent * Δy)
        else
            a, b, c, d = bounding_box
        end
        for e in DelTri.each_ghost_edge(tri)
            u, v = DelTri.edge_indices(e)
            if DelTri.is_boundary_index(v)
                u, v = v, u # Make sure that u is the boundary index 
            end
            curve_index = DelTri.get_curve_index(tri, u)
            representative_coordinates = DelTri.get_representative_point_coordinates(tri, curve_index)
            rx, ry = DelTri.getxy(representative_coordinates)
            @assert a ≤ rx ≤ b && c ≤ ry ≤ d "The representative point is not in the bounding box."
            p = DelTri.get_point(tri, v)
            px, py = DelTri.getxy(p)
            if DelTri.is_interior_curve(curve_index)
                ex, ey = rx, ry
            else
                e = DelTri.intersection_of_ray_with_bounding_box(representative_coordinates, p, a, b, c, d)
                ex, ey = DelTri.getxy(e)
            end
            push!(ghost_edges, Point2f(px, py), Point2f(ex, ey))
        end
        return ghost_edges
    end
    function get_triangulation_convex_hull!(convex_hull, tri)
        idx = DelTri.get_convex_hull_indices(tri)
        empty!(convex_hull)
        sizehint!(convex_hull, length(idx))
        for i in idx
            p = DelTri.get_point(tri, i)
            x, y = DelTri.getxy(p)
            push!(convex_hull, Point2f(x, y))
        end
        return convex_hull
    end
    function get_triangulation_constrained_edges!(constrained_edges, tri)
        empty!(constrained_edges)
        sizehint!(constrained_edges, DelTri.num_edges(DelTri.get_all_constrained_edges(tri)))
        for e in DelTri.each_constrained_edge(tri)
            u, v = DelTri.edge_indices(e)
            p = DelTri.get_point(tri, u)
            q = DelTri.get_point(tri, v)
            px, py = DelTri.getxy(p)
            qx, qy = DelTri.getxy(q)
            push!(constrained_edges, Point2f(px, py), Point2f(qx, qy))
        end
        return constrained_edges
    end
    function Makie.plot!(p::Triplot)
        points_2f = Observable(Point2f[])
        present_points_2f = Observable(Point2f[]) # Points might not be in the triangulation yet, so points_2f is not what we want for scatter
        triangles_3f = Observable(Makie.TriangleFace{Int}[])
        ghost_edges_2f = Observable(Point2f[])
        convex_hull_2f = Observable(Point2f[])
        constrained_edges_2f = Observable(Point2f[])
        function update_plot(tri)
            map(p.recompute_centers) do rc
                return rc && DelTri.compute_representative_points!(tri)
            end
            map(points_2f) do pts
                return get_all_triangulation_points!(pts, tri)
            end
            map(p.show_points, present_points_2f) do sp, pts
                return sp && get_present_triangulation_points!(pts, tri)
            end
            map(triangles_3f) do tris
                return get_triangulation_triangles!(tris, tri)
            end
            map(p.show_ghost_edges, p.ghost_edge_extension_factor, ghost_edges_2f, p.bounding_box) do sge, extent, ge, bbox
                return sge && get_triangulation_ghost_edges!(ge, extent, tri, bbox)
            end
            map(p.show_convex_hull, convex_hull_2f) do sch, ch
                return sch && get_triangulation_convex_hull!(ch, tri)
            end
            map(p.show_constrained_edges, constrained_edges_2f) do sce, ce
                return sce && get_triangulation_constrained_edges!(ce, tri)
            end
            for obs in (points_2f, triangles_3f, ghost_edges_2f, convex_hull_2f, constrained_edges_2f)
                notify(obs)
            end
            return nothing
        end
        onany(update_plot, p[1])
        update_plot(p[1][])

        poly!(p, points_2f, triangles_3f; strokewidth=p.strokewidth, strokecolor=p.strokecolor,
            color=p.triangle_color)
        linesegments!(p, ghost_edges_2f; color=p.ghost_edge_color, linewidth=p.ghost_edge_linewidth,
            linestyle=p.ghost_edge_linestyle, xautolimits=false, yautolimits=false)
        lines!(p, convex_hull_2f; color=p.convex_hull_color, linewidth=p.convex_hull_linewidth,
            linestyle=p.convex_hull_linestyle)
        linesegments!(p, constrained_edges_2f; color=p.constrained_edge_color,
            linewidth=p.constrained_edge_linewidth, linestyle=p.constrained_edge_linestyle)
        map(p.show_points) do sp
            return sp &&
                   scatter!(p, present_points_2f; markersize=p.markersize, color=p.point_color,
                strokecolor=p.strokecolor, marker=p.marker)
        end # Do last so that points go over the lines
        return p
    end
    @recipe(Voronoiplot, vorn) do scene
        th = default_theme(scene, Mesh)
        sc = default_theme(scene, Scatter)
        return Attributes(;
            show_generators=true,
            markersize=4,
            marker=sc.marker,
            point_color=sc.color,
            strokecolor=theme(scene, :patchstrokecolor),
            strokewidth=1.0,
            polygon_color=Makie.automatic,
            unbounded_edge_extension_factor=0.1,
            bounding_box=Makie.automatic,
            colormap=th.colormap,
            colorrange=th.colorrange,
            cycle=th.cycle)
    end
    function get_voronoi_tiles!(generators, polygons, vorn, bbox)
        empty!(generators)
        empty!(polygons)
        sizehint!(generators, DelTri.num_generators(vorn))
        sizehint!(polygons, DelTri.num_polygons(vorn))
        for i in DelTri.each_generator(vorn)
            g = DelTri.get_generator(vorn, i)
            x, y = DelTri.getxy(g)
            push!(generators, Point2f(x, y))
            polygon_coords = DelTri.get_polygon_coordinates(vorn, i, bbox)
            polygon_coords_2f = map(polygon_coords) do coords
                x, y = DelTri.getxy(coords)
                return Point2f(x, y)
            end
            push!(polygons, Makie.Polygon(polygon_coords_2f))
        end
        return generators, polygons
    end
    function get_voronoi_colors!(colors, vorn, cmap)
        empty!(colors)
        sizehint!(colors, DelTri.num_polygons(vorn))
        F = DelTri.number_type(vorn)
        gtr = [DelTri.get_generator(vorn, i) for i in DelTri.each_generator(vorn)]
        reverse!(gtr) # For some reason this is needed to get distinct colors for the tiles
        gtr_mat = reinterpret(reshape, F, gtr)
        _colors = get(cgrad(cmap), gtr_mat, :extrema)
        for c in eachcol(_colors)
            a, b = c
            push!(colors, (a + b) / 2)
        end
        return colors
    end
    function Makie.plot!(p::Voronoiplot)
        generators_2f = Observable(Point2f[])
        PolyType = typeof(Makie.Polygon(Point2f[], [Point2f[]]))
        polygons = Observable(PolyType[])
        colors = map(p.polygon_color) do polycol
            if polycol == Makie.automatic
                Makie.RGBA{Float64}[]
            else
                polycol
            end
        end
        function update_plot(vorn)
            bbox = map(p.unbounded_edge_extension_factor, p.bounding_box) do extent, bnd
                if bnd === Makie.automatic
                    return DelTri.polygon_bounds(vorn, extent)
                else
                    return p.bounding_box
                end
            end
            map(generators_2f, polygons, bbox) do gens, polys, box
                return get_voronoi_tiles!(gens, polys, vorn, box)
            end
            map(colors, p.polygon_color, p.colormap) do cols, polycol, cmap
                return polycol == Makie.automatic && get_voronoi_colors!(cols, vorn, cmap)
            end
            for obs in (generators_2f, polygons, colors)
                notify(obs)
            end
        end
        onany(update_plot, p[1])
        update_plot(p[1][])

        poly!(p, polygons; color=colors,
            strokecolor=p.strokecolor,
            strokewidth=p.strokewidth,
            colormap=p.colormap,
            colorrange=p.colorrange,
            cycle=p.cycle)
        map(p.show_generators) do sg
            return sg && scatter!(p, generators_2f;
                markersize=p.markersize,
                marker=p.marker,
                color=p.point_color)
        end
        return p
    end
    function Makie.needs_tight_limits(p::Triplot)
        return p.show_ghost_edges[]
    end
    function Makie.needs_tight_limits(p::Voronoiplot)
        return !isempty(DelTri.get_unbounded_polygons(p[1][]))
    end
end