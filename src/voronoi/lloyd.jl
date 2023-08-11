function move_generator_to_centroid!(points, vorn::VoronoiTessellation, generator)
    c = get_centroid(vorn, generator)
    cx, cy = _getxy(c)
    p = get_generator(vorn, generator)
    px, py = _getxy(p)
    dist = sqrt((cx - px)^2 + (cy - py)^2)
    set_point!(points, generator, cx, cy)
    return dist
end

function default_displacement_tolerance(vorn::VoronoiTessellation)
    xmin, xmax, ymin, ymax = polygon_bounds(vorn)
    max_extent = max(xmax - xmin, ymax - ymin)
    return 1e-4max_extent
end

function _centroidal_smooth_itr(vorn::VoronoiTessellation, set_boundary_nodes, points, edges, boundary_nodes,
    I, E, V, Es, Ts, F, rng; kwargs...)
    max_dist = zero(F)
    for i in each_generator(vorn)
        if i ∉ set_boundary_nodes
            dist = move_generator_to_centroid!(points, vorn, i)
            max_dist = max(max_dist, dist)
        end
    end
    _tri = triangulate(points; edges, boundary_nodes,
        IntegerType=I, EdgeType=E, TriangleType=V, EdgesType=Es, TrianglesType=Ts, delete_ghosts=false,
        delete_empty_features=false, check_arguments=false, rng=rng, kwargs...)
    vorn = voronoi(_tri, true)
    return vorn, max_dist
end
