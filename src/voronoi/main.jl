"""
    voronoi(tri::Triangulation, clip=has_boundary_nodes(tri)) -> VoronoiTessellation

Construct the Voronoi tessellation of the points in `tri`. If `clip` is `true` then the
Voronoi tessellation will be clipped to the convex hull of the points in `tri`. 

!!! warning 

    Exact predicates are used only for classifying intersections, but no special methods 
    are used for computing the intersections themselves. 

!!! warning 

    Clipping is not yet guaranteed to work for constrained triangulations. Unfortunate,
    because that's really what I actually care about. If you have any thoughts for that, 
    or just any thoughts on the reference clipping paper (by Yan, Wang, Levy, and Liu: 
    "Efficient Computation of Clipped Voronoi Diagram for Mesh Generation), please let me know via 
    an issue. There is certainly a lot of room for improvement in the clipping code too.
"""
function voronoi(tri::Triangulation, clip=has_boundary_nodes(tri))
    has_ghost = has_ghost_triangles(tri)
    !has_ghost && add_ghost_triangles!(tri)
    vorn = initialise_voronoi_tessellation(tri)
    has_bnds = has_boundary_nodes(tri)
    is_convex = has_bnds
    is_true(clip) && !has_bnds && lock_convex_hull!(tri)
    for i in each_generator(vorn)
        add_voronoi_polygon!(vorn, i)
    end
    is_true(clip) && clip_voronoi_tessellation!(vorn, is_convex)
    is_true(clip) && !has_bnds && unlock_convex_hull!(tri)
    !has_ghost && delete_ghost_triangles!(tri)
    return vorn
end

"""
    clip_voronoi_tessellation!(vorn::VoronoiTessellation)

Clip the Voronoi tessellation `vorn` to the convex hull of the points in `vorn`. 
"""
function clip_voronoi_tessellation!(vorn::VoronoiTessellation, is_convex=true)
    boundary_sites, segment_intersections, exterior_circumcenters, equal_circumcenter_mapping = find_all_intersections(vorn)
    n = add_intersection_points!(vorn, segment_intersections)
    clip_all_polygons!(vorn, n, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping, is_convex)
    add_all_boundary_polygons!(vorn, boundary_sites)
    return nothing
end

"""
    centroidal_smooth(vorn::VoronoiTessellation; maxiters=1000, tol=default_displacement_tolerance(vorn), rng=Random.default_rng(), kwargs...) -> VoronoiTessellation

Smooth the Voronoi tessellation `vorn` by moving each generator to the centroid of its Voronoi polygon. Uses Lloyd's algorithm.
The smoothing will stop when the maximum displacement of any generator is less than `tol` or when the maximum number of iterations `maxiters` is reached.

The refinement is not done in place, and instead the smoothed tessellation is returned with no modifications to `vorn`.

This function only works on clipped tessellations, and generators on the boundary of the underlying triangulation are fixed in-place.
"""
function centroidal_smooth(vorn::VoronoiTessellation{Tr}; maxiters=1000, tol=default_displacement_tolerance(vorn), rng=Random.default_rng(), kwargs...) where {P,Ts,I,E,Es,BN,BNM,B,BIR,BPL,Tr<:Triangulation{P,Ts,I,E,Es,BN,BNM,B,BIR,BPL}}
    iter = 0
    F = number_type(vorn)
    max_dist = typemax(F)
    tri = get_triangulation(vorn)
    has_ghost = has_ghost_triangles(tri)
    !has_ghost && add_ghost_triangles!(tri)
    has_bnds = has_boundary_nodes(tri)
    !has_bnds && lock_convex_hull!(tri)
    set_boundary_nodes = get_all_boundary_nodes(tri)
    points = (deepcopy ∘ get_points)(tri)
    boundary_nodes = get_boundary_nodes(tri)
    edges = get_constrained_edges(tri)
    if isempty(edges)
        edges = nothing
    end
    V = triangle_type(Ts)
    while iter < maxiters && max_dist > tol
        vorn, max_dist = _centroidal_smooth_itr(vorn, set_boundary_nodes, points, edges, boundary_nodes,
            I, E, V, Es, Ts, F, rng; kwargs...)
        iter += 1
    end
    !has_bnds && unlock_convex_hull!(tri)
    !has_ghost && delete_ghost_triangles!(tri)
    !has_bnds && empty!(get_all_constrained_edges(get_triangulation(vorn)))
    !has_ghost && delete_ghost_triangles!(get_triangulation(vorn))
    return vorn
end