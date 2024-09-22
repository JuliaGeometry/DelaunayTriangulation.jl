"""
    SphericalDelaunay

This module contains functions for computing Delaunay triangulations and Voronoi tessellations on the 
surface of a sphere. Currently, all spheres are assumed to be centered at the origin and have a unit radius.

This module is not yet entirely developed, and so you should consider it experimental.
"""
module SphericalDelaunay

import ..DelaunayTriangulation:
    number_type,
    getxyz,
    getxy,
    getx,
    gety,
    get_voronoi_vertex,
    Triangulation,
    triangle_circumcenter,
    get_points,
    triangle_vertices,
    sort_triangle,
    triangle_area,
    triangle_circumradius,
    det,
    each_point,
    triangulate,
    set_point!,
    push_point!,
    convex_hull!,
    AbstractParametricCurve,
    norm_sqr,
    get_representative_point_list,
    RepresentativeCoordinates,
    integer_type,
    empty_representative_points!,
    new_representative_point!,
    delete_ghost_triangles!,
    num_points,
    delete_ghost_vertices_from_graph!,
    each_ghost_triangle,
    add_triangle!,
    each_solid_triangle,
    num_solid_triangles,
    num_triangles,
    each_triangle,
    is_planar,
    _quality_statistics,
    RefinementArguments,
    triangle_radius_edge_ratio,
    triangle_sink,
    triangle_centroid,
    AbstractPredicateKernel,
    num_solid_vertices,
    AdaptiveKernel,
    get_steiner_point,
    Certificate,
    enqueue_triangle!,
    RefinementQueue,
    RefinementConstraints,
    _build_queues,
    check_steiner_point_precision,
    finalise!,
    _each_solid_triangle,
    has_ghost_triangles,
    _is_ghost_triangle,
    locate_steiner_point,
    voronoi,
    VoronoiTessellation,
    is_ghost_vertex,
    get_circumcenter_to_triangle,
    push_polygon_point!,
    get_triangle_to_circumcenter,
    num_polygons,
    get_polygon_points,
    get_triangulation,
    num_polygon_vertices,
    each_polygon,
    each_polygon_index,
    get_polygon,
    get_polygon_point

import Random:
    Random,
    rand,
    AbstractRNG,
    SamplerTrivial

export SphericalPoint, 
    SphericalTriangulation, 
    SphericalTessellation,
    spherical_triangulate,
    spherical_voronoi,  
    UnitSphere,
    SphericalLine,
    SphericalTriangle,
    SphericalPolygon,
    solidify!

include("geometry.jl")
include("utils.jl")
include("spherical_triangles.jl")
include("triangulation.jl")
include("refine.jl")
include("voronoi.jl")
include("plotting.jl")

end # module