struct TriangulationStatistics{T,V,I}
    num_vertices::I
    num_solid_vertices::I
    num_ghost_vertices::I
    num_edges::I
    num_solid_edges::I
    num_ghost_edges::I
    num_triangles::I
    num_solid_triangles::I
    num_ghost_triangles::I
    num_boundary_segments::I
    num_interior_segments::I
    num_segments::I
    num_convex_hull_vertices::I
    smallest_angle::V
    largest_angle::V
    smallest_area::V
    largest_area::V
    smallest_radius_edge_ratio::V
    largest_radius_edge_ratio::V
    area::V
    individual_statistics::Dict{T,IndividualTriangleStatistics{V}}
end

"""
    statistics(tri::Triangulation) -> TriangulationStatistics

Returns a [`TriangulationStatistics`](@ref) object containing statistics about the triangulation `tri`.
"""
function statistics(tri::Triangulation)
    F = number_type(tri)
    I = integer_type(tri)
    nverts = num_vertices(tri)
    nsolid_verts = num_solid_vertices(tri)
    nghost_verts = num_ghost_vertices(tri)
    nedges = num_edges(tri)
    nsolid_edges = num_solid_edges(tri)
    nghost_edges = num_ghost_edges(tri)
    ntris = num_triangles(tri)
    nsolid_tris = num_solid_triangles(tri)
    nghost_tris = num_ghost_triangles(tri)
    boundary_segments = each_boundary_edge(tri)
    nboundary_segments = num_edges(boundary_segments)
    interior_segments = get_interior_segments(tri)
    ninterior_segments = num_edges(interior_segments)
    segments = get_all_segments(tri)
    nsegments = num_edges(segments)
    convex_hull_vertices = get_convex_hull_vertices(tri)
    nconvex_hull_vertices = max(0, length(convex_hull_vertices) - 1) # -1 because the last index is the same as the first 
    individual_statistics = Dict{NTuple{3,I},IndividualTriangleStatistics{F}}()
    sizehint!(individual_statistics, nsolid_tris)
    smallest_angle = typemax(F)
    largest_angle = typemin(F)
    smallest_area = typemax(F)
    largest_area = typemin(F)
    smallest_radius_edge_ratio = typemax(F)
    largest_radius_edge_ratio = typemin(F)
    total_area = zero(F)
    for T in each_solid_triangle(tri)
        T = sort_triangle(T)
        u, v, w = triangle_vertices(T)
        p, q, r = get_point(tri, u, v, w)
        individual_statistics[T] = IndividualTriangleStatistics(p, q, r, triangle_sink(tri, T))
        smallest_angle = min(F(smallest_angle), individual_statistics[T].angles[1])
        largest_angle = max(F(largest_angle), individual_statistics[T].angles[3])
        smallest_area = min(F(smallest_area), individual_statistics[T].area)
        largest_area = max(F(largest_area), individual_statistics[T].area)
        smallest_radius_edge_ratio = min(F(smallest_radius_edge_ratio), individual_statistics[T].radius_edge_ratio)
        largest_radius_edge_ratio = max(F(largest_radius_edge_ratio), individual_statistics[T].radius_edge_ratio)
        total_area += F(individual_statistics[T].area)
    end
    return TriangulationStatistics(
        I(nverts),
        I(nsolid_verts),
        I(nghost_verts),
        I(nedges),
        I(nsolid_edges),
        I(nghost_edges),
        I(ntris),
        I(nsolid_tris),
        I(nghost_tris),
        I(nboundary_segments),
        I(ninterior_segments),
        I(nsegments),
        I(nconvex_hull_vertices),
        smallest_angle,
        largest_angle,
        smallest_area,
        largest_area,
        smallest_radius_edge_ratio,
        largest_radius_edge_ratio,
        total_area,
        individual_statistics,
    )
end

gen_stat_doc(name, prefix="") = """
    $(prefix)$(isempty(prefix) ? prefix : "")$name(stats::TriangulationStatistics)

Returns the `$name` field from the `stats`.
"""

@doc gen_stat_doc("num_vertices") @inline num_vertices(stats::TriangulationStatistics) = stats.num_vertices
@doc gen_stat_doc("num_solid_vertices") @inline num_solid_vertices(stats::TriangulationStatistics) = stats.num_solid_vertices
@doc gen_stat_doc("num_ghost_vertices") @inline num_ghost_vertices(stats::TriangulationStatistics) = stats.num_ghost_vertices
@doc gen_stat_doc("num_edges") @inline num_edges(stats::TriangulationStatistics) = stats.num_edges
@doc gen_stat_doc("num_solid_edges") @inline num_solid_edges(stats::TriangulationStatistics) = stats.num_solid_edges
@doc gen_stat_doc("num_ghost_edges") @inline num_ghost_edges(stats::TriangulationStatistics) = stats.num_ghost_edges
@doc gen_stat_doc("num_triangles") @inline num_triangles(stats::TriangulationStatistics) = stats.num_triangles
@doc gen_stat_doc("num_solid_triangles") @inline num_solid_triangles(stats::TriangulationStatistics) = stats.num_solid_triangles
@doc gen_stat_doc("num_ghost_triangles") @inline num_ghost_triangles(stats::TriangulationStatistics) = stats.num_ghost_triangles
@doc gen_stat_doc("num_boundary_segments") @inline num_boundary_segments(stats::TriangulationStatistics) = stats.num_boundary_segments
@doc gen_stat_doc("num_interior_segments") @inline num_interior_segments(stats::TriangulationStatistics) = stats.num_interior_segments
@doc gen_stat_doc("num_segments") @inline num_segments(stats::TriangulationStatistics) = stats.num_segments
@doc gen_stat_doc("num_convex_hull_vertices") @inline num_convex_hull_vertices(stats::TriangulationStatistics) = stats.num_convex_hull_vertices
@doc gen_stat_doc("smallest_angle", "get") @inline get_smallest_angle(stats::TriangulationStatistics) = stats.smallest_angle
@doc gen_stat_doc("largest_angle", "get") @inline get_largest_angle(stats::TriangulationStatistics) = stats.largest_angle
@doc gen_stat_doc("smallest_area", "get") @inline get_smallest_area(stats::TriangulationStatistics) = stats.smallest_area
@doc gen_stat_doc("largest_area", "get") @inline get_largest_area(stats::TriangulationStatistics) = stats.largest_area
@doc gen_stat_doc("smallest_radius_edge_ratio", "get") @inline get_smallest_radius_edge_ratio(stats::TriangulationStatistics) = stats.smallest_radius_edge_ratio
@doc gen_stat_doc("largest_radius_edge_ratio", "get") @inline get_largest_radius_edge_ratio(stats::TriangulationStatistics) = stats.largest_radius_edge_ratio
@doc gen_stat_doc("area", "get") @inline get_area(stats::TriangulationStatistics) = stats.area

@inline get_individual_statistics(stats::TriangulationStatistics) = stats.individual_statistics
@inline get_individual_statistics(stats::TriangulationStatistics, T) = get_individual_statistics(stats)[sort_triangle(T)]


gen_istat_doc(name) = """
    get_$name(stats::TriangulationStatistics, T)

Returns the `$name` field from the individual triangle statistics for the triangle `T` in the `stats`.
"""

@doc gen_istat_doc("area") @inline get_area(stats::TriangulationStatistics, T) = get_area(get_individual_statistics(stats, T))
@doc gen_istat_doc("lengths") @inline get_lengths(stats::TriangulationStatistics, T) = get_lengths(get_individual_statistics(stats, T))
@doc gen_istat_doc("circumcenter") @inline get_circumcenter(stats::TriangulationStatistics, T) = get_circumcenter(get_individual_statistics(stats, T))
@doc gen_istat_doc("circumradius") @inline get_circumradius(stats::TriangulationStatistics, T) = get_circumradius(get_individual_statistics(stats, T))
@doc gen_istat_doc("angles") @inline get_angles(stats::TriangulationStatistics, T) = get_angles(get_individual_statistics(stats, T))
@doc gen_istat_doc("radius_edge_ratio") @inline get_radius_edge_ratio(stats::TriangulationStatistics, T) = get_radius_edge_ratio(get_individual_statistics(stats, T))
@doc gen_istat_doc("edge_midpoints") @inline get_edge_midpoints(stats::TriangulationStatistics, T) = get_edge_midpoints(get_individual_statistics(stats, T))
@doc gen_istat_doc("aspect_ratio") @inline get_aspect_ratio(stats::TriangulationStatistics, T) = get_aspect_ratio(get_individual_statistics(stats, T))
@doc gen_istat_doc("inradius") @inline get_inradius(stats::TriangulationStatistics, T) = get_inradius(get_individual_statistics(stats, T))
@doc gen_istat_doc("perimeter") @inline get_perimeter(stats::TriangulationStatistics, T) = get_perimeter(get_individual_statistics(stats, T))
@doc gen_istat_doc("centroid") @inline get_centroid(stats::TriangulationStatistics, T) = get_centroid(get_individual_statistics(stats, T))
@doc gen_istat_doc("offcenter") @inline get_offcenter(stats::TriangulationStatistics, T) = get_offcenter(get_individual_statistics(stats, T))
@doc gen_istat_doc("sink") @inline get_sink(stats::TriangulationStatistics, T) = get_sink(get_individual_statistics(stats, T))

@inline get_minimum_angle(stats::TriangulationStatistics, T) = get_angles(stats, T)[1]
@inline get_maximum_angle(stats::TriangulationStatistics, T) = get_angles(stats, T)[3]

@inline function get_all_stat(stats::TriangulationStatistics, stat::Symbol)
    indiv_stats = get_individual_statistics(stats)
    T = first(keys(indiv_stats))
    F = typeof(getfield(indiv_stats[T], stat))
    stats = Vector{F}(undef, length(indiv_stats))
    return get_all_stat!(stats, indiv_stats, stat)
end
@inline function get_all_stat!(stats::Vector{F}, indiv_stats::Dict, stat::Symbol) where {F}
    for (i, T) in enumerate(keys(indiv_stats))
        stats[i] = getfield(indiv_stats[T], stat)::F
    end
    return stats
end
