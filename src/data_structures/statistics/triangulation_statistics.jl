"""
    TriangulationStatistics{T,V,I}

A struct containing statistics about a triangulation. 

# Fields
- `num_vertices::I`: The number of vertices in the triangulation.
- `num_solid_vertices::I`: The number of solid vertices in the triangulation.
- `num_ghost_vertices::I`: The number of ghost vertices in the triangulation.
- `num_edges::I`: The number of edges in the triangulation.
- `num_solid_edges::I`: The number of solid edges in the triangulation.
- `num_ghost_edges::I`: The number of ghost edges in the triangulation.
- `num_triangles::I`: The number of triangles in the triangulation.
- `num_solid_triangles::I`: The number of solid triangles in the triangulation.
- `num_ghost_triangles::I`: The number of ghost triangles in the triangulation.
- `num_boundary_segments::I`: The number of boundary segments in the triangulation.
- `num_interior_segments::I`: The number of interior segments in the triangulation.
- `num_segments::I`: The number of segments in the triangulation.
- `num_convex_hull_vertices::I`: The number of vertices on the convex hull of the triangulation.
- `smallest_angle::V`: The smallest angle in the triangulation.
- `largest_angle::V`: The largest angle in the triangulation.
- `smallest_area::V`: The smallest area of a triangle in the triangulation.
- `largest_area::V`: The largest area of a triangle in the triangulation.
- `smallest_radius_edge_ratio::V`: The smallest radius-edge ratio of a triangle in the triangulation.
- `largest_radius_edge_ratio::V`: The largest radius-edge ratio of a triangle in the triangulation.
- `area::V`: The total area of the triangulation.
- `individual_statistics::Dict{T,IndividualTriangleStatistics{V}}`: A map from triangles in the triangulation to their individual statistics. See [`IndividualTriangleStatistics`](@ref).

# Constructors
To construct these statistics, use [`statistics`](@ref), which you call as `statistics(tri::Triangulation)`.
"""
struct TriangulationStatistics{T, V, I}
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
    individual_statistics::Dict{T, IndividualTriangleStatistics{V}}
end
function Base.show(io::IO, ::MIME"text/plain", stats::TriangulationStatistics)
    println(io, "Delaunay Triangulation Statistics.")
    println(io, "   Triangulation area: ", get_area(stats))
    println(io, "   Number of vertices: ", num_vertices(stats))
    println(io, "   Number of solid vertices: ", num_solid_vertices(stats))
    println(io, "   Number of ghost vertices: ", num_ghost_vertices(stats))
    println(io, "   Number of edges: ", num_edges(stats))
    println(io, "   Number of solid edges: ", num_solid_edges(stats))
    println(io, "   Number of ghost edges: ", num_ghost_edges(stats))
    println(io, "   Number of triangles: ", num_triangles(stats))
    println(io, "   Number of solid triangles: ", num_solid_triangles(stats))
    println(io, "   Number of ghost triangles: ", num_ghost_triangles(stats))
    println(io, "   Number of boundary segments: ", num_boundary_segments(stats))
    println(io, "   Number of interior segments: ", num_interior_segments(stats))
    println(io, "   Number of segments: ", num_segments(stats))
    println(io, "   Number of convex hull vertices: ", num_convex_hull_vertices(stats))
    println(io, "   Smallest angle: ", rad2deg(get_smallest_angle(stats)), "°")
    println(io, "   Largest angle: ", rad2deg(get_largest_angle(stats)), "°")
    println(io, "   Smallest area: ", get_smallest_area(stats))
    println(io, "   Largest area: ", get_largest_area(stats))
    println(io, "   Smallest radius-edge ratio: ", get_smallest_radius_edge_ratio(stats))
    print(io, "   Largest radius-edge ratio: ", get_largest_radius_edge_ratio(stats))
end


"""
    statistics(tri::Triangulation) -> TriangulationStatistics

Returns a [`TriangulationStatistics`](@ref) object containing statistics about the triangulation `tri`.
"""
function statistics(tri::Triangulation)
    F = number_type(tri)
    V = triangle_type(tri)
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
    individual_statistics = Dict{V, IndividualTriangleStatistics{F}}()
    sizehint!(individual_statistics, nsolid_tris)
    smallest_angle = typemax(F)
    largest_angle = typemin(F)
    smallest_area = typemax(F)
    largest_area = typemin(F)
    smallest_radius_edge_ratio = typemax(F)
    largest_radius_edge_ratio = typemin(F)
    total_area = zero(F)
    for T in each_solid_triangle(tri)
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
for n in fieldnames(TriangulationStatistics)
    name = String(n)
    if !contains(name, "num")
        @eval begin
            @doc """
            get_$($(name))(stats::TriangulationStatistics)

        Returns the $($name) field from the [`TriangulationStatistics`](@ref) `stats`.
        """ ($(Symbol("get_$n")))(stats::TriangulationStatistics) = stats.$n
        end
    else
        @eval begin
            @doc """
            $($(name))(stats::TriangulationStatistics)

        Returns the $($name) field from the [TriangulationStatistics`](@ref) `stats`.
        """ ($(Symbol("$n")))(stats::TriangulationStatistics) = stats.$n
        end
    end
end
for n in fieldnames(IndividualTriangleStatistics)
    name = String(n)
    @eval begin
        @doc """
        get_$($(name))(stats::TriangulationStatistics, T)

    Returns the $($name) field from the individual triangle statistics for the triangle `T` in the [`TriangulationStatistics`](@ref) `stats`.
    """ ($(Symbol("get_$n")))(stats::TriangulationStatistics, T) =
            let indiv_stats = get_individual_statistics(stats)
            T = contains_triangle(T, keys(indiv_stats))
            !T[2] && throw(BoundsError(indiv_stats, T))
            return indiv_stats[T[1]].$n
        end
    end
end


"""
    get_minimum_angle(stats::TriangulationStatistics, T) -> Float64

Returns the minimum angle of `T` from `stats`.
"""
get_minimum_angle(stats::TriangulationStatistics, T) = get_angles(stats, T)[1]


"""
    get_maximum_angle(stats::TriangulationStatistics, T) -> Float64

Returns the maximum angle of `T` from `stats`.
"""
get_maximum_angle(stats::TriangulationStatistics, T) = get_angles(stats, T)[3]


"""
    get_median_angle(stats::TriangulationStatistics, T) -> Float64

Returns the median angle of `T` from `stats`.
"""
get_median_angle(stats::TriangulationStatistics, T) = get_angles(stats, T)[2]


"""
    get_all_stat(stats::TriangulationStatistics, stat::Symbol) -> Vector 

Returns a vector of the statistic `stat` for each triangle in `stats`.
"""
function get_all_stat(stats::TriangulationStatistics, stat::Symbol)
    indiv_stats = get_individual_statistics(stats)
    T = first(keys(indiv_stats))
    F = typeof(getfield(indiv_stats[T], stat))
    stats = Vector{F}(undef, length(indiv_stats))
    return get_all_stat!(stats, indiv_stats, stat)
end
function get_all_stat!(stats::Vector{F}, indiv_stats::Dict, stat::Symbol) where {F}
    for (i, T) in enumerate(keys(indiv_stats))
        stats[i] = getfield(indiv_stats[T], stat)::F
    end
    return stats
end