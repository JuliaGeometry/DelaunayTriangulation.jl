# Some formulas in this file come from https://perso.uclouvain.be/jean-francois.remacle/LMECA2170/robnotes.pdf.
# Unfortunately, I can't use the robust forms of the formulas because we don't have implementations for them (only the versions that take the sign of the orient predicate).

"""
    IndividualTriangleStatistics{T}

Statistics for a single triangle.

# Fields
- `area`: The area of the triangle.
- `lengths`: The sorted lengths of the edges of the triangle.
- `circumcenter`: The circumcenter of the triangle.
- `circumradius`: The circumradius of the triangle.
- `angles`: The sorted angles of the triangle.
- `radius_edge_ratio`: The radius-edge ratio of the triangle.
- `edge_midpoints`: The midpoints of the edges of the triangle.
- `aspect_ratio`: The aspect ratio of the triangle.
- `inradius`: The inradius of the triangle.
- `perimeter`: The perimeter of the triangle.
- `centroid`: The centroid of the triangle.
"""
struct IndividualTriangleStatistics{T}
    area::T
    lengths::NTuple{3,T}
    circumcenter::NTuple{2,T}
    circumradius::T
    angles::NTuple{3, T}
    radius_edge_ratio::T
    edge_midpoints::NTuple{3,NTuple{2,T}}
    aspect_ratio::T
    inradius::T
    perimeter::T
    centroid::NTuple{2,T}
end
function IndividualTriangleStatistics(p, q, r)
    ℓmin², ℓmed², ℓmax² = squared_triangle_lengths(p, q, r)
    ℓmin, ℓmed, ℓmax = sqrt(ℓmin²), sqrt(ℓmed²), sqrt(ℓmax²)
    A = triangle_area(ℓmin², ℓmed², ℓmax²)
    circumcenter = triangle_circumcenter(p, q, r, A)
    circumradius = triangle_circumradius(A, ℓmin², ℓmed², ℓmax²)
    radius_edge_ratio = triangle_radius_edge_ratio(circumradius, ℓmin)
    edge_midpoints = triangle_edge_midpoints(p, q, r)
    perimeter = triangle_perimeter(ℓmin, ℓmed, ℓmax)
    inradius = triangle_inradius(A, perimeter)
    aspect_ratio = triangle_aspect_ratio(inradius, circumradius)
    centroid = triangle_centroid(p, q, r)
    angles = triangle_angles(p, q, r)
    return IndividualTriangleStatistics(
        A,
        (ℓmin, ℓmed, ℓmax),
        circumcenter,
        circumradius,
        angles,
        radius_edge_ratio,
        edge_midpoints,
        aspect_ratio,
        inradius,
        perimeter,
        centroid
    )
end

"""
    TriangulationStatistics{T,V,I}

Statistics for a triangulation.

The constructor for this is `statistics(tri::Triangulation)`.

# Fields 
- `num_vertices`: The number of vertices in the triangulation.
- `num_solid_vertices`: The number of solid vertices in the triangulation.
- `num_ghost_vertices`: The number of ghost vertices in the triangulation.
- `num_edges`: The number of edges in the triangulation.
- `num_solid_edges`: The number of solid edges in the triangulation.
- `num_ghost_edges`: The number of ghost edges in the triangulation.
- `num_triangles`: The number of triangles in the triangulation.
- `num_solid_triangles`: The number of solid triangles in the triangulation.
- `num_ghost_triangles`: The number of ghost triangles in the triangulation.
- `num_constrained_boundary_edges`: The number of constrained boundary edges in the triangulation.
- `num_constrained_interior_edges`: The number of constrained interior edges in the triangulation.
- `num_constrained_edges`: The number of constrained edges in the triangulation.
- `num_convex_hull_points`: The number of points on the convex hull of the triangulation.
- `smallest_angle`: The smallest angle in the triangulation.
- `largest_angle`: The largest angle in the triangulation.
- `smallest_area`: The smallest area in the triangulation.
- `largest_area`: The largest area in the triangulation.
- `smallest_radius_edge_ratio`: The smallest radius-edge ratio in the triangulation.
- `largest_radius_edge_ratio`: The largest radius-edge ratio in the triangulation.
- `total_area`: The total area of the triangulation.
- `individual_statistics`: A dictionary mapping triangles to their individual statistics.
"""
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
    num_constrained_boundary_edges::I
    num_constrained_interior_edges::I
    num_constrained_edges::I
    num_convex_hull_points::I
    smallest_angle::V
    largest_angle::V
    smallest_area::V
    largest_area::V
    smallest_radius_edge_ratio::V
    largest_radius_edge_ratio::V
    total_area::V
    individual_statistics::Dict{T,IndividualTriangleStatistics{V}}
end
function Base.show(io::IO, ::MIME"text/plain", stats::TriangulationStatistics)
    println(io, "Delaunay Triangulation Statistics.")
    println(io, "   Triangulation area: $(get_total_area(stats))")
    println(io, "   Number of vertices: $(num_vertices(stats))")
    println(io, "   Number of solid vertices: $(num_solid_vertices(stats))")
    println(io, "   Number of ghost vertices: $(num_ghost_vertices(stats))")
    println(io, "   Number of edges: $(num_edges(stats))")
    println(io, "   Number of solid edges: $(num_solid_edges(stats))")
    println(io, "   Number of ghost edges: $(num_ghost_edges(stats))")
    println(io, "   Number of triangles: $(num_triangles(stats))")
    println(io, "   Number of solid triangles: $(num_solid_triangles(stats))")
    println(io, "   Number of ghost triangles: $(num_ghost_triangles(stats))")
    println(io, "   Number of constrained boundary edges: $(num_constrained_boundary_edges(stats))")
    println(io, "   Number of constrained interior edges: $(num_constrained_interior_edges(stats))")
    println(io, "   Number of constrained edges: $(num_constrained_edges(stats))")
    println(io, "   Number of convex hull points: $(num_convex_hull_points(stats))")
    println(io, "   Smallest angle: $(rad2deg(get_smallest_angle(stats)))°")
    println(io, "   Largest angle: $(rad2deg(get_largest_angle(stats)))°")
    println(io, "   Smallest area: $(get_smallest_area(stats))")
    println(io, "   Largest area: $(get_largest_area(stats))")
    println(io, "   Smallest radius-edge ratio: $(get_smallest_radius_edge_ratio(stats))")
    print(io, "   Largest radius-edge ratio: $(get_largest_radius_edge_ratio(stats))")
end

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
    constrained_boundary_edges = keys(get_boundary_edge_map(tri))
    nconstrained_boundary_edges = length(constrained_boundary_edges)
    constrained_interior_edges = get_constrained_edges(tri)
    nconstrained_interior_edges = num_edges(constrained_interior_edges)
    constrained_edges = get_all_constrained_edges(tri)
    nconstrained_edges = num_edges(constrained_edges)
    convex_hull_indices = get_convex_hull_indices(tri)
    nconvex_hull_points = max(0, length(convex_hull_indices) - 1) # -1 because the last index is the same as the first 
    individual_statistics = Dict{V,IndividualTriangleStatistics{F}}()
    sizehint!(individual_statistics, nsolid_tris)
    smallest_angle = typemax(F)
    largest_angle = typemin(F)
    smallest_area = typemax(F)
    largest_area = typemin(F)
    smallest_radius_edge_ratio = typemax(F)
    largest_radius_edge_ratio = typemin(F)
    total_area = zero(F)
    for T in each_solid_triangle(tri)
        u, v, w = indices(T)
        p, q, r = get_point(tri, u, v, w)
        individual_statistics[T] = IndividualTriangleStatistics(p, q, r)
        smallest_angle = min(smallest_angle, individual_statistics[T].angles[1])
        largest_angle = max(largest_angle, individual_statistics[T].angles[3])
        smallest_area = min(smallest_area, individual_statistics[T].area)
        largest_area = max(largest_area, individual_statistics[T].area)
        smallest_radius_edge_ratio = min(smallest_radius_edge_ratio, individual_statistics[T].radius_edge_ratio)
        largest_radius_edge_ratio = max(largest_radius_edge_ratio, individual_statistics[T].radius_edge_ratio)
        total_area += individual_statistics[T].area
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
        I(nconstrained_boundary_edges),
        I(nconstrained_interior_edges),
        I(nconstrained_edges),
        I(nconvex_hull_points),
        smallest_angle,
        largest_angle,
        smallest_area,
        largest_area,
        smallest_radius_edge_ratio,
        largest_radius_edge_ratio,
        total_area,
        individual_statistics
    )
end
for n in fieldnames(TriangulationStatistics)
    name = String(n)
    if !contains(name, "num")
        @eval begin
            @doc """
            get_$($(name))(stats::TriangulationStatistics)

        Returns the $($name) field from the triangulation statistics `stats`.
        """ ($(Symbol("get_$n")))(stats::TriangulationStatistics) = stats.$n
        end
    else
        @eval begin
            @doc """
            $($(name))(stats::TriangulationStatistics)

        Returns the $($name) field from the triangulation statistics `stats`.
        """ ($(Symbol("$n")))(stats::TriangulationStatistics) = stats.$n
        end
    end
end
for n in fieldnames(IndividualTriangleStatistics)
    name = String(n)
    @eval begin
        @doc """
        get_$($(name))(stats::TriangulationStatistics, T)

    Returns the $($name) field from the individual triangle statistics for the triangle `T` in the triangulation statistics `stats`.
    """ ($(Symbol("get_$n")))(stats::TriangulationStatistics, T) =
            let indiv_stats = get_individual_statistics(stats)
                T = contains_triangle(T, keys(indiv_stats))
                !T[2] && throw(BoundsError(indiv_stats, T))
                return indiv_stats[T[1]].$n
            end
    end
end
"""
    get_minimum_angle(stats::TriangulationStatistics, T)

Returns the smallest angle in the triangle `T`.
"""
get_minimum_angle(stats::TriangulationStatistics, T) = get_angles(stats, T)[1]
"""
    get_maximum_angle(stats::TriangulationStatistics, T)

Returns the largest angle in the triangle `T`.
"""
get_maximum_angle(stats::TriangulationStatistics, T) = get_angles(stats, T)[3]
"""
    get_median_angle(stats::TriangulationStatistics, T)

Returns the median angle in the triangle `T`.
"""
get_median_angle(stats::TriangulationStatistics, T) = get_angles(stats, T)[2]

"""
    get_all_stat(stats::TriangulationStatistics, stat::Symbol)

Returns all of the values of the statistic `stat` for all triangles in the triangulation statistics `stats`. The possible 
values for `stat` come from [`IndividualTriangleStatistics`](@ref), namely:

- `area`
- `lengths`
- `circumcenter`
- `circumradius`
- `angles`
- `radius_edge_ratio`
- `edge_midpoints`
- `aspect_ratio`
- `inradius`
- `perimeter`
- `centroid`
"""
function get_all_stat(stats::TriangulationStatistics, stat::Symbol)
    indiv_stats = get_individual_statistics(stats)
    return [getfield(indiv_stats[T], stat) for T in keys(indiv_stats)]
end

## We could simplify some of this even further, e.g. keeping some more things as their square temporarily.
## Is it really worth the effort, though?

function triangle_area(ℓ₁²::Number, ℓ₂²::Number, ℓ₃²::Number)
    A² = squared_triangle_area(ℓ₁², ℓ₂², ℓ₃²)
    A = A² < 0.0 ? zero(A²) : sqrt(A²) # needed e.g. if A² is like -1.084020e-19
    return A
end
squared_triangle_area(ℓ₁²::Number, ℓ₂²::Number, ℓ₃²::Number) = (4ℓ₁² * ℓ₂² - (ℓ₁² + ℓ₂² - ℓ₃²)^2) / 16 # Heron's formula
triangle_circumradius(A, ℓmin², ℓmed², ℓmax²) = sqrt(ℓmin² * ℓmed² * ℓmax²) / (4A)
triangle_perimeter(ℓmin::Number, ℓmed::Number, ℓmax::Number) = ℓmin + ℓmed + ℓmax
triangle_inradius(A, perimeter) = 2A / perimeter
triangle_aspect_ratio(inradius::Number, circumradius::Number) = inradius / circumradius
triangle_radius_edge_ratio(circumradius::Number, ℓmin::Number) = circumradius / ℓmin
triangle_centroid(p, q, r) = ((_getx(p) + _getx(q) + _getx(r)) / 3, (_gety(p) + _gety(q) + _gety(r)) / 3)

function triangle_angles(p, q, r)
    ℓ₁², ℓ₂², ℓ₃² = squared_triangle_lengths(p, q, r)
    A = triangle_area(ℓ₁², ℓ₂², ℓ₃²)
    px, py = _getxy(p)
    qx, qy = _getxy(q)
    rx, ry = _getxy(r)
    ax, by = px - qx, py - qy
    bx, ay = px - rx, py - ry
    dotab = ax * bx + ay * by
    θ₁ = if iszero(dotab)
        one(dotab) * π / 2 
    else 
        atan(2A / dotab) 
    end
    if θ₁ < 0
        θ₁ += π
    end
    ax, by = qx - px, qy - py
    bx, ay = qx - rx, qy - ry
    dotab = ax * bx + ay * by
    θ₂ = if iszero(dotab)
        one(dotab) * π / 2 
    else 
        atan(2A / dotab)
    end
    if θ₂ < 0
        θ₂ += π
    end
    ax, by = rx - px, ry - py
    bx, ay = rx - qx, ry - qy
    dotab = ax * bx + ay * by
    θ₃ = if iszero(dotab)
        one(dotab) * π / 2 
    else 
        atan(2A / dotab) 
    end
    if θ₃ < 0
        θ₃ += π
    end
    θ₁, θ₂, θ₃ = min_med_max(θ₁, θ₂, θ₃)
    return θ₁, θ₂, θ₃
end

function squared_triangle_area(p, q, r)
    ℓ₁², ℓ₂², ℓ₃² = squared_triangle_lengths(p, q, r)
    return squared_triangle_area(ℓ₁², ℓ₂², ℓ₃²)
end
function triangle_area(p, q, r)
    A² = squared_triangle_area(p, q, r)
    return A² < zero(A²) ? zero(A²) : sqrt(A²)
end

function squared_triangle_lengths(p, q, r)
    ℓ₁², ℓ₂², ℓ₃², _ = squared_triangle_lengths_and_smallest_index(p, q, r)
    return ℓ₁², ℓ₂², ℓ₃²
end

function squared_triangle_lengths_and_smallest_index(p, q, r)
    px, py = _getxy(p)
    qx, qy = _getxy(q)
    rx, ry = _getxy(r)
    ℓ₁² = (qx - px)^2 + (qy - py)^2
    ℓ₂² = (rx - qx)^2 + (ry - qy)^2
    ℓ₃² = (px - rx)^2 + (py - ry)^2
    ℓmin², ℓmed², ℓmax² = min_med_max(ℓ₁², ℓ₂², ℓ₃²)
    ℓmin² == ℓ₁² && return ℓmin², ℓmed², ℓmax², 1
    ℓmin² == ℓ₂² && return ℓmin², ℓmed², ℓmax², 2
    ℓmin² == ℓ₃² && return ℓmin², ℓmed², ℓmax², 3
end

function triangle_lengths(p, q, r)
    ℓmin², ℓmed², ℓmax² = squared_triangle_lengths(p, q, r)
    return sqrt(ℓmin²), sqrt(ℓmed²), sqrt(ℓmax²)
end

function triangle_circumcenter(p, q, r, A=triangle_area(p, q, r))
    px, py = _getxy(p)
    qx, qy = _getxy(q)
    rx, ry = _getxy(r)
    d11 = (px - rx)^2 + (py - ry)^2
    d12 = py - ry
    d21 = (qx - rx)^2 + (qy - ry)^2
    d22 = qy - ry
    ox = rx + (d11 * d22 - d12 * d21) / (4A)
    e11 = px - rx
    e12 = d11
    e21 = qx - rx
    e22 = d21
    oy = ry + (e11 * e22 - e12 * e21) / (4A)
    return (ox, oy)
end
function triangle_circumcenter(tri::Triangulation, T)
    i, j, k = indices(T)
    p, q, r = get_point(tri, i, j, k)
    return triangle_circumcenter(p, q, r)
end

function triangle_circumradius(p, q, r)
    ℓ₁², ℓ₂², ℓ₃² = squared_triangle_lengths(p, q, r)
    A = triangle_area(ℓ₁², ℓ₂², ℓ₃²)
    return triangle_circumradius(A, ℓ₁², ℓ₂², ℓ₃²)
end

function triangle_perimeter(p, q, r)
    ℓmin, ℓmed, ℓmax = triangle_lengths(p, q, r)
    return triangle_perimeter(ℓmin, ℓmed, ℓmax)
end

function triangle_inradius(p, q, r)
    ℓmin², ℓmed², ℓmax² = squared_triangle_lengths(p, q, r)
    ℓmin, ℓmed, ℓmax = sqrt(ℓmin²), sqrt(ℓmed²), sqrt(ℓmax²)
    A = triangle_area(ℓmin², ℓmed², ℓmax²)
    perimeter = triangle_perimeter(ℓmin, ℓmed, ℓmax)
    return triangle_inradius(A, perimeter)
end

function triangle_aspect_ratio(p, q, r)
    ℓmin², ℓmed², ℓmax² = squared_triangle_lengths(p, q, r)
    ℓmin, ℓmed, ℓmax = sqrt(ℓmin²), sqrt(ℓmed²), sqrt(ℓmax²)
    A = triangle_area(ℓmin², ℓmed², ℓmax²)
    perimeter = triangle_perimeter(ℓmin, ℓmed, ℓmax)
    inradius = triangle_inradius(A, perimeter)
    circumradius = triangle_circumradius(A, ℓmin², ℓmed², ℓmax²)
    return triangle_aspect_ratio(inradius, circumradius)
end

function triangle_radius_edge_ratio(p, q, r)
    ℓmin², ℓmed², ℓmax² = squared_triangle_lengths(p, q, r)
    ℓmin = sqrt(ℓmin²)
    A = triangle_area(ℓmin², ℓmed², ℓmax²)
    circumradius = triangle_circumradius(A, ℓmin², ℓmed², ℓmax²)
    return triangle_radius_edge_ratio(circumradius, ℓmin)
end

function triangle_edge_midpoints(p, q, r)
    px, py = _getxy(p)
    qx, qy = _getxy(q)
    rx, ry = _getxy(r)
    mx = (px + qx) / 2
    my = (py + qy) / 2
    nx = (qx + rx) / 2
    ny = (qy + ry) / 2
    ox = (rx + px) / 2
    oy = (ry + py) / 2
    return (mx, my), (nx, ny), (ox, oy)
end

function get_total_area(tri::Triangulation)
    F = number_type(tri)
    A = zero(F)
    for T in each_solid_triangle(tri)
        u,v,w=indices(T)
        p,q,r=get_point(tri,u,v,w)
        A += triangle_area(p, q, r)
    end
    return A
end