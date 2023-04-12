# Some formulas in this file come from https://perso.uclouvain.be/jean-francois.remacle/LMECA2170/robnotes.pdf.
# Unfortunately, I can't use the robust forms of the formulas because we don't have implementations for them (only the versions that take the sign of the orient predicate).

struct IndividualTriangleStatistics{T}
    area::T
    lengths::NTuple{3,T}
    squared_lengths::NTuple{3,T}
    circumcenter::NTuple{2,T}
    circumradius::T
    sine_minimum_angle::T
    sine_maximum_angle::T
    minimum_angle::T
    maximum_angle::T
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
    sine_minimum_angle = triangle_sine_minimum_angle(A, ℓmed², ℓmax²)
    sine_maximum_angle = triangle_sine_maximum_angle(A, ℓmin², ℓmed²)
    minimum_angle = triangle_minimum_angle(sine_minimum_angle)
    maximum_angle = triangle_maximum_angle(sine_maximum_angle)
    radius_edge_ratio = triangle_radius_edge_ratio(circumradius, ℓmin)
    edge_midpoints = triangle_edge_midpoints(p, q, r)
    perimeter = triangle_perimeter(ℓmin, ℓmed, ℓmax)
    inradius = triangle_inradius(A, perimeter)
    aspect_ratio = triangle_aspect_ratio(inradius, circumradius)
    centroid = triangle_centroid(p, q, r)
    return IndividualTriangleStatistics(
        A,
        (ℓmin, ℓmed, ℓmax),
        (ℓmin², ℓmed², ℓmax²),
        circumcenter,
        circumradius,
        sine_minimum_angle,
        sine_maximum_angle,
        minimum_angle,
        maximum_angle,
        radius_edge_ratio,
        edge_midpoints,
        aspect_ratio,
        inradius,
        perimeter,
        centroid
    )
end

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
    individual_statistics::Dict{T,IndividualTriangleStatistics{V}}
end
function Base.show(io::IO, ::MIME"text/plain", stats::TriangulationStatistics)
    println(io, "Delaunay Triangulation Statistics.")
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
    for T in each_solid_triangle(tri)
        u, v, w = indices(T)
        p, q, r = get_point(tri, u, v, w)
        individual_statistics[T] = IndividualTriangleStatistics(p, q, r)
        smallest_angle = min(smallest_angle, individual_statistics[T].minimum_angle)
        largest_angle = max(largest_angle, individual_statistics[T].maximum_angle)
        smallest_area = min(smallest_area, individual_statistics[T].area)
        largest_area = max(largest_area, individual_statistics[T].area)
        smallest_radius_edge_ratio = min(smallest_radius_edge_ratio, individual_statistics[T].radius_edge_ratio)
        largest_radius_edge_ratio = max(largest_radius_edge_ratio, individual_statistics[T].radius_edge_ratio)
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

## We could simplify some of this even further, e.g. keeping some more things as their square temporarily.
## Is it really worth the effort, though?

function triangle_area(ℓ₁²::Number, ℓ₂²::Number, ℓ₃²::Number)
    A² = squared_triangle_area(ℓ₁², ℓ₂², ℓ₃²)
    A = A² < 0.0 ? zero(A²) : sqrt(A²) # needed e.g. if A² is like -1.084020e-19
    return A
end
squared_triangle_area(ℓ₁²::Number, ℓ₂²::Number, ℓ₃²::Number) = (4ℓ₁² * ℓ₂² - (ℓ₁² + ℓ₂² - ℓ₃²)^2) / 16 # Heron's formula
triangle_minimum_angle(sine_minimum_angle) = asin(sine_minimum_angle)
triangle_maximum_angle(sine_maximum_angle) = asin(sine_maximum_angle)
triangle_circumradius(A, ℓmin², ℓmed², ℓmax²) = sqrt(ℓmin² * ℓmed² * ℓmax²) / (4A)
triangle_perimeter(ℓmin::Number, ℓmed::Number, ℓmax::Number) = ℓmin + ℓmed + ℓmax
triangle_inradius(A, perimeter) = 2A / perimeter
triangle_aspect_ratio(inradius::Number, circumradius::Number) = inradius / circumradius
triangle_radius_edge_ratio(circumradius::Number, ℓmin::Number) = circumradius / ℓmin
triangle_centroid(p, q, r) = ((getx(p) + getx(q) + getx(r)) / 3, (gety(p) + gety(q) + gety(r)) / 3)
triangle_sine_minimum_angle_squared(A, ℓmed², ℓmax²) = 4A^2 / (ℓmed² * ℓmax²)
triangle_sine_maximum_angle_squared(A, ℓmin², ℓmed²) = 4A^2 / (ℓmin² * ℓmed²)
triangle_sine_minimum_angle(A, ℓmed², ℓmax²) = sqrt(triangle_sine_minimum_angle_squared(A, ℓmed², ℓmax²))
triangle_sine_maximum_angle(A, ℓmin², ℓmed²) = sqrt(triangle_sine_maximum_angle_squared(A, ℓmin², ℓmed²))

function squared_triangle_area(p, q, r)
    ℓ₁², ℓ₂², ℓ₃² = squared_triangle_lengths(p, q, r)
    return squared_triangle_area(ℓ₁², ℓ₂², ℓ₃²)
end
function triangle_area(p, q, r)
    A² = squared_triangle_area(p, q, r)
    return sqrt(A²)
end

function squared_triangle_lengths(p, q, r)
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    ℓ₁² = (qx - px)^2 + (qy - py)^2
    ℓ₂² = (rx - qx)^2 + (ry - qy)^2
    ℓ₃² = (px - rx)^2 + (py - ry)^2
    ℓmin², ℓmed², ℓmax² = min_med_max(ℓ₁², ℓ₂², ℓ₃²)
    return ℓmin², ℓmed², ℓmax²
end
function triangle_lengths(p, q, r)
    ℓmin², ℓmed², ℓmax² = squared_triangle_lengths(p, q, r)
    return sqrt(ℓmin²), sqrt(ℓmed²), sqrt(ℓmax²)
end

function triangle_circumcenter(p, q, r, A=triangle_area(p, q, r))
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
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
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    mx = (px + qx) / 2
    my = (py + qy) / 2
    nx = (qx + rx) / 2
    ny = (qy + ry) / 2
    ox = (rx + px) / 2
    oy = (ry + py) / 2
    return (mx, my), (nx, ny), (ox, oy)
end