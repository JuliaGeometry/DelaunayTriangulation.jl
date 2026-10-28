"""
    initialise_voronoi_tessellation(tri::Triangulation, predicates::AbstractPredicateKernel=AdaptiveKernel()) -> VoronoiTessellation

Initialise a `VoronoiTessellation` from the triangulation `tri`.

# Arguments 
- `tri`: The [`Triangulation`](@ref).
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates, used for finding triangles with the same circumcircle; see [`assign_cocircular_neighbours!`](@ref). Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Output
- `vorn`: The [`VoronoiTessellation`](@ref). This tessellation is not yet filled in, as all the polygons and other fields need to be properly defined. This simply defines all the initial objects that will be pushed into.
"""
function initialise_voronoi_tessellation(tri::Tr, predicates::AbstractPredicateKernel = AdaptiveKernel()) where {Tr <: Triangulation}
    I = integer_type(tri)
    T = triangle_type(tri)
    F = number_type(tri)
    P = NTuple{2, F}
    polygon_points = Vector{P}()
    circumcenter_to_triangle = Dict{I, T}()
    triangle_to_circumcenter = Dict{T, I}()
    sizehint!(polygon_points, num_triangles(tri))
    sizehint!(circumcenter_to_triangle, num_triangles(tri))
    sizehint!(triangle_to_circumcenter, num_triangles(tri))
    cur_ghost_idx = I(0)
    cocircular_circumcenters = Set{I}()
    for V in each_triangle(tri)
        V = sort_triangle(V)
        haskey(triangle_to_circumcenter, V) && continue # it was given the circumcenter of a neighbour with the same circumcircle
        if !is_ghost_triangle(V)
            u, v, w = triangle_vertices(V)
            cx, cy = get_voronoi_vertex(tri, V)
            if any(isinf, (cx, cy)) && INF_WARN[]
                @warn "The triangle $((u, v, w)) has a degenerate circumcenter, $((cx, cy)). You may encounter issues with this tessellation. You can disable this warning using toggle_inf_warn!()."
            end
            push_point!(polygon_points, cx, cy)
            idx = I(num_points(polygon_points))
            circumcenter_to_triangle[idx] = V
            triangle_to_circumcenter[V] = idx
            assign_cocircular_neighbours!(triangle_to_circumcenter, cocircular_circumcenters, tri, V, idx, (cx, cy), predicates)
        else
            circumcenter_to_triangle[I(𝒢) - cur_ghost_idx] = V
            triangle_to_circumcenter[V] = I(𝒢) - cur_ghost_idx
            cur_ghost_idx += I(1)
        end
    end
    polygons = Dict{I, Vector{I}}()
    sizehint!(polygons, num_solid_vertices(tri))
    unbounded_polygons = Set{I}()
    sizehint!(unbounded_polygons, num_ghost_edges(tri))
    generators = Dict{I, P}()
    sizehint!(generators, num_solid_vertices(tri))
    for i in each_solid_vertex(tri)
        generators[i] = get_point(tri, i)
    end
    E = edge_type(tri)
    adj = Adjacent{I, E}()
    boundary_polygons = Set{I}()
    return VoronoiTessellation{Tr, P, I, T, typeof(cocircular_circumcenters), E}(tri, generators, polygon_points, polygons, circumcenter_to_triangle, triangle_to_circumcenter, unbounded_polygons, cocircular_circumcenters, adj, boundary_polygons)
end

"""
    get_voronoi_vertex(tri::Triangulation, V) -> NTuple{2, Number}

Returns the vertex of the Voronoi tessellation of `tri` that is dual to the solid triangle `V`, i.e. the circumcenter of `V`, or
its orthocenter if `tri` is weighted. It has the same number type as the points of `tri`; the orthocenter has the type of the
weights if they are wider than the points, e.g. `Float64` weights with `Float32` points, so it is converted.
"""
function get_voronoi_vertex(tri::Triangulation, V)
    F = number_type(tri)
    if !is_weighted(tri)
        u, v, w = triangle_vertices(V)
        p, q, r = get_point(tri, u, v, w)
        A = triangle_area(p, q, r)
        cx, cy = triangle_circumcenter(p, q, r, A)
    else
        cx, cy = triangle_orthocenter(tri, V)
    end
    return F(cx), F(cy)
end

"""
    assign_cocircular_neighbours!(triangle_to_circumcenter, cocircular_circumcenters, tri::Triangulation, V, idx, c, predicates::AbstractPredicateKernel=AdaptiveKernel())

Gives the vertex `idx`, at the coordinates `c`, of the Voronoi tessellation of `tri` that is dual to the triangle `V` to every
triangle of `tri` that is connected to `V` through neighbouring triangles with the same dual vertex, and adds `idx` to
`cocircular_circumcenters` if there are any. A neighbour has the same dual vertex if its opposite vertex is exactly on the
circumcircle of the triangle it neighbours (or, for a weighted triangulation, on its witness plane), as decided by
[`point_position_relative_to_circumcircle`](@ref), or if its computed vertex has the same coordinates `c`.

The exact test is needed since the circumcenters of triangles with the same circumcircle can be rounded differently, e.g.
`(0.15, 0.85)` and `(0.15000000000000002, 0.85)` for four points on a circle. Comparing coordinates also merges neighbours
whose circumcenters are only equal after rounding, so that the tessellation does not get an edge of length zero. Only
neighbours are merged: triangles that do not neighbour each other can have circumcenters that are rounded to the same
coordinates while being different, and merging those used to give a polygon the same vertex twice.
"""
function assign_cocircular_neighbours!(triangle_to_circumcenter, cocircular_circumcenters, tri::Triangulation, V, idx, c, predicates::AbstractPredicateKernel = AdaptiveKernel())
    stack = [V]
    while !isempty(stack)
        W = pop!(stack)
        for (a, b) in triangle_edges(W)
            ℓ = get_adjacent(tri, b, a)
            (!edge_exists(ℓ) || is_ghost_vertex(ℓ)) && continue
            N = sort_triangle(construct_triangle(triangle_type(tri), b, a, ℓ))
            haskey(triangle_to_circumcenter, N) && continue
            if is_on(point_position_relative_to_circumcircle(predicates, tri, W, ℓ)) || get_voronoi_vertex(tri, N) == c
                triangle_to_circumcenter[N] = idx
                push!(cocircular_circumcenters, idx)
                push!(stack, N)
            end
        end
    end
    return triangle_to_circumcenter
end

"""
    prepare_add_voronoi_polygon(vorn::VoronoiTessellation, i) -> (Vector, Vector)

Prepare to add a Voronoi polygon for the vertex `i` to the Voronoi tessellation `vorn`.

# Arguments 
- `vorn`: The [`VoronoiTessellation`](@ref).
- `i`: The vertex. 

# Outputs
- `S`: The surrounding polygon of `i`. See [`get_surrounding_polygon`](@ref).
- `B`: The buffer for the circumcenters. This is an empty `Vector{I}`, where `I = integer_type(tri)`.
"""
function prepare_add_voronoi_polygon(vorn::VoronoiTessellation, i)
    I = integer_type(vorn)
    S = get_surrounding_polygon(vorn, i)
    B = I[]
    sizehint!(B, length(S))
    return S, B
end

"""
    get_next_triangle_for_voronoi_polygon(vorn::VoronoiTessellation, i, k, S, m) -> (Vertex, Vertex)

Get the next triangle for the Voronoi polygon for the point `i` in the `VoronoiTessellation`. 

# Arguments 
- `vorn`: The [`VoronoiTessellation`](@ref).
- `i`: The polygon index.
- `k`: The vertex to add.
- `S`: The surrounding polygon of `i`. See [`get_surrounding_polygon`](@ref).
- `m`: The index of the next vertex in `S`.

# Outputs
- `ci`: The index for the circumcenter of the next triangle.
- `k`: The next vertex in `S` after the input `k`.
"""
function get_next_triangle_for_voronoi_polygon(vorn::VoronoiTessellation, i, k, S, m)
    T = triangle_type(vorn)
    j = k
    k = S[m]
    V = sort_triangle(construct_triangle(T, i, j, k))
    ci = get_triangle_to_circumcenter(vorn, V)
    return ci, k
end

"""
    connect_circumcenters!(B, ci)

Add the circumcenter index `ci` to the array `B`.
"""
function connect_circumcenters!(B, ci)
    push!(B, ci)
    return B
end

"""
    add_edge_to_voronoi_polygon!(B, vorn::VoronoiTessellation, i, k, S, m, encountered_duplicate_circumcenter) -> (Vertex, Bool, Vertex)

Add the next edge to the Voronoi polygon for the point `i` in the `VoronoiTessellation` `vorn`.

# Arguments
- `B`: The vector of circumcenters defining the polygon. 
- `vorn`: The [`VoronoiTessellation`](@ref).
- `i`: The polygon index.
- `k`: The vertex to add.
- `S`: The surrounding polygon of `i`. See [`get_surrounding_polygon`](@ref).
- `m`: The index of the next vertex in `S`.
- `encountered_duplicate_circumcenter`: Whether or not a duplicate circumcenter has been encountered.

# Outputs
- `ci`: The index for the circumcenter of the triangle considered.
- `encountered_duplicate_circumcenter`: Whether or not a duplicate circumcenter has been encountered.
- `k`: The next vertex in `S` after the input `k`.
"""
function add_edge_to_voronoi_polygon!(B, vorn::VoronoiTessellation, i, k, S, m, encountered_duplicate_circumcenter)
    ci, k = get_next_triangle_for_voronoi_polygon(vorn, i, k, S, m)
    is_ghost_vertex(ci) && add_unbounded_polygon!(vorn, i)
    (encountered_duplicate_circumcenter || ci ∈ get_cocircular_circumcenters(vorn)) && (encountered_duplicate_circumcenter = true)
    connect_circumcenters!(B, ci)
    return ci, encountered_duplicate_circumcenter, k
end

"""
    close_voronoi_polygon!(vorn::VoronoiTessellation, B, i, encountered_duplicate_circumcenter, prev_ci)

Close the Voronoi polygon for the point `i` in the `VoronoiTessellation` `vorn`.

# Arguments 
- `vorn`: The [`VoronoiTessellation`](@ref).
- `B`: The vector of circumcenters defining the polygon.
- `i`: The polygon index.
- `encountered_duplicate_circumcenter`: Whether or not a duplicate circumcenter has been encountered.
- `prev_ci`: The previous circumcenter index.

# Outputs 
There are no outputs, as `vorn` and `B` are modified in-place.
"""
function close_voronoi_polygon!(vorn::VoronoiTessellation, B, i, encountered_duplicate_circumcenter, prev_ci)
    encountered_duplicate_circumcenter && unique!(B)
    connect_circumcenters!(B, B[begin])
    add_adjacent!(vorn, prev_ci, B[begin], i)
    add_polygon!(vorn, B, i)
    return nothing
end

"""
    add_voronoi_polygon!(vorn::VoronoiTessellation, i) -> Vector

Add the Voronoi polygon for the point `i` to the `VoronoiTessellation` `vorn`.

# Arguments
- `vorn`: The [`VoronoiTessellation`](@ref).
- `i`: The polygon index.

# Outputs
- `B`: The vector of circumcenters defining the polygon. This is a circular vector, i.e. `B[begin] == B[end]`.
"""
function add_voronoi_polygon!(vorn::VoronoiTessellation, i)
    S, B = prepare_add_voronoi_polygon(vorn, i)
    m = firstindex(S) + 1
    k = S[begin]
    encountered_duplicate_circumcenter = false
    prev_ci, encountered_duplicate_circumcenter, k = add_edge_to_voronoi_polygon!(B, vorn, i, k, S, m, encountered_duplicate_circumcenter)
    for m in (firstindex(S) + 2):lastindex(S)
        ci, encountered_duplicate_circumcenter, k = add_edge_to_voronoi_polygon!(B, vorn, i, k, S, m, encountered_duplicate_circumcenter)
        add_adjacent!(vorn, prev_ci, ci, i)
        prev_ci = ci
    end
    close_voronoi_polygon!(vorn, B, i, encountered_duplicate_circumcenter, prev_ci)
    return B
end
