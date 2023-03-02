function prepare_convex_triangulation_vectors(S::AbstractArray{I}) where {I}
    k = length(S)
    next = zeros(I, k)
    prev = zeros(I, k)
    seen = Set{I}()
    sizehint!(seen, k - 3)
    shuffled_indices = collect(eachindex(S))
    return k, next, prev, seen, shuffled_indices
end

function reset_convex_triangulation_vectors!(next, prev, seen, k)
    empty!(seen)
    sizehint!(seen, k - 3)
    for i in 1:k
        next[i] = mod1(i + 1, k)
        prev[i] = mod1(i - 1, k)
    end
    return nothing
end

function delete_convex_polygon_vertices_in_random_order!(shuffled_indices, next::AbstractArray{I}, prev, seen, _S, points, k, rng) where {I}
    ## Shuffle the points
    shuffle!(rng, shuffled_indices)
    for i in k:-1:4
        next[prev[shuffled_indices[i]]] = next[shuffled_indices[i]]
        prev[next[shuffled_indices[i]]] = prev[shuffled_indices[i]]
    end

    ## Check if the three surviving vertices that survived are collinear 
    # Note that we do not need to check for negative orientation here, 
    # since the polygon is provided in a counter-clockwise order already.
    init = shuffled_indices[begin]
    u, v, w = _S[init], _S[next[init]], _S[prev[init]]
    p, q, r = get_point(points, u, v, w)
    degenerate_cert = triangle_orientation(p, q, r)
    while is_degenerate(degenerate_cert)
        reset_convex_triangulation_vectors!(next, prev, seen, k)
        delete_convex_polygon_vertices_in_random_order!(shuffled_indices, next, prev, seen, _S, points, k, rng)
    end
    return nothing
end

function initialise_convex_triangulation(points::P, S, shuffled_indices, next, prev;
    IntegerType::Type{I}=Int64,
    EdgeType::Type{E}=NTuple{2,IntegerType},
    TriangleType::Type{V}=NTuple{3,IntegerType},
    EdgesType::Type{Es}=Set{EdgeType},
    TrianglesType::Type{Ts}=Set{TriangleType}) where {P,I,E,V,Es,Ts}
    tri = Triangulation(points; IntegerType, EdgeType, TriangleType, EdgesType,
        TrianglesType)
    πᵢ = shuffled_indices[begin]
    u = S[πᵢ]
    v = S[next[πᵢ]]
    w = S[prev[πᵢ]] # Orientation has already been checked in delete_convex_polygon_vertices_in_random_order!
    initial_triangle = construct_triangle(V, u, v, w)
    add_triangle!(tri, initial_triangle; update_ghost_edges=false)
    return tri
end

function triangulate_convex(points, S;
    IntegerType::Type{I}=Int64,
    EdgeType::Type{E}=NTuple{2,IntegerType},
    TriangleType::Type{V}=NTuple{3,IntegerType},
    EdgesType::Type{Es}=Set{EdgeType},
    TrianglesType::Type{Ts}=Set{TriangleType},
    rng::AbstractRNG=Random.default_rng(),
    add_ghost_triangles=false) where {I,E,V,Es,Ts}
    Base.require_one_based_indexing(S)
    if is_circular(S)
        S = @views S[begin:(end-1)]
    end
    k, next, prev, seen, shuffled_indices = prepare_convex_triangulation_vectors(S)
    reset_convex_triangulation_vectors!(next, prev, seen, k) # initialise
    delete_convex_polygon_vertices_in_random_order!(shuffled_indices, next, prev, seen, S, points, k, rng)

    # Now initialise the triangulation with the three surviving vertices 
    tri = initialise_convex_triangulation(points, S, shuffled_indices, next, prev;
        IntegerType,
        EdgeType,
        TriangleType,
        EdgesType,
        TrianglesType)

    ## Now add in the triangles
    set_S = Set(S) # for O(1) lookup with ∈
    for i in 4:k
        πᵢ = shuffled_indices[i]
        u = S[πᵢ]
        v = S[next[πᵢ]]
        w = S[prev[πᵢ]]
        add_point_convex_triangulation!(tri, u, v, w, set_S)
    end
    add_ghost_triangles && add_ghost_triangles!(tri)
    return tri
end

function add_point_convex_triangulation!(tri::Triangulation, u, v, w, S)
    x = get_adjacent(tri, w, v)
    @show x
    # In the statement below, we check x ∈ S so that we do not leave the 
    # polygon. If we are just triangulation a convex polygon in isolation, 
    # this is just checking that edge_exists(tri, w, v). If instead 
    # this is a polygon representing a cavity of a triangulation, this 
    # is just checking that x is a vertex of the cavity. This check should 
    # be O(1), hence why we use Set(S) in triangulate_convex.
    if x ∈ S && is_inside(point_position_relative_to_circumcircle(tri, u, v, w, x))
        # uvw and wvx are not Delaunay
        @show tri.triangles, w, v, x
        delete_triangle!(tri, w, v, x)
        add_point_convex_triangulation!(tri, u, v, x, S)
        add_point_convex_triangulation!(tri, u, x, w, S)
        return nothing
    else
        # vw is a Delaunay edge
        add_triangle!(tri, u, v, w; protect_boundary=true)
        return nothing
    end
end