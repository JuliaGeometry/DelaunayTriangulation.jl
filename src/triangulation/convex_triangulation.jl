"""
    triangulate_convex(points, S;
        IntegerType::Type{I}=Int64,
        EdgeType::Type{E}=NTuple{2,IntegerType},
        TriangleType::Type{V}=NTuple{3,IntegerType},
        EdgesType::Type{Es}=Set{EdgeType},
        TrianglesType::Type{Ts}=Set{TriangleType},
        rng::AbstractRNG=Random.default_rng(),
        add_ghost_triangles=false,
        add_convex_hull=false,
        compute_centers=false,
        delete_empty_features=false) where {I,E,V,Es,Ts}

Given some `points` and a counter-clockwise list of indices `S` for points in `points` 
that define a convex polygon, triangulates it with Chew's algorithm.


# Arguments 
- `points::P`: The set of points.
- `S`: A counter-clockwise list of vertices defining a convex polygon from the corresponding points in `points`.

# Keyword Arguments 
- `IntegerType::Type{I}=Int64`: The integer type to use for indexing. 
- `EdgeType::Type{E}=NTuple{2,IntegerType}`: The type to use for representing edges. 
- `TriangleType::Type{V}=NTuple{3,IntegerType}`: The type to use for representing triangles. 
- `EdgesType::Type{Es}=Set{EdgeType}`: The type to use for representing collections of edges. 
- `TrianglesType::Type{Ts}=Set{TriangleType}`: The type to use for representing collections of triangles. 
- `representative_point_list = get_empty_representative_points(IntegerType, number_type(points))`: The representative point list to use.    
- `rng::AbstractRNG=Random.default_rng()`: The RNG to use.
- `add_ghost_triangles=true`: Whether to add the ghost triangles at the end of the triangulation. 
- `add_convex_hull=true`: Whether to populate the convex hull field of `tri` with `S` at the end of the triangulation.
- `compute_centers=true`: Whether to recompute the `RepresentativePointList` at the end of the triangulation.
- `delete_empty_features=true`: Whether to delete any empty neighbourhoods and adjacencies at the end of the triangulation. 

# Outputs 
Returns a [`Triangulation`](@ref) of the convex polygon.
"""
function triangulate_convex(points, S;
    IntegerType::Type{I}=Int64,
    EdgeType::Type{E}=NTuple{2,IntegerType},
    TriangleType::Type{V}=NTuple{3,IntegerType},
    EdgesType::Type{Es}=Set{EdgeType},
    TrianglesType::Type{Ts}=Set{TriangleType},
    representative_point_list = get_empty_representative_points(IntegerType, number_type(points)),
    rng::AbstractRNG=Random.default_rng(),
    add_ghost_triangles=true,
    add_convex_hull=true,
    compute_centers=true,
    delete_empty_features=true) where {I,E,V,Es,Ts}
    tri = Triangulation(points;
        IntegerType,
        EdgeType,
        TriangleType,
        EdgesType,
        TrianglesType,
        representative_point_list)
    triangulate_convex!(tri, S;
        rng,
        add_ghost_triangles,
        add_convex_hull,
        compute_centers,
        delete_empty_features)
    return tri
end

# note that the keyword arguments here do not match those above
function triangulate_convex!(tri::Triangulation, S;
    rng,
    add_ghost_triangles=false,
    add_convex_hull=false,
    compute_centers=false,
    delete_empty_features=false,
    update_ghost_edges=false)
    next, prev, k, S, shuffled_indices = prepare_convex_triangulation_vectors(S)
    reset_convex_triangulation_vectors!(next, prev, shuffled_indices, k, rng)
    delete_vertices_in_random_order!(next, prev, shuffled_indices, tri, S, k, rng)
    u, v, w = index_shuffled_linked_list(S, next, prev, shuffled_indices, 1)
    add_triangle!(tri, u, v, w; protect_boundary=!update_ghost_edges, update_ghost_edges)
    set_S = Set(S)
    for i in 4:k
        u, v, w = index_shuffled_linked_list(S, next, prev, shuffled_indices, i)
        add_point_convex_triangulation!(tri, u, v, w, set_S, update_ghost_edges)
    end
    convex_triangulation_post_processing!(tri,
        S,
        add_ghost_triangles,
        add_convex_hull,
        compute_centers,
        delete_empty_features)
    return nothing
end

function prepare_convex_triangulation_vectors(S::AbstractArray{I}) where {I}
    Base.require_one_based_indexing(S)
    if is_circular(S)
        S = @views S[begin:(end-1)]
    end
    k = length(S)
    next = zeros(I, k)
    prev = zeros(I, k)
    shuffled_indices = collect(eachindex(S))
    return next, prev, k, S, shuffled_indices
end

function reset_convex_triangulation_vectors!(next, prev, shuffled_indices, k, rng)
    for i in eachindex(next, prev)
        next[i] = mod1(i + 1, k)
        prev[i] = mod1(i - 1, k)
    end
    shuffle!(rng, shuffled_indices)
    return nothing
end

function index_shuffled_linked_list(S, next, prev, shuffled_indices, i)
    πᵢ = shuffled_indices[i]
    u = S[πᵢ]
    v = S[next[πᵢ]]
    w = S[prev[πᵢ]]
    return u, v, w
end

function delete_vertices_in_random_order!(next, prev, shuffled_indices, tri, S, k, rng)
    ## Delete the points in a random order
    for i in k:-1:4
        πᵢ = shuffled_indices[i]
        next[prev[πᵢ]] = next[πᵢ]
        prev[next[πᵢ]] = prev[πᵢ]
    end

    ## Check if the three surviving vertices that survived are collinear 
    # Note that we do not need to check for negative orientation here, 
    # since the polygon is provided in a counter-clockwise order already. 
    u, v, w = index_shuffled_linked_list(S, next, prev, shuffled_indices, 1)
    a, b, c = get_point(tri, u, v, w) # don't special case for boundary indices 
    degenerate_cert = triangle_orientation(a, b, c)
    if !is_positively_oriented(degenerate_cert) #is_degenerate(degenerate_cert)
        reset_convex_triangulation_vectors!(next, prev, shuffled_indices, k, rng)
        delete_vertices_in_random_order!(next, prev, shuffled_indices, tri, S, k, rng)
    end
    return nothing
end

function convex_triangulation_post_processing!(tri,
    S,
    add_ghost_triangles,
    add_convex_hull,
    compute_centers,
    delete_empty_features)
    I = integer_type(tri)
    if add_ghost_triangles
        for i in 1:(length(S)-1)
            u = S[i]
            v = S[i+1]
            add_triangle!(tri, v, u, I(BoundaryIndex))
        end
        u = S[end]
        v = S[begin]
        add_triangle!(tri, v, u, I(BoundaryIndex))
    end
    if add_convex_hull || compute_centers
        idx = get_convex_hull_indices(tri)
        append!(idx, S)
        push!(idx, S[begin])
    end
    compute_centers && compute_representative_points!(tri; use_convex_hull=true)
    delete_empty_features && clear_empty_features!(tri)
    return nothing
end

function add_point_convex_triangulation!(tri::Triangulation, u, v, w, S, update_ghost_edges) 
    x = get_adjacent(tri, w, v)
    # In the statement below, we check x ∈ S so that we do not leave the 
    # polygon. If we are just triangulation a convex polygon in isolation, 
    # this is just checking that edge_exists(tri, w, v). If instead 
    # this is a polygon representing a cavity of a triangulation, this 
    # is just checking that x is a vertex of the cavity. This check should 
    # be O(1), hence why we use Set(S) in triangulate_convex.
    # When we have an outer boundary index, we replace the RepresentativePointList[1] with some other point.
    # To keep this behaviour and not treat this point as if it is inside the triangulation,
    # we avoid is_inside(point_position_relative_to_circumcircle(tri, u, v, w, x)), which special-cases 
    # for outer boundary indices. This'll work for other points too, so we just do that.
    if x ∈ S# && !is_boundary_index(x)
        a, b, c, p = get_point(tri, u, v, w, x) # do this in an awkward way (i.e. without an if statement for defining flag) since we might access get_point with a DefaultAdjacentValue index
        flag = is_inside(point_position_relative_to_circle(a, b, c, p))
    else
        flag = false
    end
    if flag
        # uvw and wvx are not Delaunay
        delete_triangle!(tri, w, v, x; protect_boundary=!update_ghost_edges, update_ghost_edges)
        add_point_convex_triangulation!(tri, u, v, x, S, update_ghost_edges)
        add_point_convex_triangulation!(tri, u, x, w, S, update_ghost_edges)
    else
        # vw is a Delaunay edge
        a, b, c = get_point(tri, u, v, w)
        add_triangle!(tri, u, v, w; protect_boundary=!update_ghost_edges, update_ghost_edges)
        return nothing
    end
end