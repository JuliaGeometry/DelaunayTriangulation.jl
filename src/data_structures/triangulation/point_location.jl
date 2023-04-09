"""
    brute_force_search(tri::Triangulation, q)

Returns the triangle in `tri` containing `q` using brute force search.
"""
function brute_force_search(tri::Triangulation, q)
    T = get_triangles(tri)
    points = get_points(tri)
    represetative_point_list = get_representative_point_list(tri)
    boundary_map = get_boundary_map(tri)
    return brute_force_search(T, q, points, represetative_point_list, boundary_map)
end

"""
    jump_and_march(tri::Triangulation, q;
        point_indices=each_point_index(tri),
        m=default_num_samples(length(point_indices)),
        try_points=(),
        k=select_initial_point(get_points(tri), q; m, point_indices, try_points),
        check_existence::C=Val(has_multiple_segments(tri)),
        store_history::F=Val(false),
        history=nothing,
        rng::AbstractRNG=Random.default_rng()) where {C,F}

Returns the triangle containing `q` using the jump-and-march algorithm.

# Arguments 
- `tri::Triangulation`: The triangulation.
- `q`: The query point.

# Keyword Arguments
- `point_indices=each_point_index(tri)`: The indices of the points in the triangulation.
- `m=default_num_samples(length(point_indices))`: The number of samples to use when sampling the point to start the algorithm at.
- `try_points=()`: Additional points to try when determining which point to start at.
- `k=select_initial_point(get_points(tri), q; m, point_indices, try_points)`: The index of the point to start the algorithm at. 
- `check_existence::C=Val(has_multiple_segments(tri))`: This is used when we want to check the existence of certain ghost triangles. See [`get_adjacent`}(@ref).
- `store_history::F=Val(false)`: Whether to record the history of the algorithm. See also [`PointLocationHistory`](@ref).
- `history=nothing`: The object to store the history into, if `is_true(store_history)`.
- `rng::AbstractRNG=Random.default_rng()`: The random number generator to use.

# Outputs 
Returns `V`, the triangle in `tri` containing `q`.
"""
function jump_and_march(tri::Triangulation, q;
    point_indices=each_solid_vertex(tri),
    m=default_num_samples(num_vertices(point_indices)),
    try_points=(),
    k=select_initial_point(get_points(tri), q; m, point_indices, try_points),
    check_existence::C=Val(has_multiple_segments(tri)),
    store_history::F=Val(false),
    history=nothing,
    rng::AbstractRNG=Random.default_rng()) where {C,F}
    points = get_points(tri)
    adjacent = get_adjacent(tri)
    adjacent2vertex = get_adjacent2vertex(tri)
    graph = get_graph(tri)
    boundary_index_ranges = get_boundary_index_ranges(tri)
    representative_point_list = get_representative_point_list(tri)
    boundary_map = get_boundary_map(tri)
    return jump_and_march(points, adjacent, adjacent2vertex, graph, boundary_index_ranges, representative_point_list, boundary_map, q;
        m,
        point_indices,
        try_points,
        k,
        check_existence,
        store_history,
        history,
        rng)
end