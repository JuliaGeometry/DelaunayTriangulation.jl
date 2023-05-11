"""
    add_point!(tri::Triangulation, new_point[, new_point_y];
        point_indices=get_vertices(tri),
        m=default_num_samples(length(point_indices)),
        try_points=(),
        rng::AbstractRNG=Random.default_rng(),
        initial_search_point=integer_type(tri)(select_initial_point(get_points(tri),new_point;point_indices,m,try_points,rng)),
        update_representative_point=false,
        store_event_history = Val(false),
        event_history = nothing,
        exterior_curve_index=1,
        V=jump_and_march(
            tri,
            new_point isa Integer ? get_point(tri, new_point) : new_point;
            m=nothing,
            point_indices=nothing,
            try_points=nothing,
            k=initial_search_point,
            rng,
            check_existence=Val(has_multiple_segments(tri)),
            exterior_curve_index
        ),
        peek = Val(false),
        )

Adds the point `new_point` to the triangulation `tri`.

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `new_point[, new_point_y]`: The point to add. This `new_point` can be an integer, in which case `get_point(tri, new_point)` is added. If `new_point` is just a set of coordinates, we add that into `tri` via [`push_point!`](@ref) and then add that index into `tri`. Lastly, if we provide `(new_point, new_point_y)`, then the point is treated as this `Tuple` and inserted.

# Keyword Arguments 
- `point_indices=each_solid_vertex(tri)`: The indices of the non-ghost points in the triangulation.
- `m=default_num_samples(length(point_indices))`: The number of points to sample from `point_indices` to use as the initial search point.
- `try_points=()`: A list of points to try as the initial search point in addition to those sampled.
- `rng::AbstractRNG=Random.default_rng()`: The random number generator to use.
- `initial_search_point=integer_type(tri)(select_initial_point(get_points(tri),new_point;point_indices,m,try_points,rng))`: The initial search point to use. If this is not provided, then we use [`select_initial_point`](@ref) to select one.
- `update_representative_point=false`: Whether to update the representative point list after adding the new point.
- `store_event_history = Val(false)`: Whether to store the event history. See [`InsertionEventHistory`](@ref).
- `event_history = nothing`: The event history to store the events in. See [`InsertionEventHistory`](@ref). Only needed if `is_true(store_event_history)`. This object is not returned, instead we just mutate it inplace.
- `exterior_curve_index=1`: The curve (or curves) corresponding to the outermost boundary.
- `V=jump_and_march(tri, new_point isa Integer ? get_point(tri, new_point) : new_point; m=nothing, point_indices=nothing, try_points=nothing, k=initial_search_point, rng, check_existence=Val(has_multiple_segments(tri)), exterior_curve_index=exterior_curve_index)`: The triangle that `q` is in.
- `peek=Val(false)`: If `is_true(peek)`, then we don't actually add the point, but all operations that update the history will be run. (So you should only really want this if you are using `event_history`.)

# Outputs 
The triangulation is updated in-place with the new point, but we also return the triangle `V` containing `new_point`.
"""
function add_point!(tri::Triangulation, new_point;
    point_indices=each_solid_vertex(tri),
    m=default_num_samples(length(point_indices)),
    try_points=(),
    rng::AbstractRNG=Random.default_rng(),
    initial_search_point=integer_type(tri)(select_initial_point(get_points(tri), new_point; point_indices, m, try_points, rng)),
    update_representative_point=false,
    store_event_history=Val(false),
    event_history=nothing,
    exterior_curve_index=1,
    V=jump_and_march(
        tri,
        new_point isa Integer ? get_point(tri, new_point) : new_point;
        m=nothing,
        point_indices=nothing,
        try_points=nothing,
        k=initial_search_point,
        rng,
        check_existence=Val(has_multiple_segments(tri)),
        exterior_curve_index
    ),
    peek=Val(false))
    int_flag = new_point isa Integer
    if !int_flag
        push_point!(tri, new_point)
        new_point = num_points(tri)
    end
    q = get_point(tri, new_point)
    flag = point_position_relative_to_triangle(tri, V, q)
    V = add_point_bowyer_watson_and_process_after_found_triangle!(tri, new_point, V, q, flag, update_representative_point, store_event_history, event_history, peek)
    if !int_flag && is_true(peek)
        pop_point!(tri)
    end 
    return V
end

function add_point!(tri::Triangulation, new_point_x, new_point_y;
    point_indices=get_vertices(tri),
    m=default_num_samples(length(point_indices)),
    try_points=(),
    rng::AbstractRNG=Random.default_rng(),
    initial_search_point=integer_type(tri)(select_initial_point(get_points(tri), (new_point_x, new_point_y); point_indices, m, try_points, rng)),
    update_representative_point=false,
    store_event_history=Val(false),
    event_history=nothing,
    exterior_curve_index=1,
    V=jump_and_march(
        tri,
        (new_point_x, new_point_y);
        m=nothing,
        point_indices=nothing,
        try_points=nothing,
        k=initial_search_point,
        rng,
        check_existence=Val(has_multiple_segments(tri)),
        exterior_curve_index
    ),
    peek=Val(false))
    push_point!(tri, new_point_x, new_point_y)
    VV = add_point!(
        tri,
        num_points(tri);
        point_indices=point_indices,
        m=m,
        try_points=try_points,
        rng=rng,
        initial_search_point=initial_search_point,
        update_representative_point=update_representative_point,
        store_event_history=store_event_history,
        event_history=event_history,
        exterior_curve_index=exterior_curve_index,
        V,
        peek)
    is_true(peek) && pop_point!(tri)
    return VV
end