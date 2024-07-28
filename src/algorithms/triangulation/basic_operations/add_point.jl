"""
    add_point!(tri::Triangulation, new_point; kwargs...) -> Triangle
    add_point!(tri::Triangulation, x, y; kwargs...) -> Triangle 

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `new_point`: The point to be added to the triangulation. The second method uses `(x, y)` to represent the new point instead. If `new_point` is an integer, then the point added is `get_point(tri, new_point)`.

# Keyword Arguments 
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.
- `point_indices=each_solid_vertex(tri)`: The indices of the points to be used in the [`find_triangle`](@ref) algorithm for selecting the initial point.
- `m=default_num_samples(length(point_indices))`: The number of samples (without replacement) to be used in the [`find_triangle`](@ref) algorithm for selecting the initial point.
- `try_points=()`: Additional points to try for selecting the initial point, in addition to the `m` sampled.
- `rng::Random.AbstractRNG=Random.default_rng()`: The random number generator to be used in [`find_triangle`](@ref).
- `initial_search_point=integer_type(tri)(select_initial_point(tri, new_point; point_indices, m, try_points, rng))`: The initial point to be used in [`find_triangle`](@ref).
- `update_representative_point=false`: Whether to update the representative point of the triangulation after adding the new point. 
- `store_event_history=Val(false)`: Whether to store the event history of the triangulation from adding the new point. 
- `event_history=nothing`: The event history of the triangulation from adding the new point. Only updated if `store_event_history` is true, in which case it needs to be an [`InsertionEventHistory`](@ref) object.
- `concavity_protection=false`: Whether to use concavity protection for finding `V` below. See [`concavity_protection_check`](@ref). This is only needed if your triangulation is not convex. 
- `V=find_triangle(tri, get_point(tri, new_point); m=nothing, point_indices=nothing, try_points=nothing, k=initial_search_point, concavity_protection, rng)`: The positively oriented triangle containing the point being added.

!!! warning "Non-convex domains"

    In cases where your triangulation is not convex and `!concavity_protection`, this `V` may not be correct, and you may encounter errors - errors either during `add_point!` or separately when 
    you try to use the triangulation. In such cases, you should set `concavity_protection=true` to ensure that `V` is correct.

- `peek=Val(false)`: Whether the point should actually be added into the triangulation, or just 'peeked' at so that the events that would occur from its addition can be added into `event_history`.

# Outputs 
The triangulation is updated in-place, but we do return 
- `V`: The triangle containing the point being added.

!!! warn "Convex hull"

    In cases where `(x, y)` is outside of the triangulation, it will be added successfully but note that 
    the `convex_hull` field of `tri` will no longer be accurate. You can use [`convex_hull!`](@ref) to fix it.
"""
function add_point!(
        tri::Triangulation, new_point;
        predicates::AbstractPredicateKernel = AdaptiveKernel(),
        point_indices = each_solid_vertex(tri),
        m = default_num_samples(length(point_indices)),
        try_points = (),
        rng::Random.AbstractRNG = Random.default_rng(),
        initial_search_point = integer_type(tri)(select_initial_point(tri, new_point; point_indices, m, try_points, rng)),
        update_representative_point = false,
        store_event_history = Val(false),
        event_history = nothing,
        concavity_protection = false,
        V = find_triangle(
            tri,
            get_point(tri, new_point);
            m = nothing,
            point_indices = nothing,
            try_points = nothing,
            k = initial_search_point,
            concavity_protection,
            predicates,
            rng,
        ),
        peek::P = Val(false),
    ) where {P}
    int_flag = new_point isa Integer
    if !int_flag && !is_true(peek)
        push_point!(tri, new_point)
        new_point = num_points(tri)
    end
    q = get_point(tri, new_point)
    if is_weighted(tri)
        cert = point_position_relative_to_circumcircle(predicates, tri, V, new_point) # redirects to point_position_relative_to_witness_plane
        is_outside(cert) && return V # If the point is submerged, then we don't add it
    end
    flag = point_position_relative_to_triangle(predicates, tri, V, q)
    V = add_point_bowyer_watson_and_process_after_found_triangle!(tri, new_point, V, q, flag, update_representative_point, store_event_history, event_history, peek, predicates)
    return V
end

function add_point!(
        tri::Triangulation, new_point_x, new_point_y;
        predicates::AbstractPredicateKernel = AdaptiveKernel(),
        point_indices = get_vertices(tri),
        m = default_num_samples(length(point_indices)),
        try_points = (),
        rng::Random.AbstractRNG = Random.default_rng(),
        initial_search_point = integer_type(tri)(select_initial_point(tri, (new_point_x, new_point_y); point_indices, m, try_points, rng)),
        update_representative_point = false,
        store_event_history = Val(false),
        event_history = nothing,
        concavity_protection = false,
        V = find_triangle(
            tri,
            (new_point_x, new_point_y);
            m = nothing,
            point_indices = nothing,
            try_points = nothing,
            k = initial_search_point,
            concavity_protection,
            predicates,
            rng,
        ),
        peek::P = Val(false),
    ) where {P}
    !is_true(peek) && push_point!(tri, new_point_x, new_point_y)
    VV = add_point!(
        tri,
        is_true(peek) ? (new_point_x, new_point_y) : num_points(tri);
        predicates,
        point_indices = point_indices,
        m = m,
        try_points = try_points,
        rng = rng,
        initial_search_point = initial_search_point,
        update_representative_point = update_representative_point,
        store_event_history = store_event_history,
        event_history = event_history,
        V,
        peek,
    )
    return VV
end

"""
    add_point!(tri::Triangulation, x, y, w; kwargs...)

Adds the point `(x, y)` into `tri` with weight `w`. This function requires that [`add_weight!`](@ref)
is defined on the weights stored in `tri`. The `kwargs` match those from [`add_point!(tri::Triangulation, ::Any)`](@ref).
"""
function add_point!(tri::Triangulation, x, y, w; kwargs...)
    add_weight!(tri, w)
    return add_point!(tri, x, y; kwargs...)
end
