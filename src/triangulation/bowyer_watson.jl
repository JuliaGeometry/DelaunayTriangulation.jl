@inline function get_point_order(points, randomise, skip_points, IntegerType::Type{I},
    rng::AbstractRNG=Random.default_rng()) where {I}
    point_indices = convert(Vector{I}, collect(each_point_index(points)))
    randomise && shuffle!(rng, point_indices)
    setdiff!(point_indices, skip_points)
    return point_indices
end

@inline function get_initial_triangle(::Type{V}, ::Type{I}, point_order, points) where {V,I}
    local initial_triangle
    i, j, k = point_order[begin], point_order[begin+1], point_order[begin+2]
    initial_triangle = construct_positively_oriented_triangle(V, i, j, k, points)
    i, j, k = indices(initial_triangle)
    p, q, r = get_point(points, i, j, k)
    degenerate_cert = triangle_orientation(p, q, r)
    itr = 0
    if length(point_order) > 3 # Do not get stuck in an infinite loop if there are just three points, the three of them being collinear
        while is_degenerate(degenerate_cert) && itr ≤ length(point_order) # If all points are collinear, this loop could go forever
            @static if VERSION ≥ v"1.8.1"
                circshift!(point_order, -1)::Vector{I}
            else
                _point_order = circshift(point_order, -1)
                copyto!(point_order, _point_order)
            end
            i, j, k = point_order[begin], point_order[begin+1], point_order[begin+2]
            initial_triangle = construct_positively_oriented_triangle(V, i, j, k, points)
            i, j, k = indices(initial_triangle)
            p, q, r = get_point(points, i, j, k)
            degenerate_cert = triangle_orientation(p, q, r)
            itr += 1
        end
    end
    return initial_triangle
end
@inline function get_initial_triangle(tri::Triangulation, point_order)
    V = triangle_type(tri)
    I = integer_type(tri)
    points = get_points(tri)
    return get_initial_triangle(V, I, point_order, points)
end

function initialise_bowyer_watson(points::P;
    IntegerType::Type{I}=Int,
    EdgeType::Type{E}=NTuple{2,IntegerType},
    TriangleType::Type{V}=NTuple{3,IntegerType},
    EdgesType::Type{Es}=Set{EdgeType},
    TrianglesType::Type{Ts}=Set{TriangleType},
    randomise=true,
    skip_points=(),
    rng::AbstractRNG=Random.default_rng(),
    point_order=get_point_order(points, randomise, skip_points, IntegerType, rng)) where {P,I,E,V,Es,Ts}
    tri = Triangulation(points; IntegerType, EdgeType, TriangleType, EdgesType, TrianglesType)
    initial_triangle = get_initial_triangle(V, I, point_order, points) # point_order could get mutated here
    add_triangle!(tri, initial_triangle; update_ghost_edges=true)
    new_representative_point!(tri, I(1))
    u, v, w = indices(initial_triangle)
    p, q, r = get_point(tri, u, v, w)
    for pt in (p, q, r)
        update_centroid_after_addition!(tri, I(1), pt)
    end
    return tri::Triangulation{P,Ts,I,E,Es,Vector{I},Dict{E,Tuple{Vector{I},I}},OrderedDict{I,Vector{I}},
        OrderedDict{I,UnitRange{I}},Dict{I,RepresentativeCoordinates{I,number_type(P)}}}
end

function triangulate_bowyer_watson(points::P;
    IntegerType::Type{I}=Int,
    EdgeType::Type{E}=NTuple{2,IntegerType},
    TriangleType::Type{V}=NTuple{3,IntegerType},
    EdgesType::Type{Es}=Set{EdgeType},
    TrianglesType::Type{Ts}=Set{TriangleType},
    randomise=true,
    delete_ghosts=true,
    delete_empty_features=true,
    try_last_inserted_point=true,
    skip_points=(),
    num_sample_rule::M=default_num_samples,
    rng::AbstractRNG=Random.default_rng(),
    point_order=get_point_order(points, randomise, skip_points, IntegerType, rng),
    recompute_representative_point=true) where {P,I,E,V,
    Es,Ts,M}
    tri = initialise_bowyer_watson(points;
        IntegerType=I,
        EdgeType=E,
        TriangleType=V,
        EdgesType=Es,
        TrianglesType=Ts,
        randomise,
        skip_points,
        rng,
        point_order)::Triangulation{
        P,Ts,I,E,Es,Vector{I},Dict{E,Tuple{Vector{I},I}},
        OrderedDict{I,Vector{I}},OrderedDict{I,UnitRange{I}},Dict{I,RepresentativeCoordinates{I,number_type(P)}}
    }
    _triangulate_bowyer_watson!(tri, point_order, num_sample_rule, delete_ghosts,
        delete_empty_features, try_last_inserted_point, rng,
        recompute_representative_point)
    return tri
end

function _triangulate_bowyer_watson!(tri::Triangulation, point_order,
    num_sample_rule::M=n -> default_num_samples(n),
    delete_ghosts=true, delete_empty_features=true,
    try_last_inserted_point=true,
    rng::AbstractRNG=Random.default_rng(),
    recompute_representative_point=true) where {M}
    remaining_points = @view point_order[(begin+3):end]
    for (num_points, new_point) in enumerate(remaining_points)
        num_currently_inserted = num_points + 3 - 1     # + 3 for the points already inserted
        last_inserted_point_index = point_order[num_currently_inserted]
        currently_inserted_points = @view point_order[begin:num_currently_inserted]
        m = num_sample_rule(num_currently_inserted)
        if try_last_inserted_point
            initial_search_point = select_initial_point(
                tri,
                new_point;
                m,
                point_indices=currently_inserted_points,
                try_points=last_inserted_point_index,
                rng
            )
        else
            initial_search_point = select_initial_point(
                tri,
                new_point;
                m,
                point_indices=currently_inserted_points,
                rng
            )
        end
        add_point_bowyer_watson!(tri, new_point, initial_search_point, rng)
    end
    convex_hull!(tri; reconstruct=false)
    recompute_representative_point &&
        compute_representative_points!(tri; use_convex_hull=true)
    delete_ghosts && delete_ghost_triangles!(tri)
    delete_empty_features && clear_empty_features!(tri)
    return nothing
end

function add_point_bowyer_watson_and_process_after_found_triangle!(
    tri::Triangulation,
    new_point,
    V,
    q,
    flag,
    update_representative_point=true,
    store_event_history=Val(false),
    event_history=nothing,
    peek::P=Val(false)) where {P}
    I = integer_type(tri)
    if !is_true(peek)
        new_point = I(new_point)
    end
    add_point_bowyer_watson_after_found_triangle!(tri, new_point, V, q, flag, update_representative_point, store_event_history, event_history, peek)
    add_point_bowyer_watson_onto_constrained_segment!(tri, new_point, V, q, flag, update_representative_point, store_event_history, event_history, peek)
    return V
end

function add_point_bowyer_watson!(
    tri::Triangulation,
    new_point,
    initial_search_point::I,
    rng::AbstractRNG=Random.default_rng(),
    update_representative_point=true,
    store_event_history=Val(false),
    event_history=nothing,
    exterior_curve_index=1,
    peek::P=Val(false)) where {I,P}
    if !is_true(peek)
        new_point = I(new_point)
    end
    q = get_point(tri, new_point)
    V = jump_and_march(tri, q; m=nothing, point_indices=nothing, try_points=nothing, k=initial_search_point, rng, exterior_curve_index)
    flag = point_position_relative_to_triangle(tri, V, q)
    add_point_bowyer_watson_and_process_after_found_triangle!(tri, new_point, V, q, flag, update_representative_point, store_event_history, event_history, peek)
    return V
end

function add_point_bowyer_watson_onto_constrained_segment!(
    tri,
    new_point,
    V,
    q,
    flag,
    update_representative_point=true,
    store_event_history=Val(false),
    event_history=nothing,
    peek::P=Val(false)) where {P}
    _new_point = is_true(peek) ? num_points(tri) + 1 : new_point # If we are peeking, then we need to use the number of points in the triangulation as the index for the new point since we don't actually insert the point
    if is_on(flag) && is_constrained(tri)
        # If the point we are adding appears on a segment, then we perform the depth-first search 
        # on each side the segment. We also need to update the 
        u, v = find_edge(tri, V, new_point)
        if contains_constrained_edge(tri, u, v)
            # First, let us search on the other side of the segment, provide we are not on a boundary edge. 
            w = get_adjacent(tri, v, u)
            T = triangle_type(tri)
            V = construct_triangle(T, v, u, w)
            add_point_bowyer_watson_dig_cavities!(tri, new_point, V, q, flag, update_representative_point, store_event_history, event_history, peek)
            # Now, we need to replace the segment by this new segment.
            E = edge_type(tri)
            constrained_edges = get_constrained_edges(tri)
            all_constrained_edges = get_all_constrained_edges(tri)
            if !is_true(peek)
                for edges in (constrained_edges, all_constrained_edges)
                    delete_edge!(edges, construct_edge(E, u, v))
                    delete_edge!(edges, construct_edge(E, v, u))
                    add_edge!(edges, construct_edge(E, u, new_point))
                    add_edge!(edges, construct_edge(E, new_point, v))
                end
            end
            if is_true(store_event_history)
                delete_edge!(event_history, construct_edge(E, u, v))
                add_edge!(event_history, construct_edge(E, u, _new_point))
                add_edge!(event_history, construct_edge(E, _new_point, v))
            end
        end
        if contains_boundary_edge(tri, u, v)
            !is_true(peek) && split_boundary_edge!(tri, u, v, new_point)
            is_true(store_event_history) && split_boundary_edge!(event_history, u, v, _new_point)
        elseif contains_boundary_edge(tri, v, u)
            !is_true(peek) && split_boundary_edge!(tri, v, u, new_point)
            is_true(store_event_history) && split_boundary_edge!(event_history, v, u, _new_point)
        end
    end
end

function add_point_bowyer_watson_after_found_triangle!(
    tri,
    new_point,
    V,
    q,
    flag,
    update_representative_point=true,
    store_event_history=Val(false),
    event_history=nothing,
    peek::P=Val(false)) where {P}
    if is_ghost_triangle(V) && is_constrained(tri)
        # When we have a constrained boundary edge, we don't want to walk into its 
        # interior. So let's just check this case now.
        V = rotate_ghost_triangle_to_standard_form(V)
        u, v, w = indices(V)
        if !contains_boundary_edge(tri, v, u)
            add_point_bowyer_watson_dig_cavities!(tri, new_point, V, q, flag, update_representative_point, store_event_history, event_history, peek)
        end
    else
        add_point_bowyer_watson_dig_cavities!(tri, new_point, V, q, flag, update_representative_point, store_event_history, event_history, peek)
    end
end

function add_point_bowyer_watson_dig_cavities!(
    tri::Triangulation{P,Ts,I},
    new_point::N,
    V,
    q,
    flag,
    update_representative_point=true,
    store_event_history=Val(false),
    event_history=nothing,
    peek::F=Val(false)) where {P,Ts,I,N,F}
    _new_point = is_true(peek) ? num_points(tri) + 1 : new_point # If we are peeking, then we need to use the number of points in the triangulation as the index for the new point since we don't actually insert the point
    i, j, k = indices(V)
    ℓ₁ = get_adjacent(tri, j, i)
    ℓ₂ = get_adjacent(tri, k, j)
    ℓ₃ = get_adjacent(tri, i, k)
    !is_true(peek) && delete_triangle!(tri, V; protect_boundary=true, update_ghost_edges=false)
    is_true(store_event_history) && delete_triangle!(event_history, V)
    dig_cavity!(tri, new_point, i, j, ℓ₁, flag, V, store_event_history, event_history, peek)
    dig_cavity!(tri, new_point, j, k, ℓ₂, flag, V, store_event_history, event_history, peek)
    dig_cavity!(tri, new_point, k, i, ℓ₃, flag, V, store_event_history, event_history, peek)
    if is_on(flag) && (is_boundary_triangle(tri, V) || is_ghost_triangle(V) && !is_boundary_node(tri, new_point)[1])
        # ^ Need to fix the ghost edges if the point is added onto an existing boundary edge. Note that the last 
        #   condition is in case the ghost edges were already correctly added.
        # Also: Why isn't some of this just done using the split_edge! operation?
        u, v = find_edge(tri, V, new_point)
        if is_boundary_edge(tri, u, v) || is_boundary_edge(tri, v, u) # If the edge is not itself a boundary edge, no need to worry.
            if !is_boundary_edge(tri, u, v)
                u, v = v, u
            end
            g = get_adjacent(tri, u, v)
            if !is_true(peek)
                delete_triangle!(tri, v, u, new_point; protect_boundary=true, update_ghost_edges=false)
                delete_triangle!(tri, u, v, g; protect_boundary=true, update_ghost_edges=false)
                add_triangle!(tri, new_point, v, g; update_ghost_edges=false)
                add_triangle!(tri, u, new_point, g; update_ghost_edges=false)
            end
            if is_true(store_event_history)
                trit = triangle_type(tri)
                delete_triangle!(event_history, construct_triangle(trit, u, v, g))
                add_triangle!(event_history, construct_triangle(trit, _new_point, v, g))
                add_triangle!(event_history, construct_triangle(trit, u, _new_point, g))
            end
            if is_constrained(tri) && (contains_boundary_edge(tri, u, v) || contains_boundary_edge(tri, v, u)) # If we don't do this here now, then when we try and do it later we will get a KeyError since we've already modified the boundary edge but we wouldn't have updated the constrained fields
                # We also only do this if contains_boundary_edge since we want to assume that (u, v) does not appear in constrained_edges
                if contains_boundary_edge(tri, u, v)
                    !is_true(peek) && split_boundary_edge!(tri, u, v, new_point)
                    is_true(store_event_history) && split_boundary_edge!(event_history, u, v, _new_point)
                elseif contains_boundary_edge(tri, v, u)
                    !is_true(peek) && split_boundary_edge!(tri, v, u, new_point)
                    is_true(store_event_history) && split_boundary_edge!(event_history, v, u, _new_point)
                end
                E = edge_type(tri)
                # constrained_edges = get_constrained_edges(tri) < -- Don't need this actually, we're just looking at boundary edges here, so no need to consider individually constrained edges
                all_constrained_edges = get_all_constrained_edges(tri)
                if !is_true(peek)
                    delete_edge!(all_constrained_edges, construct_edge(E, u, v))
                    delete_edge!(all_constrained_edges, construct_edge(E, v, u))
                    add_edge!(all_constrained_edges, construct_edge(E, u, new_point))
                    add_edge!(all_constrained_edges, construct_edge(E, new_point, v))
                end
                if is_true(store_event_history)
                    delete_edge!(event_history, construct_edge(E, u, v))
                    add_edge!(event_history, construct_edge(E, u, _new_point))
                    add_edge!(event_history, construct_edge(E, _new_point, v))
                end
            end
        end
    end
    update_representative_point && !is_true(peek) && update_centroid_after_addition!(tri, I(1), q) # How do we efficiently determine which curve to update for a given q when adding into an existing triangulation?
    return nothing
end

@inline function dig_cavity!(tri::Triangulation{P,Ts,I}, r, i, j, ℓ, flag, V, store_event_history=Val(false), event_history=nothing, peek::F=Val(false)) where {P,Ts,I,F}
    if !edge_exists(ℓ)
        # The triangle has already been deleted in this case. 
        return nothing
    end
    _r = is_true(peek) ? num_points(tri) + 1 : r # If we are peeking, then we need to use the number of points in the triangulation as the index for the new point since we don't actually insert the point
    if !contains_constrained_edge(tri, i, j) &&
       is_inside(point_position_relative_to_circumcircle(tri, r, i, j, ℓ)) &&
       !is_boundary_index(ℓ)
        ℓ₁ = get_adjacent(tri, ℓ, i)
        ℓ₂ = get_adjacent(tri, j, ℓ)
        !is_true(peek) && delete_triangle!(tri, j, i, ℓ; protect_boundary=true, update_ghost_edges=false)
        dig_cavity!(tri, r, i, ℓ, ℓ₁, flag, V, store_event_history, event_history, peek)
        dig_cavity!(tri, r, ℓ, j, ℓ₂, flag, V, store_event_history, event_history, peek)
        if is_true(store_event_history)
            trit = triangle_type(tri)
            delete_triangle!(event_history, construct_triangle(trit, j, i, ℓ))
        end
    else
        # If we are here, then this means that we are on an edge of the polygonal cavity. 
        # Note that this also covers the is_boundary_index(ℓ) case.
        if is_on(flag) && contains_constrained_edge(tri, i, j)
            # When we have a constrained edge (i, j) and we add a point r that lies on the edge, 
            # we can run into an issue where we add a degenerate triangle (i, j, r). We need to 
            # check this. There is probably a much smarter way to do this, but this should only 
            # be done very rarely anyway, so I'm not too concerned about the performance hit here. 
            u, v = find_edge(tri, V, r)
            if u == i && v == j
                return nothing
            else
                !is_true(peek) && add_triangle!(tri, r, i, j; update_ghost_edges=false)
                if is_true(store_event_history)
                    trit = triangle_type(tri)
                    add_triangle!(event_history, construct_triangle(trit, _r, i, j))
                end
            end
        else
            !is_true(peek) && add_triangle!(tri, r, i, j; update_ghost_edges=false)
            if is_true(store_event_history)
                trit = triangle_type(tri)
                add_triangle!(event_history, construct_triangle(trit, _r, i, j))
            end
        end
    end
    return nothing
end
