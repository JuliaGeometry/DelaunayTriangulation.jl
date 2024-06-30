"""
    triangulate_curve_bounded(points::P;
    segments=nothing,
    boundary_nodes=nothing,
    IntegerType::Type{I}=Int,
    polygonise_n=4096,
    coarse_n=0,
    check_arguments=true,
    delete_ghosts=false,
    delete_empty_features=true,
    recompute_representative_points=true,
    rng::AbstractRNG=Random.default_rng(),
    insertion_order=nothing, 
    kwargs...) where {P,I} -> Triangulation

Triangulates a curve-bounded domain defined by `(points, segments, boundary_nodes)`. Please
see [`triangulate`](@ref) for a description of the arguments. The only differences are:

- `insertion_order=nothing`: This argument is ignored for curve-bounded domains.
- `polygonise_n=4096`: For generating a high-resolution discretisation of a boundary initially for 
   the construction of a [`PolygonHierarchy`](@ref), many points are needed. This number of points is defined 
   by `polygonise_n`, and must be a power of 2 (otherwise, the next highest power of 2 is used). See [`polygonise`](@ref).
- `coarse_n=0`: This is the number of points to use for initialising a curve-bounded domain via [`coarse_discretisation!`](@ref).
   The default `coarse_n=0` means the discretisation is performed until the maximum variation over any subcurve is 
   less than `π/2`.
- `skip_points`: This is still used, but it is ignored during the enrichment phase (see [`enrich_boundary!`](@ref)).

See also [`BoundaryEnricher`](@ref) and [`enrich_boundary!`](@ref).

!!! note "Refinement"

    To refine the mesh further beyond its initial coarse discretisation, as produced from this function, 
    please see [`refine!`](@ref).
"""
@unstable function triangulate_curve_bounded(points::P;
    segments=nothing,
    boundary_nodes=nothing,
    IntegerType::Type{I}=Int,
    polygonise_n=4096,
    coarse_n=0,
    check_arguments=true,
    delete_ghosts=false,
    delete_empty_features=true,
    recompute_representative_points=true,
    rng::AbstractRNG=Random.default_rng(),
    insertion_order=nothing, # use this so that it gets ignored by the kwargs
    kwargs...) where {P,I}
    enricher = BoundaryEnricher(points, boundary_nodes, segments; IntegerType, n=polygonise_n, coarse_n)
    return _triangulate_curve_bounded(points, enricher;
        IntegerType,
        check_arguments,
        delete_ghosts,
        delete_empty_features,
        recompute_representative_points,
        rng,
        insertion_order,
        kwargs...)
end
@unstable function _triangulate_curve_bounded(points::P, enricher;
    IntegerType::Type{I}=Int,
    check_arguments=true,
    delete_ghosts=false,
    delete_empty_features=true,
    recompute_representative_points=true,
    rng::AbstractRNG=Random.default_rng(),
    insertion_order=nothing, # use this so that it gets ignored by the kwargs
    kwargs...) where {P,I}
    check_arguments && check_args(enricher)
    enrich_boundary!(enricher)
    new_boundary_nodes = get_boundary_nodes(enricher)
    new_segments = get_segments(enricher)
    full_polygon_hierarchy = get_polygon_hierarchy(enricher)
    boundary_curves = get_boundary_curves(enricher)
    tri = triangulate(points;
        IntegerType,
        segments=new_segments,
        boundary_nodes=new_boundary_nodes,
        full_polygon_hierarchy,
        boundary_curves,
        boundary_enricher=enricher,
        check_arguments=false,
        delete_ghosts=false,
        rng,
        kwargs...)
    postprocess_triangulate!(tri; delete_ghosts, delete_empty_features, recompute_representative_points)
    return tri
end

const M_INF = 16 # 2^16
"""
    coarse_discretisation!(points, boundary_nodes, boundary_curve; n=0) 

Constructs an initial coarse discretisation of a curve-bounded domain with bonudary defines by 
`(points, boundary_nodes, boundary_curves)`, where `boundary_nodes` and `boundary_curves` should 
come from [`convert_boundary_curves!`](@ref). The argument `n` is the amount of times to split an edge.
If non-zero, this should be a power of two (otherwise it will be rounded up to the next power of two). If it is 
zero, then the splitting will continue until the maximum total variation over any subcurve is less than π/2.
"""
function coarse_discretisation!(points, boundary_nodes, boundary_curves; n::I=0) where {I}
    !is_curve_bounded(boundary_curves) && return points, boundary_nodes
    if n > 0
        n = max(n, I(4))
        n = !ispow2(n) ? nextpow(2, n) : n
        m = ceil(I, log2(n))
    else
        m = I(M_INF)
    end
    if has_multiple_curves(boundary_nodes)
        _coarse_discretisation_multiple_curves!(points, boundary_nodes, boundary_curves, m)
    elseif has_multiple_sections(boundary_nodes)
        _coarse_discretisation_multiple_sections!(points, boundary_nodes, boundary_curves, m)
    else
        _coarse_discretisation_contiguous!(points, boundary_nodes, boundary_curves, 1, m)
    end
    return points, boundary_nodes
end
function _coarse_discretisation_multiple_curves!(points, boundary_nodes, boundary_curves, m)
    curve_index = 1
    for boundary_curve_index in 1:num_curves(boundary_nodes)
        curve_nodes = get_boundary_nodes(boundary_nodes, boundary_curve_index)
        for section_index in 1:num_sections(curve_nodes)
            section_nodes = get_boundary_nodes(curve_nodes, section_index)
            _coarse_discretisation_contiguous!(points, section_nodes, boundary_curves, curve_index, m)
            curve_index += 1
        end
    end
    return nothing
end
function _coarse_discretisation_multiple_sections!(points, boundary_nodes, boundary_curves, m)
    for section_index in 1:num_sections(boundary_nodes)
        section_nodes = get_boundary_nodes(boundary_nodes, section_index)
        _coarse_discretisation_contiguous!(points, section_nodes, boundary_curves, section_index, m)
    end
    return nothing
end
function _coarse_discretisation_contiguous!(points, boundary_nodes, boundary_curves, curve_index, m)
    is_piecewise_linear(boundary_curves, curve_index) && return nothing
    is_while = m == M_INF
    for _ in 1:m
        nn = num_boundary_edges(boundary_nodes)
        t₁ = 0.0
        u = get_boundary_nodes(boundary_nodes, 1)
        ctr = 1
        max_variation = -Inf
        for _ in 1:nn # not using split_subcurve! here since we have the advantage of knowing that we are splitting consecutive edges
            j = ctr + 1
            v = get_boundary_nodes(boundary_nodes, j)
            q = get_point(points, v)
            t₂ = get_inverse(boundary_curves, curve_index, q)
            if iszero(t₂)
                t₂ = 1.0 # fix for periodic curves 
            end
            t, Δθ = get_equivariation_split(boundary_curves, curve_index, t₁, t₂)
            ct = eval_boundary_curve(boundary_curves, curve_index, t)
            push_point!(points, ct)
            r = num_points(points)
            insert_boundary_node!(boundary_nodes, (boundary_nodes, j), r)
            t₁ = t₂
            u = v
            ctr += 2 # because we inserted a point 
            max_variation = max(Δθ, max_variation)
        end
        is_while && nn ≥ 2 && max_variation < π / 2 && break
    end
    return nothing
end

"""
    enrich_boundary!(enricher::BoundaryEnricher)

Enriches the initial boundary defined inside `enricher`, implementing the algorithm of Gosselin and Ollivier-Gooch (2007).
At the termination of the algorithm, all edges will contain no other points inside their 
diametral circles.
"""
function enrich_boundary!(enricher::BoundaryEnricher)
    queue = get_queue(enricher)
    points = get_points(enricher)
    enqueue_all!(queue, each_point_index(points))
    while !isempty(queue)
        _enrich_boundary_itr!(enricher)
    end
    return enricher
end
function _enrich_boundary_itr!(enricher::BoundaryEnricher)
    queue = get_queue(enricher)
    points = get_points(enricher)
    spatial_tree = get_spatial_tree(enricher)
    v = popfirst!(queue)
    r = get_point(points, v)
    intersections = get_intersections(spatial_tree, v; cache_id=1)
    requeued = false
    for bbox in intersections
        i, j = get_edge(bbox)
        i, j = reorient_edge(enricher, i, j)
        p, q = get_point(points, i, j)
        in_cert = point_position_relative_to_diametral_circle(p, q, r)
        if is_inside(in_cert)
            vis_cert = test_visibility(enricher, i, j, v)
            if is_visible(vis_cert)
                has_precision_issues = split_subcurve!(enricher, i, j)
                has_precision_issues && continue
                !requeued && push!(queue, v)
                requeued = true
            end
        end
    end
    return enricher
end

"""
    split_subcurve!(enricher::BoundaryEnricher, i, j) -> Bool

Splits the curve associated with the edge `(i, j)` into two subcurves by inserting a point `r` between `(i, j)` such that the 
total variation of the subcurve is equal on `(i, r)` and `(r, j)`. The returned value is a `flag` that is `true` 
if there was a precision issue, and `false` otherwise.
"""
function split_subcurve!(enricher::BoundaryEnricher, i, j)
    flag, apex, complex_id, _ = is_small_angle_complex_member(enricher, i, j)
    if !flag
        return _split_subcurve_standard!(enricher, i, j)
    else
        return _split_subcurve_complex!(enricher, apex, complex_id)
    end
end
function _split_subcurve_standard!(enricher::BoundaryEnricher, i, j)
    points = get_points(enricher)
    t, Δθ, ct = compute_split_position(enricher, i, j)
    if isnan(Δθ)
        return true
    end
    push_point!(points, ct)
    r = num_points(points)
    split_edge!(enricher, i, j, r)
    queue = get_queue(enricher)
    push!(queue, r)
    return false
end
function _split_subcurve_complex!(enricher::BoundaryEnricher, apex, complex_id)
    # All members of the complex get split at their intersection with the same circular shell centered at apex 
    complexes = get_small_angle_complexes(enricher, apex)
    complex = complexes[complex_id]
    points = get_points(enricher)
    emin = get_minimum_edge_length(complex, points)
    circle_radius = balanced_power_of_two_ternary_split(emin)
    members = get_members(complex)
    p = get_point(points, apex)
    num_precision_issues = 0
    num_members = length(members)
    for (member_id, member) in enumerate(members)
        split_point = _compute_split_position_complex(enricher, apex, member, circle_radius)
        next_edge = get_next_edge(member)
        q = get_point(points, next_edge)
        if check_split_subsegment_precision(getx(split_point), gety(split_point), p, q)
            num_precision_issues += 1
            continue
        end
        push_point!(points, split_point)
        r = num_points(points)
        i, j = reorient_edge(enricher, apex, next_edge)
        split_edge!(enricher, i, j, r)
        replace_next_edge!(enricher, apex, complex_id, member_id, r)
    end
    return num_precision_issues == num_members
end

"""
    compute_split_position(enricher::BoundaryEnricher, i, j) -> (Float64, Float64, NTuple{2,Float64})

Gets the point to split the edge `(i, j)` at.

# Arguments 
- `enricher::BoundaryEnricher`: The enricher.
- `i`: The first point of the edge.
- `j`: The second point of the edge.

# Outputs
- `t`: The parameter value of the split point.
- `Δθ`: The total variation of the subcurve `(i, t)`. If a split was created due to a small angle, this will be set to zero.
- `ct`: The point to split the edge at.
"""
function compute_split_position(enricher::BoundaryEnricher, i, j)
    num_adjoin, adjoin_vert = has_acute_neighbouring_angles(enricher, i, j)
    if num_adjoin == 0
        t, Δθ, ct = _compute_split_position_standard(enricher, i, j)
    else
        t, Δθ, ct = _compute_split_position_acute(enricher, i, j, num_adjoin, adjoin_vert)
    end
    points = get_points(enricher)
    p, q = get_point(points, i, j)
    if check_split_subsegment_precision(getx(ct), gety(ct), p, q)
        #=
        This happens when the points p and q are so close to each other, that the interval 
        (i-1, i, i+1) found for the binary splitting used in get_closest_point are the exact same, 
        and the splitting point ends up being exactly the same as one of the endpoints. For example,
        take 
        p = (14.000000000006803, 0.0)
        q = (13.999999999998183, 0.0)
        with curve
        curve = CatmullRomSpline([(10.0, -3.0), (20.0, 0.0), (18.0, 0.0), (10.0, 0.0)])
        then 
        t₁ = get_inverse(curve, p) 
        t₂ = get_inverse(curve, q)
        which gives t₁ == t₂.
        =#
        t = NaN
        Δθ = NaN
    end
    return t, Δθ, ct
end
function _compute_split_position_acute(enricher::BoundaryEnricher, i, j, num_adjoin, adjoin_vert)
    points = get_points(enricher)
    p, q = get_point(points, i, j)
    if num_adjoin == 2
        flipped = i ≥ j
        a, b = flipped ? (q, p) : (p, q)
        t = compute_concentric_shell_quarternary_split_position(a, b)
        if flipped
            t = one(t) - t
        end
    else # num_adjoin == 1 
        t = compute_concentric_shell_ternary_split_position(p, q)
        if adjoin_vert ≠ i
            t = one(t) - t
        end
    end
    if abs(t - 1 / 2) < MIDPOINT_TOLERANCE
        t = 1 / 2
    end
    px, py = getxy(p)
    qx, qy = getxy(q)
    parent_curve = get_parent(enricher, i, j)
    if parent_curve == ∅ || is_piecewise_linear(enricher, parent_curve) # == ∅ means the member is associated with an interior segment
        ct = (px + t * (qx - px), py + t * (qy - py))
    else
        ℓ = dist((px, py), (qx, qy))
        circle_radius = t * ℓ
        t₁ = get_inverse(enricher, parent_curve, p)
        t₂ = get_inverse(enricher, parent_curve, q)
        if adjoin_vert ≠ i
            t₁, t₂ = t₂, t₁ # shell needs to be centered appropriately 
        end
        t, ct = get_circle_intersection(enricher, parent_curve, t₁, t₂, circle_radius)
    end
    Δθ = 0.0
    return t, Δθ, ct
end
function _compute_split_position_standard(enricher::BoundaryEnricher, i, j)
    points = get_points(enricher)
    parent_curve = get_parent(enricher, i, j)
    p, q = get_point(points, i, j)
    if parent_curve == ∅ || is_piecewise_linear(enricher, parent_curve) # == ∅ means the member is associated with an interior segment
        px, py = getxy(p)
        qx, qy = getxy(q)
        ct = midpoint((px, py), (qx, qy))
        Δθ = 0.0
        t = NaN # there's no sensible value to apply here. it's not 1/2 since that's for LineSegments, not for PiecewiseLinear
    else
        t₁ = get_inverse(enricher, parent_curve, p)
        t₂ = get_inverse(enricher, parent_curve, q)
        if iszero(t₂) # fix for periodic curves 
            t₂ = one(t₂)
        end
        t, Δθ = get_equivariation_split(enricher, parent_curve, t₁, t₂)
        ct = eval_boundary_curve(enricher, parent_curve, t)
    end
    return t, Δθ, ct
end
function _compute_split_position_complex(enricher::BoundaryEnricher, apex, member, circle_radius)
    points = get_points(enricher)
    p = get_point(points, apex)
    parent_curve = get_parent_curve(member)
    next_edge = get_next_edge(member)
    q = get_point(points, next_edge)
    if parent_curve == ∅ || is_piecewise_linear(enricher, parent_curve) # == ∅ means the member is associated with an interior segment
        px, py = getxy(p)
        qx, qy = getxy(q)
        ℓ = dist((px, py), (qx, qy))
        t = circle_radius / ℓ
        split_point = (px + t * (qx - px), py + t * (qy - py))
    else
        t₁ = get_inverse(enricher, parent_curve, p)
        t₂ = get_inverse(enricher, parent_curve, q)
        t, split_point = get_circle_intersection(enricher, parent_curve, t₁, t₂, circle_radius)
    end
    return split_point
end

"""
    has_acute_neighbouring_angles(enricher::BoundaryEnricher, i, j) -> Int, Vertex

Given a boundary edge `(i, j)`, tests if the neighbouring angles are acute. The first returned value 
is the number of angles adjoining `(i, j)` that are acute (0, 1, or 2). The second returned value is the
vertex that adjoins the edge `(i, j)` that is acute. If there are no such angles, or if there are two, then this 
returned vertex is `$∅`.

(The purpose of this function is similar to [`segment_vertices_adjoin_other_segments_at_acute_angle`](@ref).)
"""
function has_acute_neighbouring_angles(enricher::BoundaryEnricher{P,B,C,I}, i, j) where {P,B,C,I}
    is_segment(enricher, i, j) && return 0, I(∅)
    points = get_points(enricher)
    p, q = get_point(points, i, j)
    pos, idx = get_boundary_edge_map(enricher, i, j) # If idx = 1, then this edge would have already been placed into a small angle complex, so we can assume idx > 1, i.e. indexing at idx - 1 is not a problem.
    boundary_nodes = get_boundary_nodes(enricher)
    section_nodes = get_boundary_nodes(boundary_nodes, pos)
    n = num_boundary_edges(section_nodes)
    (idx == 1 || idx == n) && return 0, I(∅)
    u, v = get_boundary_nodes(section_nodes, idx - 1), get_boundary_nodes(section_nodes, idx + 2) # j is at idx + 1 
    r, s = get_point(points, u), get_point(points, v)
    angle_pqr = opposite_angle(r, q, p)
    angle_qsp = opposite_angle(p, s, q)
    num_adjoin = is_acute(angle_pqr) + is_acute(angle_qsp)
    if num_adjoin == 0 || num_adjoin == 2
        return num_adjoin, I(∅)
    else
        return num_adjoin, is_acute(angle_pqr) ? i : j
    end
end

"""
    test_visibility(enricher::BoundaryEnricher, i, j, k) -> Certificate

Tests if the vertex `k` is visible from the edge `(i, j)`. Returns a [`Certificate`](@ref) which is

- `Invisible`: If `k` is not visible from `(i, j)`.
- `Visible`: If `k` is visible from `(i, j)`.  

For this function, `k` should be inside the diametral circle of `(i, j)`.

We say that `k` is invisibile from `(i, j)` if the edges `(i, k)` or `(j, k)` intersect any other 
boundary edges, or there is a hole between `(i, j)` and `k`.

!!! danger "Definition incompatibility"

    This is not the same definition used in defining constrained Delaunay triangulations, 
    where visibility means visible from ANY point on the edge instead of only from the endpoints.
"""
function test_visibility(enricher::BoundaryEnricher, i, j, k)
    parent_curve_index = get_parent(enricher, i, j)
    points = get_points(enricher)
    boundary_nodes = get_boundary_nodes(enricher)
    boundary_curves = get_boundary_curves(enricher)
    spatial_tree = get_spatial_tree(enricher)
    curve_index_map = get_curve_index_map(enricher)
    polygon_hierarchy = get_polygon_hierarchy(enricher)
    return test_visibility(points, boundary_nodes, boundary_curves, parent_curve_index, spatial_tree, curve_index_map, polygon_hierarchy, i, j, k)
end
function test_visibility(points, boundary_nodes, boundary_curves, parent_curve_index, spatial_tree, curve_index_map, polygon_hierarchy, i, j, k)
    p, q, a = get_point(points, i, j, k)
    if parent_curve_index ≠ ∅ # interior segments don't bound holes, so this first check is useless
        side_e = point_position_relative_to_line(p, q, a)
        is_pl = is_piecewise_linear(boundary_curves, parent_curve_index)
        if !is_pl
            side_c = point_position_relative_to_curve(boundary_curves, parent_curve_index, a)
        else
            # Two possibilities: 
            #   If side_e == Cert.Right, then yes the point is outside of the domain.
            #   If side_e == Cert.Left, then well this still doesn't guarantee that the point is inside the domain. For example,
            #       a take an edge ((0.0, 0.0), (0.25, 0.0))) on a triangle ((0.0,0.0), (1.0,0.0), (1.0,1.0)). The point 
            #       (0.0, 0.05) is outside of the domain, but it is left of (1, 2). 
            if is_right(side_e)
                side_c = side_e
            else
                boundary_curve = curve_index_map[parent_curve_index]
                δ = if has_multiple_curves(boundary_nodes)
                    distance_to_polygon(a, points, get_boundary_nodes(boundary_nodes, boundary_curve))
                else
                    distance_to_polygon(a, points, boundary_nodes)
                end
                is_pos = get_polygon_orientation(polygon_hierarchy, boundary_curve)
                if !is_pos
                    δ *= -one(δ)
                end
                side_c = if δ > zero(δ)
                    Cert.Left
                elseif δ < zero(δ)
                    Cert.Right
                else
                    Cert.On
                end
            end
        end
        if is_collinear(side_e)
            return Cert.Visible
        elseif !is_pl && side_c ≠ side_e
            return Cert.Visible
        elseif is_right(side_c)
            return Cert.Invisible
        end
    end
    int₁ = false
    int₂ = false
    intersections = get_intersections(spatial_tree, i, j, k; cache_id=2)
    for box in intersections
        u, v = get_edge(box)
        !edges_are_disjoint((i, j), (u, v)) && continue
        p′, q′ = get_point(points, u, v)
        if !int₁
            cert = line_segment_intersection_type(p, a, p′, q′)
            int₁ = !has_no_intersections(cert) && !is_touching(cert)
        end
        if !int₂
            cert = line_segment_intersection_type(a, q, p′, q′)
            int₂ = !has_no_intersections(cert) && !is_touching(cert)
        end
        int₁ && int₂ && return Cert.Invisible
    end
    return Cert.Visible
end