struct BoundaryEnricher{P,B,C,I,T,S,E}
    points::P
    boundary_nodes::B
    segments::S
    boundary_curves::C
    polygon_hierarchy::PolygonHierarchy{I}
    parent_map::Dict{E,I}
    curve_index_map::Dict{I,I}
    boundary_edge_map::BoundaryEdgeMap{E,T}
    spatial_tree::BoundaryRTree{P}
    queue::Queue{I}
    small_angle_complexes::Dict{I,Vector{SmallAngleComplex{I}}}
end

@inline Base.copy(enricher::BoundaryEnricher) = enrcopy(enricher)
@inline enrcopy(::Nothing; kwargs...) = nothing

function enrcopy(enricher::BoundaryEnricher; 
        points = copy(get_points(enricher)), 
        boundary_nodes  = copy(get_boundary_nodes(enricher)),
        segments = copy(get_segments(enricher)),
        boundary_curves = _plcopy.(get_boundary_curves(enricher); points),
        polygon_hierarchy = copy(get_polygon_hierarchy(enricher)),
        boundary_edge_map = copy(get_boundary_edge_map(enricher))) 
    parent_map = copy(get_parent_map(enricher))
    curve_index_map = copy(get_curve_index_map(enricher))
    spatial_tree = copy(get_spatial_tree(enricher))
    queue = copy(get_queue(enricher))
    small_angle_complexes = copy(get_small_angle_complexes(enricher))
    return BoundaryEnricher(points, boundary_nodes, segments, boundary_curves,
        polygon_hierarchy, parent_map, curve_index_map,
        boundary_edge_map, spatial_tree, queue, small_angle_complexes)
end

function BoundaryEnricher(
    points::P, boundary_nodes::B, segments=nothing; n=4096, coarse_n=0,
    IntegerType::Type{I}=Int,
    EdgeType::Type{E}=isnothing(segments) ? NTuple{2,IntegerType} : (edge_type ∘ typeof)(segments),
    EdgesType::Type{Es}=isnothing(segments) ? Set{EdgeType} : typeof(segments),
) where {P,B,I,E,Es}
    boundary_curves, new_boundary_nodes = convert_boundary_curves!(points, boundary_nodes, I)
    polygon_hierarchy = construct_polygon_hierarchy(points, new_boundary_nodes, boundary_curves; IntegerType, n)
    return _construct_boundary_enricher(points, new_boundary_nodes, boundary_curves, polygon_hierarchy, segments, n, coarse_n, I, E, Es)
end
function _construct_boundary_enricher(points, boundary_nodes, boundary_curves, polygon_hierarchy, segments, n, coarse_n, ::Type{I}, ::Type{E}, ::Type{Es}) where {I,E,Es}
    expand_bounds!(polygon_hierarchy, ε(number_type(points)))
    coarse_discretisation!(points, boundary_nodes, boundary_curves; n=coarse_n)
    boundary_edge_map = construct_boundary_edge_map(boundary_nodes, I)
    parent_map = Dict{E,I}()
    curve_index_map = Dict{I,I}()
    spatial_tree = BoundaryRTree(points)
    queue = Queue{I}()
    small_angle_complexes = get_small_angle_complexes(points, boundary_nodes, boundary_curves, segments; IntegerType=I)
    _segments = isnothing(segments) ? Es() : segments
    enricher = BoundaryEnricher(points, boundary_nodes, _segments, boundary_curves, polygon_hierarchy, parent_map, curve_index_map, boundary_edge_map, spatial_tree, queue, small_angle_complexes)
    construct_parent_map!(enricher)
    construct_curve_index_map!(enricher)
    construct_tree!(enricher)
    return enricher
end

function Base.:(==)(enricher1::BoundaryEnricher, enricher2::BoundaryEnricher)
    get_points(enricher1) ≠ get_points(enricher2) && return false
    get_boundary_nodes(enricher1) ≠ get_boundary_nodes(enricher2) && return false
    get_boundary_curves(enricher1) ≠ get_boundary_curves(enricher2) && return false
    get_polygon_hierarchy(enricher1) ≠ get_polygon_hierarchy(enricher2) && return false
    get_parent_map(enricher1) ≠ get_parent_map(enricher2) && return false
    get_curve_index_map(enricher1) ≠ get_curve_index_map(enricher2) && return false
    get_boundary_edge_map(enricher1) ≠ get_boundary_edge_map(enricher2) && return false
    get_spatial_tree(enricher1) ≠ get_spatial_tree(enricher2) && return false
    get_queue(enricher1) ≠ get_queue(enricher2) && return false
    get_small_angle_complexes(enricher1) ≠ get_small_angle_complexes(enricher2) && return false
    return true
end

function check_args(enricher::BoundaryEnricher)
    points = get_points(enricher)
    boundary_nodes = get_boundary_nodes(enricher)
    hierarchy = get_polygon_hierarchy(enricher)
    boundary_curves = get_boundary_curves(enricher)
    return check_args(points, boundary_nodes, hierarchy, boundary_curves)
end

@inline get_points(boundary_enricher::BoundaryEnricher) = boundary_enricher.points
@inline get_boundary_nodes(boundary_enricher::BoundaryEnricher) = boundary_enricher.boundary_nodes
@inline get_boundary_curves(boundary_enricher::BoundaryEnricher) = boundary_enricher.boundary_curves
@inline get_segments(boundary_enricher::BoundaryEnricher) = boundary_enricher.segments
@inline function has_segments(boundary_enricher::BoundaryEnricher)
    segments = get_segments(boundary_enricher)
    isnothing(segments) && return false
    return !isempty(segments)
end
@inline function is_segment(enricher::BoundaryEnricher, i, j)
    !has_segments(enricher) && return false
    segments = get_segments(enricher)
    E = edge_type(segments)
    e = construct_edge(E, i, j)
    return contains_unoriented_edge(e, segments)
end

@inline get_boundary_curve(boundary_enricher::BoundaryEnricher, curve_index) = get_boundary_curves(boundary_enricher)[curve_index]
@inline get_polygon_hierarchy(boundary_enricher::BoundaryEnricher) = boundary_enricher.polygon_hierarchy
@inline get_parent_map(boundary_enricher::BoundaryEnricher) = boundary_enricher.parent_map
@inline get_curve_index_map(boundary_enricher::BoundaryEnricher) = boundary_enricher.curve_index_map
@inline map_curve_index(boundary_enricher::BoundaryEnricher, curve_index) = get_curve_index_map(boundary_enricher)[curve_index]
@inline get_boundary_edge_map(boundary_enricher::BoundaryEnricher) = boundary_enricher.boundary_edge_map
@inline get_boundary_edge_map(boundary_enricher::BoundaryEnricher, i, j) = get_boundary_edge_map(boundary_enricher)[(i, j)]
@inline get_spatial_tree(boundary_enricher::BoundaryEnricher) = boundary_enricher.spatial_tree
@inline get_queue(boundary_enricher::BoundaryEnricher) = boundary_enricher.queue
@inline get_small_angle_complexes(boundary_enricher::BoundaryEnricher) = boundary_enricher.small_angle_complexes
@inline is_small_angle_complex_apex(boundary_enricher::BoundaryEnricher, apex) = haskey(get_small_angle_complexes(boundary_enricher), apex)
@inline function get_small_angle_complexes(boundary_enricher::BoundaryEnricher, apex)
    complexes = get_small_angle_complexes(boundary_enricher)
    return complexes[apex]
end
@inline get_parent(boundary_enricher::BoundaryEnricher{P,B,C,I}, i, j) where {P,B,C,I} = Base.get(get_parent_map(boundary_enricher), (i, j), I(∅))

@inline function set_parent!(boundary_enricher::BoundaryEnricher, i, j, k)
    get_parent_map(boundary_enricher)[(i, j)] = k
    return nothing
end

@inline function delete_edge!(boundary_enricher::BoundaryEnricher, i, j)
    delete!(get_parent_map(boundary_enricher), (i, j))
    return nothing
end

@inline function update_parent_map!(boundary_enricher::BoundaryEnricher, i, j, k)
    parent = get_parent(boundary_enricher, i, j)
    delete_edge!(boundary_enricher, i, j)
    set_parent!(boundary_enricher, i, k, parent)
    set_parent!(boundary_enricher, k, j, parent)
    return nothing
end

@inline function each_boundary_edge(enricher::BoundaryEnricher)
    return keys(get_parent_map(enricher))
end

function construct_parent_map!(enricher::BoundaryEnricher)
    parent_map = get_parent_map(enricher)
    empty!(parent_map)
    boundary_nodes = get_boundary_nodes(enricher)
    if has_multiple_curves(boundary_nodes)
        _construct_parent_map_multiple_curves!(enricher)
    elseif has_multiple_sections(boundary_nodes)
        _construct_parent_map_multiple_sections!(enricher)
    else
        _construct_parent_map_contiguous!(enricher)
    end
    return enricher
end
function _construct_parent_map_multiple_curves!(enricher::BoundaryEnricher)
    ctr = 1
    boundary_nodes = get_boundary_nodes(enricher)
    for curve_index in 1:num_curves(boundary_nodes)
        curve_nodes = get_boundary_nodes(boundary_nodes, curve_index)
        for section_index in 1:num_sections(curve_nodes)
            section_nodes = get_boundary_nodes(curve_nodes, section_index)
            _construct_parent_map_contiguous!(enricher, section_nodes, ctr)
            ctr += 1
        end
    end
    return enricher
end
function _construct_parent_map_multiple_sections!(enricher::BoundaryEnricher, boundary_nodes=get_boundary_nodes(enricher), ctr=1)
    ns = num_sections(boundary_nodes)
    for i in 1:ns
        section_nodes = get_boundary_nodes(boundary_nodes, i)
        _construct_parent_map_contiguous!(enricher, section_nodes, ctr)
        ctr += 1
    end
    return enricher
end
function _construct_parent_map_contiguous!(enricher::BoundaryEnricher, boundary_nodes=get_boundary_nodes(enricher), ctr=1)
    n = num_boundary_edges(boundary_nodes)
    for i in 1:n
        u = get_boundary_nodes(boundary_nodes, i)
        v = get_boundary_nodes(boundary_nodes, i + 1)
        set_parent!(enricher, u, v, ctr)
    end
    return enricher
end

function construct_curve_index_map!(enricher::BoundaryEnricher)
    boundary_nodes = get_boundary_nodes(enricher)
    if has_multiple_curves(boundary_nodes)
        _construct_curve_index_map_multiple_curves!(enricher)
    elseif has_multiple_sections(boundary_nodes)
        _construct_curve_index_map_multiple_sections!(enricher)
    else
        _construct_curve_index_map_contiguous!(enricher)
    end
    return enricher
end
function _construct_curve_index_map_multiple_curves!(enricher::BoundaryEnricher)
    boundary_nodes = get_boundary_nodes(enricher)
    curve_index_map = get_curve_index_map(enricher)
    ctr = 1
    for j in 1:num_curves(boundary_nodes)
        curve_nodes = get_boundary_nodes(boundary_nodes, j)
        for _ in 1:num_sections(curve_nodes)
            curve_index_map[ctr] = j
            ctr += 1
        end
    end
    return enricher
end
function _construct_curve_index_map_multiple_sections!(enricher::BoundaryEnricher)
    boundary_nodes = get_boundary_nodes(enricher)
    curve_index_map = get_curve_index_map(enricher)
    for i in 1:num_sections(boundary_nodes)
        curve_index_map[i] = 1
    end
    return enricher
end
function _construct_curve_index_map_contiguous!(enricher::BoundaryEnricher)
    curve_index_map = get_curve_index_map(enricher)
    curve_index_map[1] = 1
    return enricher
end

@inline function is_piecewise_linear(enricher::BoundaryEnricher, curve_index)
    boundary_curves = get_boundary_curves(enricher)
    return is_piecewise_linear(boundary_curves, curve_index)
end
@inline function is_piecewise_linear(boundary_curves::C, curve_index) where {C<:Tuple}
    isempty(boundary_curves) && return true
    return eval_fnc_at_het_tuple_element(is_piecewise_linear, boundary_curves, curve_index)
end

@inline function is_linear(enricher::BoundaryEnricher, curve_index)
    boundary_curves = get_boundary_curves(enricher)
    return is_linear(boundary_curves, curve_index)
end
@inline function is_linear(boundary_curves::C, curve_index) where {C<:Tuple}
    isempty(boundary_curves) && return true
    return eval_fnc_at_het_tuple_element(is_linear, boundary_curves, curve_index)
end

@inline function get_inverse(enricher::BoundaryEnricher, curve_index, q)
    boundary_curves = get_boundary_curves(enricher)
    return get_inverse(boundary_curves, curve_index, q)
end
@inline function get_inverse(boundary_curves::C, curve_index, q) where {C<:Tuple}
    return eval_fnc_at_het_tuple_element_with_arg(get_inverse, boundary_curves, (q,), curve_index)
end

@inline function get_equivariation_split(enricher::BoundaryEnricher, curve_index, t₁, t₂)
    boundary_curves = get_boundary_curves(enricher)
    return get_equivariation_split(boundary_curves, curve_index, t₁, t₂)
end
@inline function get_equivariation_split(boundary_curves::C, curve_index, t₁, t₂) where {C<:Tuple}
    return eval_fnc_at_het_tuple_element_with_arg(get_equivariation_split, boundary_curves, (t₁, t₂), curve_index)
end

@inline function get_equidistant_split(enricher::BoundaryEnricher, curve_index, t₁, t₂)
    boundary_curves = get_boundary_curves(enricher)
    return get_equidistant_split(boundary_curves, curve_index, t₁, t₂)
end
@inline function get_equidistant_split(boundary_curves::C, curve_index, t₁, t₂) where {C<:Tuple}
    return eval_fnc_at_het_tuple_element_with_arg(get_equidistant_split, boundary_curves, (t₁, t₂), curve_index)
end

@inline function eval_boundary_curve(enricher::BoundaryEnricher, curve_index, t)
    boundary_curves = get_boundary_curves(enricher)
    return eval_boundary_curve(boundary_curves, curve_index, t)
end
@inline function eval_boundary_curve(boundary_curves::C, curve_index, t) where {C<:Tuple}
    return eval_fnc_in_het_tuple(boundary_curves, t, curve_index)
end

@inline @optarg1 DEFAULT_KERNEL function point_position_relative_to_curve(kernel::AbstractPredicateKernel, enricher::BoundaryEnricher, curve_index, p)
    boundary_curves = get_boundary_curves(enricher)
    return point_position_relative_to_curve(kernel, boundary_curves, curve_index, p)
end
@inline @optarg1 DEFAULT_KERNEL function point_position_relative_to_curve(kernel::AbstractPredicateKernel, boundary_curves::C, curve_index, p) where {C<:Tuple}
    return eval_fnc_at_het_tuple_element_with_arg_and_prearg(point_position_relative_to_curve, boundary_curves, kernel, (p,), curve_index)
end

function angle_between(enricher::BoundaryEnricher, curve_index1, curve_index2)
    boundary_curves = get_boundary_curves(enricher)
    return angle_between(boundary_curves, curve_index1, curve_index2)
end
function angle_between(boundary_curves::Tuple, curve_index1, curve_index2)
    return eval_fnc_at_het_tuple_two_elements(angle_between, boundary_curves, curve_index1, curve_index2)
end

function get_circle_intersection(enricher::BoundaryEnricher, curve_index, t₁, t₂, r)
    boundary_curves = get_boundary_curves(enricher)
    return get_circle_intersection(boundary_curves, curve_index, t₁, t₂, r)
end
function get_circle_intersection(boundary_curves::Tuple, curve_index, t₁, t₂, r)
    return eval_fnc_at_het_tuple_element_with_arg(get_circle_intersection, boundary_curves, (t₁, t₂, r), curve_index)
end

function polygonise(points, boundary_nodes, boundary_curves; n=4096)
    new_points = deepcopy(points)
    new_boundary_nodes = deepcopy(boundary_nodes)
    coarse_discretisation!(new_points, new_boundary_nodes, boundary_curves; n)
    return new_points, new_boundary_nodes
end

function construct_tree!(enricher::BoundaryEnricher)
    tree = get_spatial_tree(enricher)
    parent_map = get_parent_map(enricher)
    for (i, j) in keys(parent_map)
        insert!(tree, i, j)
    end
    segments = get_segments(enricher)
    if !isnothing(segments)
        for e in segments
            i, j = edge_vertices(e)
            insert!(tree, i, j)
        end
    end
    return enricher
end

function reorient_edge(enricher::BoundaryEnricher, i, j)
    boundary_edge_map = get_boundary_edge_map(enricher)
    if haskey(boundary_edge_map, (i, j))
        return (i, j)
    else
        return (j, i)
    end
end

function split_boundary_edge!(enricher::BoundaryEnricher, i, j, r, update_boundary_nodes=Val(true))
    boundary_nodes = get_boundary_nodes(enricher)
    boundary_edge_map = get_boundary_edge_map(enricher)
    spatial_tree = get_spatial_tree(enricher)
    if is_true(update_boundary_nodes)
        pos = get_boundary_edge_map(enricher, i, j)
        new_pos = (pos[1], pos[2] + 1)
        insert_boundary_node!(boundary_nodes, new_pos, r)
        split_boundary_edge_map!(boundary_edge_map, boundary_nodes, pos, i, j)
    end
    split_edge!(spatial_tree, i, j, r)
    update_parent_map!(enricher, i, j, r)
    return enricher
end

function split_interior_segment!(enricher::BoundaryEnricher, i, j, r, update_segments=Val(true))
    segments = get_segments(enricher)
    spatial_tree = get_spatial_tree(enricher)
    E = edge_type(segments)
    e = construct_edge(E, i, j)
    e1 = construct_edge(E, i, r)
    e2 = construct_edge(E, r, j)
    if is_true(update_segments)
        delete_unoriented_edge!(segments, e)
        add_edge!(segments, e1, e2)
    end
    split_edge!(spatial_tree, i, j, r)
    return enricher
end

function split_edge!(enricher::BoundaryEnricher, i, j, r, update_boundary_nodes=Val(true), update_segments=Val(true), is_interior=is_segment(enricher, i, j))
    if is_interior
        split_interior_segment!(enricher, i, j, r, update_segments)
    else
        split_boundary_edge!(enricher, i, j, r, update_boundary_nodes)
    end
    return enricher
end

function is_small_angle_complex_member(enricher::BoundaryEnricher, i, j)
    if !is_small_angle_complex_apex(enricher, i) && !is_small_angle_complex_apex(enricher, j)
        return false, oftype(i, 0), 0, 0
    end
    for r in (i, j)
        r′ = r == i ? j : i
        !is_small_angle_complex_apex(enricher, r) && continue
        r_complexes = get_small_angle_complexes(enricher, r)
        for (complex_id, complex) in enumerate(r_complexes)
            for (member_id, member) in enumerate(get_members(complex))
                if get_next_edge(member) == r′
                    return true, r, complex_id, member_id
                end
            end
        end
    end
    return false, oftype(i, 0), 0, 0
end

function replace_next_edge!(enricher::BoundaryEnricher, apex, complex_id, member_id, next_edge)
    complexes = get_small_angle_complexes(enricher, apex)
    complex = complexes[complex_id]
    replace_next_edge!(complex, member_id, next_edge)
    return enricher
end
