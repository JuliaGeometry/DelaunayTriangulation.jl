# Comparisons that do not depend on the order of hashed collections, e.g. polygons that may start at a different vertex.

export compare_trees, sort_edge_vector, compare_edge_vectors, _compare_pairs, get_child_from_tree, traverse_tree
export is_rotation, same_elements, same_points, same_polygon, same_cycle, polygon_coordinates, polygon_triangles
export why_not_equal

# Many collections are hashed, so the order in which things are found (and hence e.g. the numbering
# of Voronoi circumcenters, or where a polygon's vertex list starts) depends on the Julia version. These helpers compare
# results in a way that does not depend on that order.

"Returns `true` if `a` is a rotation of `b`."
is_rotation(a, b) = length(a) == length(b) && (isempty(a) || any(k -> circshift(b, k) == a, eachindex(b)))

"Returns `true` if `a` and `b` contain the same elements with the same multiplicities, in any order."
same_elements(a, b) = length(a) == length(b) && all(x -> count(==(x), a) == count(==(x), b), a)

"Returns `true` if the points in `A` and `B` are approximately equal, in any order."
same_points(A, B) = length(A) == length(B) && sort(collect.(A)) ≈ sort(collect.(B))

"Returns `true` if the closed polygons `A` and `B`, given as coordinates, are approximately equal up to where they start."
same_polygon(A, B) = DT.circular_equality(collect.(A), collect.(B), ≈)

"Returns `true` if the points in `A` are approximately a rotation of the points in `B` (neither list repeating its first point at the end)."
same_cycle(A, B) = same_polygon([A; A[1:1]], [B; B[1:1]])

"Returns the coordinates of the vertices of the `i`th polygon of `vorn`."
polygon_coordinates(vorn, i) = [get_polygon_point(vorn, j) for j in get_polygon(vorn, i)]

"""
    polygon_triangles(vorn, i)

Returns, for each vertex of the `i`th polygon of `vorn`, the (sorted) triangle whose circumcenter it is. Unlike the
circumcenter indices themselves, these do not depend on hash order, and they also identify the ghost vertices of
unbounded polygons.
"""
polygon_triangles(vorn, i) = [DT.sort_triangle(DT.get_circumcenter_to_triangle(vorn, c)) for c in get_polygon(vorn, i)]

function sort_edge_vector(E)
    sorted_E = similar(E)
    for i in eachindex(E)
        u, v = E[i]
        e = (min(u, v), max(u, v))
        sorted_E[i] = e
    end
    return sort(sorted_E)
end

function compare_edge_vectors(E1, E2)
    E1s = sort_edge_vector(collect(E1))
    E2s = sort_edge_vector(collect(E2))
    return E1s == E2s
end

function Base.:(==)(stats1::DT.TriangulationStatistics, stats2::DT.TriangulationStatistics)
    for f in fieldnames(DT.TriangulationStatistics)
        if f ≠ :individual_statistics
            if getfield(stats1, f) ≠ getfield(stats2, f)
                return false
            end
        else
            indiv_dict1 = stats1.individual_statistics
            indiv_dict2 = stats2.individual_statistics
            if length(indiv_dict1) ≠ length(indiv_dict2)
                return false
            end
            for (T, v) in indiv_dict1
                V, flag = DT.contains_triangle(T, keys(indiv_dict2))
                if !flag
                    return false
                end
                if v ≠ indiv_dict2[V]
                    return false
                end
            end
        end
    end
    return true
end

function _lexicographically_sort_pair_vector(pairs)
    lt = (x, y) -> let (key_x, value_x) = x, (key_y, value_y) = y
        if value_x < value_y
            return true
        elseif value_x > value_y
            return false
        else
            return key_x < key_y
        end
    end
    return sort(pairs, lt=lt)
end

function _compare_pairs(pairs, dt_pairs)
    @test length(pairs) == length(dt_pairs)
    @test _lexicographically_sort_pair_vector(pairs) == _lexicographically_sort_pair_vector(dt_pairs)
end

function why_not_equal(tri1, tri2)
    !DT.has_ghost_triangles(tri1) && DT.has_ghost_triangles(tri2) && println("!has_ghost_triangles(tri1) && has_ghost_triangles(tri2)")
    !DT.has_ghost_triangles(tri2) && DT.has_ghost_triangles(tri1) && println("!has_ghost_triangles(tri2) && has_ghost_triangles(tri1)")
    DT.get_points(tri1) ≠ DT.get_points(tri2) && println("get_points(tri1) ≠ get_points(tri2)")
    !DT.compare_triangle_collections(DT.get_triangles(tri1), DT.get_triangles(tri2)) && println("!compare_triangle_collections(get_triangles(tri1), get_triangles(tri2))")
    !DT.compare_unoriented_edge_collections(DT.get_interior_segments(tri1), DT.get_interior_segments(tri2)) && println("!compare_unoriented_edge_collections(get_interior_segments(tri1), get_interior_segments(tri2))")
    !DT.compare_unoriented_edge_collections(DT.get_all_segments(tri1), DT.get_all_segments(tri2)) && println("!compare_unoriented_edge_collections(get_all_segments(tri1), get_all_segments(tri2))")
    DT.get_adjacent(tri1) ≠ get_adjacent(tri2) && println("get_adjacent(tri1) ≠ get_adjacent(tri2)")
    DT.get_adjacent2vertex(tri1) ≠ DT.get_adjacent2vertex(tri2) && println("get_adjacent2vertex(tri1) ≠ get_adjacent2vertex(tri2)")
    DT.get_graph(tri1) ≠ DT.get_graph(tri2) && println("get_graph(tri1) ≠ get_graph(tri2)")
    DT.get_boundary_edge_map(tri1) ≠ DT.get_boundary_edge_map(tri2) && println("get_boundary_edge_map(tri1) ≠ get_boundary_edge_map(tri2)")
    DT.get_ghost_vertex_map(tri1) ≠ DT.get_ghost_vertex_map(tri2) && println("get_ghost_vertex_map(tri1) ≠ get_ghost_vertex_map(tri2)")
    DT.get_ghost_vertex_ranges(tri1) ≠ DT.get_ghost_vertex_ranges(tri2) && println("get_ghost_vertex_ranges(tri1) ≠ get_ghost_vertex_ranges(tri2)")
    DT.get_convex_hull(tri1) ≠ DT.get_convex_hull(tri2) && println("get_convex_hull(tri1) ≠ get_convex_hull(tri2)")
    rep1 = DT.get_representative_point_list(tri1)
    rep2 = DT.get_representative_point_list(tri2)
    length(rep1) ≠ length(rep2) && println("length(rep1) ≠ length(rep2)")
    for i in 1:DT.num_curves(tri1)
        p1 = DT.get_representative_point_coordinates(tri1, i)
        p2 = DT.get_representative_point_coordinates(tri2, i)
        !([getx(p1), gety(p1)] ≈ [getx(p2), gety(p2)]) && println("!([getx(p1), gety(p1)] ≈ [getx(p2), gety(p2)]) for curve $i")
    end
    DT.get_polygon_hierarchy(tri1) ≠ DT.get_polygon_hierarchy(tri2) && println("get_polygon_hierarchy(tri1) ≠ get_polygon_hierarchy(tri2)")
    DT.get_boundary_nodes(tri1) ≠ DT.get_boundary_nodes(tri2) && println("get_boundary_nodes(tri1) ≠ get_boundary_nodes(tri2)")
    DT.get_weights(tri1) ≠ DT.get_weights(tri2) && println("get_weights(tri1) ≠ get_weights(tri2)")
    DT.has_ghosts(tri1) ≠ DT.has_ghosts(tri2) && println("has_ghosts(tri1) ≠ has_ghosts(tri2)")
    DT.get_boundary_vertex_to_ghost(tri1) ≠ DT.get_boundary_vertex_to_ghost(tri2) && println("get_boundary_vertex_to_ghost(tri1) ≠ get_boundary_vertex_to_ghost(tri2)")
end

function get_child_from_tree(tree::DT.PolygonTree, i)
    return first(filter(subtree -> DT.get_index(subtree) == i, DT.get_children(tree)))
end

function traverse_tree(tree::DT.PolygonTree)
    height = DT.get_height(tree)
    index = DT.get_index(tree)
    parent_index = DT.has_parent(tree) ? DT.get_index(DT.get_parent(tree)) : 0
    tup = (height, index, parent_index)
    children = DT.get_children(tree)
    schildren = sort(OrderedSet(children), by=x -> DT.get_index(x))
    for child in schildren
        tup = (tup..., traverse_tree(child))
    end
    return tup
end

function traverse_tree(hierarchy::DT.PolygonHierarchy)
    dict = Dict()
    trees = DT.get_trees(hierarchy)
    strees = sort(OrderedDict(trees), by=x -> x[1])
    for (index, tree) in strees
        dict[index] = traverse_tree(tree)
    end
    return dict
end

function compare_trees(hierarchy1::DT.PolygonHierarchy, hierarchy2::DT.PolygonHierarchy)
    DT.get_bounding_boxes(hierarchy1) ≠ DT.get_bounding_boxes(hierarchy2) && return false
    DT.get_polygon_orientations(hierarchy1) ≠ DT.get_polygon_orientations(hierarchy2) && return false
    trav1 = traverse_tree(hierarchy1)
    trav2 = traverse_tree(hierarchy2)
    trav1 ≠ trav2 && return false
    return true
end

function Base.:(==)(cache1::DT.TriangulationCache, cache2::DT.TriangulationCache)
    DT.get_triangulation(cache1) ≠ DT.get_triangulation(cache2) && return false
    DT.get_triangulation_2(cache1) ≠ DT.get_triangulation_2(cache2) && return false
    DT.get_marked_vertices(cache1) ≠ DT.get_marked_vertices(cache2) && return false
    DT.get_interior_segments_on_hull(cache1) ≠ DT.get_interior_segments_on_hull(cache2) && return false
    DT.get_surrounding_polygon(cache1) ≠ DT.get_surrounding_polygon(cache2) && return false
    DT.get_fan_triangles(cache1) ≠ DT.get_fan_triangles(cache2) && return false
    length(DT.get_incircle_cache(cache1)) ≠ length(DT.get_incircle_cache(cache2)) && return false # no point comparing elements 
    typeof(DT.get_incircle_cache(cache1)) ≠ typeof(DT.get_incircle_cache(cache2)) && return false
    length(DT.get_orient3_cache(cache1)) ≠ length(DT.get_orient3_cache(cache2)) && return false # no point comparing elements
    typeof(DT.get_orient3_cache(cache1)) ≠ typeof(DT.get_orient3_cache(cache2)) && return false
    length(DT.get_insphere_cache(cache1)) ≠ length(DT.get_insphere_cache(cache2)) && return false # no point comparing elements
    typeof(DT.get_insphere_cache(cache1)) ≠ typeof(DT.get_insphere_cache(cache2)) && return false
    return true
end
