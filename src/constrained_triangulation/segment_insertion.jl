"""
    vertex_is_closer_than_neighbours(tri::Triangulation, u, v, jᵢ, jᵢ₋₁, jᵢ₊₁)
    vertex_is_closer_than_neighbours(tri::Triangulation, V, u, v, j, shuffled_indices, prev, next)

Given a triangulation `tri` and a line through points with indices `u` and `v`,
tests if the point with index `jᵢ` is closer to the line than those with index 
`jᵢ₋₁` and `jᵢ₊₁`, assuming all these points are to the left of the line. The second 
method extracts these latter two indices using the linked list `(prev, next)` of vertices 
and a shuffled set of indices `shuffled_indices` together with the original vertex list `V`.

This function is useful for constrained triangulations since the algorithm 
used will not work if a point being inserted on the cavity has interior angle 
of 360° or greater. This is possible only if a vertex is closer to the line than 
its neighbours on the polygon.
"""
function vertex_is_closer_than_neighbours(tri::Triangulation, u, v, jᵢ, jᵢ₋₁, jᵢ₊₁)
    prev_comp = point_closest_to_line(tri, u, v, jᵢ, jᵢ₋₁)
    is_further(prev_comp) && return false
    next_comp = point_closest_to_line(tri, u, v, jᵢ, jᵢ₊₁)
    is_further(next_comp) && return false
    return true
end
function vertex_is_closer_than_neighbours(tri::Triangulation, V, u, v, j, shuffled_indices, prev, next)
    jᵢ, jᵢ₋₁, jᵢ₊₁ = index_shuffled_linked_list(V, next, prev, shuffled_indices, j)
    return vertex_is_closer_than_neighbours(tri, u, v, jᵢ, jᵢ₋₁, jᵢ₊₁)
end

"""
    select_random_vertex(tri::Triangulation, V, shuffled_indices, prev, next, range, u, v, rng::AbstractRNG=Random.default_rng())

Given a triangulation `tri`, a line through points with indices `u` and `v`, a shuffled set of 
indices `shuffled_indices`, and a doubly-linked list `(prev, next)` of vertex indices, 
selects a random vertex `j ∈ range` that is not closer to the line than both of its 
neighbours. `V` is the original vertex list.
"""
function select_random_vertex(tri::Triangulation, V, shuffled_indices, prev, next, range, u, v, rng::AbstractRNG=Random.default_rng())
    j = rand(rng, range)
    while vertex_is_closer_than_neighbours(tri, V, u, v, j, shuffled_indices, prev, next)
        j = rand(rng, range)
    end
    return j
end

"""
    update_vertex_linked_list!(shuffled_indices, prev, next, i, j)

Let `π = shuffled_indices`. This function replaces `next[prev[π[j]]]` with `next[π[j]]`,
`prev[next[π[j]]]` with `prev[π[j]]`, and interchanges `π[i]` and `π[j]`. This has the act 
of deleting `V[π[j]]` from the polygon, where `V` is the list of polygon vertices of the 
polygon being evacuated during segment insertion for a constrained triangulation.
"""
function update_vertex_linked_list!(shuffled_indices, prev, next, i, j)
    next[prev[shuffled_indices[j]]] = next[shuffled_indices[j]]
    prev[next[shuffled_indices[j]]] = prev[shuffled_indices[j]]
    shuffled_indices[i], shuffled_indices[j] = shuffled_indices[j], shuffled_indices[i] # move π[j] to follow the live vertices
end

"""
    prepare_vertex_linked_list(V::AbstractArray{I}) where {I}

Given a list of polygon vertices `V`, defines a linked list `(prev, next)`
of polygon vertices so that `(prev[i], i, next[i])` define a trio of polygon vertices 
in counter-clockwise order, and defines, and returns `shuffled_indices` which is currently 
unshuffled.

The first and last entries of the returned values `(prev, next, shuffled_indices)` will 
not be populated, instead being `0`.
"""
function prepare_vertex_linked_list(V::AbstractArray{I}) where {I}
    Base.require_one_based_indexing(V)
    m = length(V)
    next = zeros(I, m)
    prev = zeros(I, m)
    shuffled_indices = zeros(I, m)
    for i = 2:(m-1)
        next[i] = i + 1
        prev[i] = i - 1
        shuffled_indices[i] = i
    end
    return prev, next, shuffled_indices
end

"""
    delete_polygon_vertices_in_random_order!(tri::Triangulation, V, shuffled_indices, prev, next, u, v, rng::AbstractRNG=Random.default_rng())

Given a triangulation `tri`, a vertex list `V`, a set of `shuffled_indices`, a linked list `(prev, next)` for the 
poylgon vertices, and a segment `(u, v)` that was inserted in order to define the polygon `V`, deletes vertices of `V`,
via their representation in `(prev, next, shuffled_indices)`, in a random order.
"""
function delete_polygon_vertices_in_random_order!(tri::Triangulation, V, shuffled_indices, prev, next, u, v, rng::AbstractRNG=Random.default_rng())
    m = length(V)
    for i in (m-1):-1:3
        j = select_random_vertex(tri, V, shuffled_indices, prev, next, 2:i, u, v, rng)
        update_vertex_linked_list!(shuffled_indices, prev, next, i, j)
    end
    return nothing
end

"""
    triangulate_cavity_cdt(points, V;
        IntegerType::Type{I}=Int64,
        EdgeType::Type{E}=NTuple{2,IntegerType},
        TriangleType::Type{Vs}=NTuple{3,IntegerType},
        EdgesType::Type{Es}=Set{EdgeType},
        TrianglesType::Type{Ts}=Set{TriangleType},
        rng::AbstractRNG=Random.default_rng()) where {I,E,Vs,Es,Ts}
    triangulate_cavity_cdt(tri, V; rng::AbstractRNG=Random.default_rng())

Triangulates the cavity, represented as a counter-clockwise list of 
vertices `V` with indices corresponding to those in `points`, 
left behind when deleting triangles intersected in a triangulation by an edge. 
If a triangulation is provided, the points are used from that.
"""
function triangulate_cavity_cdt(points, V;
    IntegerType::Type{I}=Int64,
    EdgeType::Type{E}=NTuple{2,IntegerType},
    TriangleType::Type{Vs}=NTuple{3,IntegerType},
    EdgesType::Type{Es}=Set{EdgeType},
    TrianglesType::Type{Ts}=Set{TriangleType},
    rng::AbstractRNG=Random.default_rng()) where {I,E,Vs,Es,Ts}
    tri = Triangulation(points;
        IntegerType,
        EdgeType,
        TriangleType,
        EdgesType,
        TrianglesType)
    triangulate_cavity_cdt!(tri, V; rng)
    return tri
end
function triangulate_cavity_cdt(tri::Triangulation{P,Ts,I,E,Es,BN,B,BIR}, V; rng::AbstractRNG=Random.default_rng()) where {P,Ts,I,E,Es,BN,B,BIR}
    return triangulate_cavity_cdt(get_points(tri), V;
        IntegerType=I,
        EdgeType=E,
        TriangleType=triangle_type(Ts),
        EdgesType=Es,
        TrianglesType=Ts,
        rng)
end

function setup_cavity_cdt(tri::Triangulation, V; rng::AbstractRNG=Random.default_rng())
    v = V[begin]
    u = V[end]
    prev, next, shuffled_indices = prepare_vertex_linked_list(V)
    delete_polygon_vertices_in_random_order!(tri, V, shuffled_indices, prev, next, u, v, rng)
    m = length(V)
    I = integer_type(tri)
    marked_vertices = nothing#I[]
    return prev, next, shuffled_indices, m, marked_vertices
end

function triangulate_cavity_cdt!(tri::Triangulation, V; rng::AbstractRNG=Random.default_rng())
    prev, next, shuffled_indices, m, marked_vertices = setup_cavity_cdt(tri, V; rng)
    add_triangle!(tri, V[begin], V[shuffled_indices[2]], V[end]; protect_boundary=true, update_ghost_edges=false)
    for i in 3:(m-1)
        a, b, c = index_shuffled_linked_list(V, next, prev, shuffled_indices, i)
        add_point_cavity_cdt!(tri, marked_vertices, a, b, c)
        #=
        if a ∈ marked_vertices
        end
        =#
    end
    return nothing
end

function add_point_cavity_cdt!(tri::Triangulation, marked_vertices, u, v, w)
    x = get_adjacent(tri, w, v)
    if !edge_exists(x)
        insert_flag = true
    else
        p, q, r, s = get_point(tri, w, v, x, u) # Don't want to deal with boundary handling here 
        incircle_test = point_position_relative_to_circle(p, q, r, s)
        orient_test = triangle_orientation(tri, u, v, w)
        insert_flag = !is_inside(incircle_test) && is_positively_oriented(orient_test)
    end
    if insert_flag
        add_triangle!(tri, u, v, w; protect_boundary=true, update_ghost_edges=false)
    else
        delete_triangle!(tri, w, v, x; protect_boundary=true, update_ghost_edges=false)
        add_point_cavity_cdt!(tri, marked_vertices, u, v, x)
        add_point_cavity_cdt!(tri, marked_vertices, u, x, w)
        #=
        if !is_inside(incircle_test)
            push!(marked_vertices, u, v, w, x)
        end
        =#
    end
    return nothing
end

function add_new_triangles!(tri_original::Triangulation, tri_lower, tri_upper)
    for tri in each_triangle(tri_lower)
        add_triangle!(tri_original, tri; protect_boundary=true, update_ghost_edges=false)
    end
    for tri in each_triangle(tri_upper)
        add_triangle!(tri_original, tri; protect_boundary=true, update_ghost_edges=false)
    end
    return nothing
end