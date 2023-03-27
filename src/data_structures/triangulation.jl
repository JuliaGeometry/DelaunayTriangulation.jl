"""
    Triangulation{P,Ts,I,E,Es,BN,B,BIR}

Struct for a triangulation.

See also [`triangulate`](@ref) and [`generate_mesh`](@ref).

# Fields 
- `points::P`

The nodes in the triangulation. Do note that this is not necessarily the same 
as all points that appear in a triangle, which could happen when the triangulation 
is incomplete or a point has been deleted. `get_vertices` may be useful if you 
are concerned about this, checking with `is_boundary_index` to ensure no 
boundary indices are involved.
- `triangles::Ts`

The triangles in the triangulation. All triangles are positively oriented.
- `adjacent::Adjacent{I, E}`

The [`Adjacent`](@ref) map, mapping edges to vertices that together form a 
positively oriented triangle in `triangles`. See also [`Adjacent`](@ref).
- `adjacent2vertex::Adjacent2Vertex{I,Es,E}`

The [`Adjacent2Vertex`](@ref) map, mapping vertices to all edges that together 
form a positively oriented triangle in `triangles`. See also [`Adjacent2Vertex`](@ref).
- `graph::Graph{I}`

The [`Graph`](@ref), mapping vertices to all other vertices that are connected with 
that vertex via an edge of a triangle in `triangles`. See also [`Graph`](@ref).
- `boundary_nodes::BN`

The boundary nodes, with the outer boundary given in counter-clockwise order
and the inner boundaries given in clockwise order. The default 
(default because you could customise it if you wish, see 
[`Interfaces`](@ref)) form of this field depends on your domain:

-- `Vector{Int64}`

If you provide no boundary curves or provide a single connected boundary curve, this 
will be a counter-clockwise vector of boundary nodes with `boundary_nodes[begin] == 
boundary_nodes[end]`. Note that if no boundary curve was provided, this gives the 
convex hull of `points`.

-- `Vector{Vector{Int64}}`

In this case, you have provided a single boundary curve but separated into different parts. 
In this case, `boundary_nodes[n]` gives the boundary nodes for the `n`th segment. Additionally,
`boundary_nodes[n][end] == boundary_nodes[n+1][begin]` and `boundary_nodes[end][end] ==
boundary_nodes[begin][begin]`. The nodes should be counter-clockwise in this case also.

-- `Vector{Vector{Vector{Int64}}}`

In this case, you have provided multiple boundary curves each separated into different parts, 
e.g. an annulus with each circle split into its lower and upper halves. The component 
`boundary_nodes[m][n]` gives the boundary nodes for the `n`th segment of the `m`th 
boundary curve, with `boundary_nodes[begin]` the outer-most boundary and `boundary_nodes[i]` the 
inner boundaries (`i > 1`). As in the previous case, `boundary_nodes[m][n][end] ==
boundary_nodes[m][n+1][begin]` and `boundary_nodes[m][end][end] == boundary_nodes[m][begin][begin]`. 
With this form, `boundary_nodes[m]` should be a counter-clockwise list of nodes for `m == 1`,
while for `m > 1` it` should be a clockwise list of nodes.

- `boundary_map::B`

This is an `OrderedDict` that maps a given boundary index to its position in `boundary_nodes`.
For example, if `boundary_map[-4] = (2, 3)`, this means that the boundary index `-4` 
corresponds to the nodes in `get_boundary_nodes(boundary_nodes, 2, 3)`. If there is just a 
single continuous curve, so that `boundary_nodes` acts like a vector of integers, then 
`boundary_map[$BoundaryIndex]` simply returns `boundary_nodes`. See also 
[`construct_boundary_map`](@ref) and [`map_boundary_index`](@ref). The ordering is such that 
the boundary index with the lowest absolute value comes first.
- `boundary_index_ranges::BIR`

This is an `OrderedDict` that maps a boundary index to a range of all other boundary indices 
that the corresponding boundary curve could correspond to. For example, for a curve with four 
segments, there are four possible boundary indices that a segment could correspond to - this 
`Dict` will extract this range. See also [`construct_boundary_index_ranges`](@ref).
- `constrained_edges::Es`

This is the set of extra edges added into the triangulation that you have provided. Note that 
this will not include any of the other constrained edges from `boundary_nodes`, for ths see 
the `all_constrained_edges` field. Moreover, you should include duplicate edges, i.e.
do not include both `(i, j)` and `(j, i)` - order is not important here.

If you have a constrained segment that happens to be collinear with another vertex, and that 
vertex is on the segment, then we mutate `constrained_edges` so that the segment is split 
at this vertex. 

- `convex_hull::ConvexHull{P,Vector{I}}`

This will be a vector of integers corresponding to indices in the points that 
together give the convex hull of the set of points, with `convex_hull[begin] == convex_hull[end]`.

- `all_constrained_edges::Es`

This is a set of all constrained edges, basically the union of the constrained edges in `constrained_edges`
and `boundary_nodes`. You shouldn't need to work with this field directly.

# Constructors 

There are several ways to construct this struct directly (though you should probably 
be using [`triangulate`](@ref) or [`generate_mesh`](@ref)).

## Default Constructor 

The default constructor is available, i.e. 

    Triangulation(points, triangles, adjacent, adjacent2vertex, graph, boundary_nodes, constrained_edges, boundary_map, convex_hull, all_constrained_edges)

## Empty Triangulation 

If you want to initialise the triangulation, you can use 

    Triangulation(points::P;
        IntegerType::Type{I}=Int64,
        EdgeType::Type{E}=NTuple{2,IntegerType},
        TriangleType::Type{V}=NTuple{3,IntegerType},
        EdgesType::Type{Es}=Set{EdgeType},
        TrianglesType::Type{Ts}=Set{TriangleType},
        boundary_nodes::BN=IntegerType[],
        constrained_edges=initialise_edges(EdgesType)) where {P,Ts,I,E,Es,BN,V}

## Existing Triangulation Objects 

The second constructor, mainly existing so that [`generate_mesh`](@ref) can convert 
its results into a `Triangulation`, takes in existing `triangles`, `points`, and 
`boundary_nodes` and creates a `Triangulation`:

    Triangulation(points::P, triangles::T, boundary_nodes::BN;
        IntegerType::Type{I}=Int64,
        EdgeType::Type{E}=NTuple{2,IntegerType},
        TriangleType::Type{V}=NTuple{3,IntegerType},
        EdgesType::Type{Es}=Set{EdgeType},
        TrianglesType::Type{Ts}=Set{TriangleType},
        add_ghost_triangles=false) where {I, E, V, Es, Ts, T}

You can set `add_ghost_triangles = true` to add ghost triangles into the structure. 

# Extended help 
There are many functions available for working with a `Triangulation`, some of which 
we list below. 

## Accessors:

- [`get_points`](@ref)
- [`get_triangles`](@ref)
- [`get_adjacent`](@ref)
- [`get_adjacent2vertex`](@ref)
- [`get_graph`](@ref)
- [`get_boundary_nodes`](@ref)
- [`get_boundary_map`](@ref)
- [`get_boundary_index_ranges`](@ref)
- [`get_constrained_edges`](@ref)
- [`get_convex_hull`](@ref)

## Operations:

- [`add_triangle!`](@ref)
- [`delete_triangle!`](@ref)
- [`add_ghost_triangles!`](@ref)
- [`delete_ghost_triangles!`](@ref)
- [`add_boundary_information!`](@ref)
- [`add_point!`](@ref)
- [`delete_point!`](@ref)
- [`flip_edge!`](@ref)
- [`legalise_edge!`](@ref)
- [`split_edge!`](@ref)
- [`split_triangle!`](@ref)

## Point location:

- [`brute_force_search`](@ref)
- [`jump_and_march`](@ref)

## Working with the [`Adjacent`](@ref) field:

- [`get_adjacent`](@ref)
- [`add_adjacent!`](@ref)
- [`delete_adjacent!`](@ref)
  
## Working with the [`Adjacent2Vertex`](@ref) field:

- [`get_adjacent2vertex`](@ref)
- [`add_adjacent2vertex!`](@ref)
- [`delete_adjacent2vertex!`](@ref)

## Working with the [`Graph`](@ref) field:

- [`get_edges`](@ref)
- [`get_vertices`](@ref)
- [`get_neighbours`](@ref)
- [`add_vertex!`](@ref)
- [`add_neighbour!`](@ref)
- [`delete_neighbour!`](@ref)
- [`delete_vertex!`](@ref)
- [`delete_boundary_vertices_from_graph!`](@ref)
- [`num_neighbours`](@ref)

## Working with the `boundary_nodes` field:

- [`has_multiple_curves`](@ref)
- [`has_multiple_segments`](@ref)
- [`num_curves`](@ref)
- [`num_segments`](@ref)
- [`get_boundary_nodes`](@ref)
- [`map_boundary_index`](@ref)
- [`get_curve_index`](@ref)
- [`get_segment_index`](@ref)
- [`num_outer_boundary_segments`](@ref)
- [`get_right_boundary_node`](@ref)
- [`get_left_boundary_node`](@ref)

## Working with the `convex_hull` field:

- [`get_convex_hull_indices`](@ref)
- [`convex_hull!`](@ref)

## Working with the `triangles` field:

- [`triangle_type`](@ref)
- [`num_triangles`](@ref)
- [`each_triangle`](@ref)
- [`each_solid_triangle`](@ref)
- [`each_ghost_triangle`](@ref)
- [`contains_triangle`](@ref)

## Working with edges:

- [`edge_type`](@ref)
- [`num_edges`](@ref)
- [`each_edge`](@ref)
- [`each_solid_edge`](@ref)
- [`each_ghost_edge`](@ref)

## Working with the `points` field:

- [`get_point`](@ref)
- [`each_point_index`](@ref)
- [`each_point`](@ref)
- [`num_points`](@ref)

## Working with the `all_constrained_edges` field:

- [`get_all_constrained_edges`](@ref)

## Predicates:

- [`is_boundary_edge`](@ref)
- [`is_boundary_triangle`](@ref)
- [`triangle_orientation`](@ref)
- [`point_position_relative_to_circumcircle`](@ref)
- [`point_position_relative_to_line`](@ref)
- [`point_position_on_line_segment`](@ref)
- [`line_segment_intersection_type`](@ref)
- [`point_position_relative_to_triangle`](@ref)
- [`is_outer_boundary_index`](@ref)
- [`is_outer_boundary_node`](@ref)
- [`is_boundary_node`](@ref)
- [`edge_exists`](@ref)
- [`has_ghost_triangles`](@ref)
- [`has_boundary_nodes`](@ref)
- [`is_legal`](@ref)

## Miscellaneous:

- [`integer_type`](@ref)
- [`number_type`](@ref)
- [`compute_representative_points!`](@ref)
- [`clear_empty_features!`](@ref)
- [`find_edge`](@ref)
- [`is_constrained`](@ref)
- [`all_boundary_indices`](@ref)
- [`get_surrounding_polygon`](@ref)
"""
struct Triangulation{P,Ts,I,E,Es,BN,B,BIR}
    points::P
    triangles::Ts
    adjacent::Adjacent{I,E}
    adjacent2vertex::Adjacent2Vertex{I,Es,E}
    graph::Graph{I}
    boundary_nodes::BN
    boundary_map::B
    boundary_index_ranges::BIR
    constrained_edges::Es
    all_constrained_edges::Es
    convex_hull::ConvexHull{P,Vector{I}}
end
function Base.show(io::IO, ::MIME"text/plain", tri::Triangulation)
    println(io, "Delaunay Triangulation.")
    println(io, "    Constrained: $(is_constrained(tri))")
    println(io, "    Has ghost triangles: $(has_ghost_triangles(tri))")
    println(io, "    Number of points: $(num_points(tri))")
    println(io, "    Number of triangles: $(num_triangles(tri))")
    print(io, "    Number of edges: $(num_edges(tri))")
end

Base.@constprop :aggressive function Triangulation(points::P;
    IntegerType::Type{I}=Int64,
    EdgeType::Type{E}=NTuple{2,IntegerType},
    TriangleType::Type{V}=NTuple{3,
        IntegerType},
    EdgesType::Type{Es}=Set{EdgeType},
    TrianglesType::Type{Ts}=Set{TriangleType},
    boundary_nodes::BN=IntegerType[],
    constrained_edges=initialise_edges(EdgesType)) where {P,
    Ts,
    I,
    E,
    Es,
    BN,
    V}
    T = initialise_triangles(Ts)
    adj = Adjacent{I,E}()
    adj2v = Adjacent2Vertex{I,Es,E}()
    graph = Graph{I}()
    bn_map = construct_boundary_map(boundary_nodes; IntegerType=I)
    bn_range = construct_boundary_index_ranges(boundary_nodes; IntegerType=I)
    ch = ConvexHull(points, I[])
    all_constrained_edges = merge_constrained_edges(bn_map, boundary_nodes, constrained_edges)
    n = num_points(points)
    sizehint!(T, 2n - 5) # maximum number of triangles
    sizehint!(adj, 3n - 6) # maximum number of edges 
    sizehint!(adj2v, n)
    sizehint!(graph, 3n - 6, n, n)
    sizehint!(ch, n)
    tri = Triangulation(points, T, adj, adj2v, graph, boundary_nodes, bn_map, bn_range,
        constrained_edges, all_constrained_edges, ch)
    return tri
end

function merge_constrained_edges(bn_map, boundary_nodes, constrained_edges::Es) where {Es}
    all_constrained = initialise_edges(Es)
    E = edge_type(Es)
    for segment_index in values(bn_map)
        bn_nodes = get_boundary_nodes(boundary_nodes, segment_index)
        nedges = num_boundary_edges(bn_nodes)
        for edge_idx in 1:nedges 
            vᵢ = get_boundary_nodes(bn_nodes, edge_idx)
            vᵢ₊₁ = get_boundary_nodes(bn_nodes,edge_idx+1)
            e = construct_edge(E, vᵢ, vᵢ₊₁)
            add_edge!(all_constrained, e)
        end
    end
    for e in each_edge(constrained_edges)
        add_edge!(all_constrained, e)
    end
    return all_constrained
end

## Accessors
for n in fieldnames(Triangulation)
    name = String(n)
    @eval begin
        @doc """
        get_$($(name))(tri::Triangulation)

    Returns the $($name) field from the triangulation `tri`.
    """ ($(Symbol("get_$n")))(tri::Triangulation) = tri.$n
    end
end

## Extending methods of adjacent, etc 
# Adjacent 
@inline function get_adjacent(tri::Triangulation, uv;
    check_existence::V=Val(has_multiple_segments(tri))) where {V}
    return get_adjacent(get_adjacent(tri), uv; check_existence,
        boundary_index_ranges=get_boundary_index_ranges(tri))
end
@inline function get_adjacent(tri::Triangulation, u, v;
    check_existence::V=Val(has_multiple_segments(tri))) where {V}
    return get_adjacent(get_adjacent(tri), u, v; check_existence,
        boundary_index_ranges=get_boundary_index_ranges(tri))
end
@inline add_adjacent!(tri::Triangulation, uv, w) = add_adjacent!(get_adjacent(tri), uv, w)
@inline function add_adjacent!(tri::Triangulation, u, v, w)
    return add_adjacent!(get_adjacent(tri), u, v, w)
end
@inline delete_adjacent!(tri::Triangulation, uv) = delete_adjacent!(get_adjacent(tri), uv)
@inline function delete_adjacent!(tri::Triangulation, u, v)
    return delete_adjacent!(get_adjacent(tri), u, v)
end

# Adjacent2Vertex
@inline function get_adjacent2vertex(tri::Triangulation, w)
    return get_adjacent2vertex(get_adjacent2vertex(tri), w)
end
@inline function add_adjacent2vertex!(tri::Triangulation, w, uv)
    return add_adjacent2vertex!(get_adjacent2vertex(tri), w, uv)
end
@inline function add_adjacent2vertex!(tri::Triangulation, w, u, v)
    return add_adjacent2vertex!(get_adjacent2vertex(tri), w, u, v)
end
@inline function delete_adjacent2vertex!(tri::Triangulation, w, uv)
    return delete_adjacent2vertex!(get_adjacent2vertex(tri), w, uv)
end
@inline function delete_adjacent2vertex!(tri::Triangulation, w, u, v)
    return delete_adjacent2vertex!(get_adjacent2vertex(tri), w, u, v)
end
@inline function delete_adjacent2vertex!(tri::Triangulation, w)
    return delete_adjacent2vertex!(get_adjacent2vertex(tri), w)
end

# Graph 
@inline get_edges(tri::Triangulation) = get_edges(get_graph(tri))
@inline get_vertices(tri::Triangulation) = get_vertices(get_graph(tri))
@inline get_neighbours(tri::Triangulation) = get_neighbours(get_graph(tri))
@inline get_neighbours(tri::Triangulation, u) = get_neighbours(get_graph(tri), u)
@inline add_vertex!(tri::Triangulation, u...) = add_vertex!(get_graph(tri), u...)
@inline num_neighbours(tri::Triangulation, u) = num_neighbours(get_graph(tri), u)
@inline function add_neighbour!(tri::Triangulation, u, v...)
    return add_neighbour!(get_graph(tri), u, v...)
end
@inline function delete_neighbour!(tri::Triangulation, u, v...)
    return delete_neighbour!(get_graph(tri), u, v...)
end
@inline delete_vertex!(tri::Triangulation, u...) = delete_vertex!(get_graph(tri), u...)
@inline function delete_boundary_vertices_from_graph!(tri::Triangulation)
    return delete_boundary_vertices_from_graph!(get_graph(tri))
end
@inline each_vertex(tri::Triangulation) = each_vertex(get_graph(tri))
@inline num_vertices(tri::Triangulation) = num_vertices(get_graph(tri))

# Boundary Nodes 
@inline function has_multiple_curves(tri::Triangulation)
    return has_multiple_curves(get_boundary_nodes(tri))
end
@inline function has_multiple_segments(tri::Triangulation)
    return has_multiple_segments(get_boundary_nodes(tri))
end
@inline num_curves(tri::Triangulation) = num_curves(get_boundary_nodes(tri))
@inline num_segments(tri::Triangulation) = num_segments(get_boundary_nodes(tri))
@inline function get_boundary_nodes(tri::Triangulation, mnℓ...)
    return get_boundary_nodes(get_boundary_nodes(tri), mnℓ...)
end
@inline function map_boundary_index(tri::Triangulation, i)
    return map_boundary_index(get_boundary_map(tri), i)
end
@inline get_curve_index(tri::Triangulation, i) = get_curve_index(get_boundary_map(tri), i)
@inline function get_segment_index(tri::Triangulation, i)
    return get_segment_index(get_boundary_map(tri), i)
end
@inline function num_outer_boundary_segments(tri::Triangulation)
    return num_outer_boundary_segments(get_boundary_nodes(tri))
end
@inline function get_right_boundary_node(tri::Triangulation, k, boundary_index)
    return get_right_boundary_node(get_adjacent(tri), k, boundary_index,
        get_boundary_index_ranges(tri),
        Val(has_multiple_segments(tri)))
end
@inline function get_left_boundary_node(tri::Triangulation, k, boundary_index)
    return get_left_boundary_node(get_adjacent(tri), k, boundary_index,
        get_boundary_index_ranges(tri),
        Val(has_multiple_segments(tri)))
end
@inline function get_boundary_index_range(tri::Triangulation, i)
    return map_boundary_index(get_boundary_index_ranges(tri), i)
end

# Convex Hull 
@inline get_convex_hull_indices(tri::Triangulation) = get_indices(get_convex_hull(tri))
function convex_hull!(tri::Triangulation; reconstruct=is_constrained(tri))
    I = integer_type(tri)
    if reconstruct
        convex_hull!(get_convex_hull(tri))
        return nothing
    elseif has_ghost_triangles(tri)
        outer_boundary_edges = get_adjacent2vertex(tri, I(BoundaryIndex)) # Not exactly all of them, e.g. the outer boundary could be made up of multiple segments, but this will be close
        indices = get_convex_hull_indices(tri)
        empty!(indices)
        sizehint!(indices, num_edges(outer_boundary_edges))
        e = first(each_edge(outer_boundary_edges))
        u = initial(e) # Need to start somewhere 
        push!(indices, u)
        start_index = u
        v = get_right_boundary_node(tri, u, I(BoundaryIndex))
        while v ≠ start_index
            push!(indices, v)
            v = get_right_boundary_node(tri, v, I(BoundaryIndex))
        end
        push!(indices, start_index)
        return nothing
    else
        add_ghost_triangles!(tri)
        convex_hull!(tri; reconstruct=false)
        delete_ghost_triangles!(tri)
    end
    return nothing
end

# Triangles 
@inline triangle_type(::Triangulation{P,Ts}) where {P,Ts} = triangle_type(Ts)
@inline num_triangles(tri::Triangulation) = num_triangles(get_triangles(tri))
@inline each_triangle(tri::Triangulation) = each_triangle(get_triangles(tri))
@inline contains_triangle(tri::Triangulation, T) = contains_triangle(T, get_triangles(tri))
@inline function contains_triangle(tri::Triangulation, i, j, k)
    return contains_triangle(i, j, k, get_triangles(tri))
end
@inline function construct_positively_oriented_triangle(tri::Triangulation, i, j, k)
    return construct_positively_oriented_triangle(triangle_type(tri), i, j, k,
        get_points(tri))
end

abstract type AbstractEachTriangle{T} end
Base.IteratorSize(::Type{<:AbstractEachTriangle}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:AbstractEachTriangle{T}}) where {T} = Base.IteratorEltype(T)
Base.eltype(::Type{<:AbstractEachTriangle{T}}) where {T} = triangle_type(T)
initialise_triangles(::Type{<:AbstractEachTriangle{T}}) where {T} = initialise_triangles(T)
each_triangle(itr::AbstractEachTriangle) = itr
struct EachSolidTriangle{T} <: AbstractEachTriangle{T}
    triangles::T
end
struct EachGhostTriangle{T} <: AbstractEachTriangle{T}
    triangles::T
end
each_solid_triangle(tri) = EachSolidTriangle(each_triangle(tri))
each_ghost_triangle(tri) = EachGhostTriangle(each_triangle(tri))
function Base.iterate(itr::EachSolidTriangle, state...)
    tri_state = iterate(itr.triangles, state...)
    tri_state === nothing && return nothing
    tri, state = tri_state
    while is_ghost_triangle(tri)
        tri_state = iterate(itr.triangles, state)
        tri_state === nothing && return nothing
        tri, state = tri_state
    end
    return tri, state
end
function Base.iterate(itr::EachGhostTriangle, state...)
    tri_state = iterate(itr.triangles, state...)
    tri_state === nothing && return nothing
    tri, state = tri_state
    while !is_ghost_triangle(tri)
        tri_state = iterate(itr.triangles, state)
        tri_state === nothing && return nothing
        tri, state = tri_state
    end
    return tri, state
end

# Edges 
@inline edge_type(::Triangulation{P,Ts,I,E}) where {P,Ts,I,E} = E
@inline num_edges(tri::Triangulation) = num_edges(get_graph(tri))
@inline each_edge(tri::Triangulation) = get_edges(get_graph(tri))

abstract type AbstractEachEdge{E} end
Base.IteratorSize(::Type{<:AbstractEachEdge}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:AbstractEachEdge{E}}) where {E} = Base.IteratorEltype(E)
Base.eltype(::Type{<:AbstractEachEdge{E}}) where {E} = edge_type(E)
initialise_edges(::Type{<:AbstractEachEdge{E}}) where {E} = initialise_edges(E)
each_edge(itr::AbstractEachEdge) = itr
struct EachSolidEdge{E} <: AbstractEachEdge{E}
    edges::E
end
struct EachGhostEdge{E} <: AbstractEachEdge{E}
    edges::E
end
each_solid_edge(tri) = EachSolidEdge(each_edge(tri))
each_ghost_edge(tri) = EachGhostEdge(each_edge(tri))
function Base.iterate(itr::EachSolidEdge, state...)
    edges_state = iterate(itr.edges, state...)
    edges_state === nothing && return nothing
    edges, state = edges_state
    while is_ghost_edge(edges)
        edges_state = iterate(itr.edges, state)
        edges_state === nothing && return nothing
        edges, state = edges_state
    end
    return edges, state
end
function Base.iterate(itr::EachGhostEdge, state...)
    edges_state = iterate(itr.edges, state...)
    edges_state === nothing && return nothing
    edges, state = edges_state
    while !is_ghost_edge(edges)
        edges_state = iterate(itr.edges, state)
        edges_state === nothing && return nothing
        edges, state = edges_state
    end
    return edges, state
end

# Points 
@inline function get_point(tri::Triangulation, i)
    return get_point(get_points(tri), get_boundary_map(tri), i)
end
@inline function get_point(tri::Triangulation, i::Vararg{Any,N}) where {N}
    return get_point(get_points(tri), get_boundary_map(tri), i...)
end
@inline each_point_index(tri::Triangulation) = each_point_index(get_points(tri))
@inline each_point(tri::Triangulation) = each_point(get_points(tri))
@inline num_points(tri::Triangulation) = num_points(get_points(tri))

abstract type AbstractEachVertex{V} end
Base.IteratorSize(::Type{<:AbstractEachVertex}) = Base.HasLength()
Base.IteratorEltype(::Type{<:AbstractEachVertex{V}}) where {V} = Base.IteratorEltype(V)
Base.eltype(::Type{<:AbstractEachVertex{V}}) where {V} = eltype(V)
each_vertex(itr::AbstractEachVertex) = itr
struct EachSolidVertex{V,T} <: AbstractEachVertex{V}
    vertices::V
    tri::T
end
Base.length(verts::EachSolidVertex) = num_vertices(verts.tri) - length(all_boundary_indices(verts.tri)) * (integer_type(verts.tri)(BoundaryIndex) ∈ verts.vertices)
struct EachGhostVertex{V,T} <: AbstractEachVertex{V}
    vertices::V
    tri::T
end
Base.length(verts::EachGhostVertex) = length(all_boundary_indices(verts.tri))
each_solid_vertex(tri) = EachSolidVertex(each_vertex(tri), tri)
each_ghost_vertex(tri) = EachGhostVertex(all_boundary_indices(tri), tri)
function Base.iterate(itr::EachSolidVertex, state...)
    vertices_state = iterate(itr.vertices, state...)
    vertices_state === nothing && return nothing
    vertices, state = vertices_state
    while is_boundary_index(vertices)
        vertices_state = iterate(itr.vertices, state...)
        vertices_state === nothing && return nothing
        vertices, state = vertices_state
    end
    return vertices, state
end
Base.iterate(itr::EachGhostVertex, state...) = Base.iterate(itr.vertices, state...)

# Predicates 
@inline is_boundary_edge(tri::Triangulation, ij) = is_boundary_edge(ij, get_adjacent(tri))
@inline function is_boundary_edge(tri::Triangulation, i, j)
    return is_boundary_edge(i, j, get_adjacent(tri))
end
@inline function is_boundary_triangle(tri::Triangulation, i, j, k)
    return is_boundary_triangle(i, j, k, get_adjacent(tri))
end
@inline function is_boundary_triangle(tri::Triangulation, T)
    return is_boundary_triangle(tri, geti(T), getj(T), getk(T))
end
@inline function triangle_orientation(tri::Triangulation, i, j, k)
    return triangle_orientation(i, j, k, get_points(tri), get_boundary_map(tri))
end
@inline function triangle_orientation(tri::Triangulation, T)
    return triangle_orientation(tri, geti(T), getj(T), getk(T))
end
@inline function point_position_relative_to_circumcircle(tri::Triangulation, i, j, k, ℓ)
    return point_position_relative_to_circumcircle(i, j, k, ℓ, get_points(tri),
        get_boundary_map(tri))
end
@inline function point_position_relative_to_circumcircle(tri::Triangulation, T, ℓ)
    return point_position_relative_to_circumcircle(tri, geti(T), getj(T), getk(T), ℓ)
end
@inline function point_position_relative_to_line(tri::Triangulation, i, j, u)
    return point_position_relative_to_line(i, j, u, get_points(tri), get_boundary_map(tri))
end
@inline function point_position_on_line_segment(tri::Triangulation, i, j, u)
    return point_position_on_line_segment(i, j, u, get_points(tri))
end
@inline function line_segment_intersection_type(tri::Triangulation, u, v, i, j)
    return line_segment_intersection_type(u, v, i, j, get_points(tri))
end
@inline function point_position_relative_to_triangle(tri::Triangulation, i, j, k, u)
    return point_position_relative_to_triangle(i, j, k, u, get_points(tri),
        get_boundary_map(tri))
end
@inline function point_position_relative_to_triangle(tri::Triangulation, T, u)
    return point_position_relative_to_triangle(tri, geti(T), getj(T), getk(T), u)
end
@inline function is_outer_boundary_index(tri::Triangulation, i)
    return is_outer_boundary_index(i, get_boundary_map(tri))
end
@inline function is_outer_boundary_node(tri::Triangulation, i)
    return is_outer_boundary_node(i, get_graph(tri), get_boundary_index_ranges(tri))
end
@inline function is_boundary_node(tri::Triangulation, i)
    return is_boundary_node(i, get_graph(tri), get_boundary_index_ranges(tri))
end
@inline edge_exists(tri::Triangulation, i, j) = edge_exists(i, j, get_adjacent(tri))
@inline edge_exists(tri::Triangulation, ij) = edge_exists(ij, get_adjacent(tri))
@inline function has_ghost_triangles(tri::Triangulation)
    return has_ghost_triangles(get_adjacent(tri), get_adjacent2vertex(tri))
end
@inline function has_boundary_nodes(tri::Triangulation)
    return has_multiple_segments(tri) || num_boundary_edges(get_boundary_nodes(tri)) ≠ 0
end
function is_legal(tri::Triangulation, i, j)
    (is_boundary_edge(tri, i, j) ||
     is_boundary_edge(tri, j, i) ||
     !edge_exists(tri, i, j) ||
     !edge_exists(tri, j, i)) && return Cert.Legal
    k = get_adjacent(tri, i, j)
    ℓ = get_adjacent(tri, j, i)
    p, q, r, s = get_point(tri, i, j, k, ℓ)
    cert = is_legal(p, q, r, s)
    return cert
end

## Point Location
function brute_force_search(tri::Triangulation, q)
    return brute_force_search(get_triangles(tri), q, get_points(tri), get_boundary_map(tri))
end
function jump_and_march(tri::Triangulation, q;
    point_indices=each_point_index(tri),
    m=default_num_samples(length(point_indices)),
    try_points=(),
    k=select_initial_point(get_points(tri), q; m, point_indices,
        try_points),
    check_existence::C=Val(has_multiple_segments(tri)),
    rng::AbstractRNG=Random.default_rng()) where {C}
    return jump_and_march(get_points(tri),
        get_adjacent(tri),
        get_adjacent2vertex(tri),
        get_graph(tri),
        get_boundary_index_ranges(tri),
        get_boundary_map(tri),
        q; m, point_indices, try_points, k,
        TriangleType=triangle_type(tri), check_existence, rng)
end

# Miscellaneous
@inline integer_type(::Triangulation{P,Ts,I}) where {P,Ts,I} = I
@inline number_type(::Triangulation{P}) where {P} = number_type(P)
@inline function compute_representative_points!(tri::Triangulation;
    use_convex_hull=!has_multiple_segments(tri) &&
                    num_boundary_edges(get_boundary_nodes(tri)) ==
                    0)
    if !use_convex_hull
        compute_representative_points!(get_points(tri), get_boundary_nodes(tri))
    else
        compute_representative_points!(get_points(tri), get_convex_hull_indices(tri))
    end
end
function clear_empty_features!(tri)
    adj = get_adjacent(tri)
    adj2v = get_adjacent2vertex(tri)
    graph = get_graph(tri)
    clear_empty_keys!(adj)
    clear_empty_keys!(adj2v)
    clear_empty_points!(graph)
    return nothing
end
@inline find_edge(tri::Triangulation, T, ℓ) = find_edge(T, get_points(tri), ℓ)
function is_constrained(tri::Triangulation)
    return !is_empty(get_constrained_edges(tri)) || has_boundary_nodes(tri)
end
all_boundary_indices(tri::Triangulation) = keys(get_boundary_index_ranges(tri))
get_surrounding_polygon(tri::Triangulation, u; skip_boundary_indices=false) = get_surrounding_polygon(get_adjacent(tri), get_graph(tri), u, get_boundary_index_ranges(tri), Val(has_multiple_segments(tri)); skip_boundary_indices)

## Triangulating points, triangles, boundary_nodes
function Triangulation(points::P, triangles::T, boundary_nodes::BN;
    IntegerType::Type{I}=Int64,
    EdgeType::Type{E}=NTuple{2,IntegerType},
    TriangleType::Type{V}=NTuple{3,IntegerType},
    EdgesType::Type{Es}=Set{EdgeType},
    TrianglesType::Type{Ts}=Set{TriangleType},
    add_ghost_triangles=false) where {P,BN,I,E,V,Es,Ts,T}
    tri = Triangulation(points; boundary_nodes, IntegerType, EdgeType, TriangleType,
        EdgesType, TrianglesType)
    adj = get_adjacent(tri)
    adj2v = get_adjacent2vertex(tri)
    graph = get_graph(tri)
    tris = get_triangles(tri)
    for τ in each_triangle(triangles)
        # add_triangle!(tri, τ)  # Not using this so that we can avoid boundary handling 
        add_triangle!(adj, τ)
        add_triangle!(adj2v, τ)
        add_triangle!(graph, τ)
        add_triangle!(tris, τ)
    end
    add_boundary_information!(tri)
    add_ghost_triangles && add_ghost_triangles!(tri)
    convex_hull!(tri; reconstruct=true)
    return tri
end

## Define all the docstrings 
@doc """
    get_adjacent(tri::Triangulation, uv; check_existence::V=Val(false)) where {V}
    get_adjacent(tri::Triangulation, u, v; check_existence::V=Val(false)) where {V}

Given the edge `(u, v)`, gets the vertex `w` such that `(u, v, w)`
is a positively oriented triangle in the triangulation `tri`.

Use `check_existence=Val(true)` to be safe against ghost edges corresponding to 
neighbouring boundary indices.
""" get_adjacent(::Triangulation, ::Any)

@doc """
    add_adjacent!(tri::Triangulation, uv, w)
    add_adjacent!(tri::Triangulation, u, v, w)

Given the edge `(u, v)`, adds the vertex `w` into the adjacent map 
of the triangulation `tri`, so that `(u, v, w)` is a positively 
oriented triangle in `tri`.
""" add_adjacent!(::Triangulation, ::Any)

@doc """
    delete_adjacent!(tri::Triangulation, uv)
    delete_adjacent!(tri::Triangulation, u, v)

Given the edge `(u, v)`, deletes the key `(u, v)` from the adjacent map 
of the triangulation `tri`.
""" delete_adjacent!(::Triangulation, ::Any)

@doc """
    get_adjacent2vertex(tri::Triangulation, w) 

Given the vertex `w`, returns the set of edges `(u, v)` from 
the adjacent2vertex map such that `(u, v, w)` is a positively oriented 
triangle in `tri`.
""" get_adjacent2vertex(::Triangulation, ::Any)

@doc """
    add_adjacent2vertex!(tri::Triangulation, w, uv)
    add_adjacent2vertex!(tri::Triangulation, w, u, v)

Given the vertex `w` and an edge `(u, v)`, pushes the edge `(u, v)` 
into the set of edges defined by the key `w` in the adjacent2vertex map 
of the triangulation `tri`.
""" add_adjacent2vertex!(::Triangulation, ::Any, ::Any)

@doc """
    delete_adjacent2vertex!(tri::Triangulation, w, uv)
    delete_adjacent2vertex!(tri::Triangulation, w, u, v)

Given the vertex `w` and an edge `(u, v)`, deletes the edge `(u, v)`
from the set of edges defined by the key `w` in the adjacent2vertex map 
of the triangulation `tri`.
""" delete_adjacent2vertex!(::Triangulation, ::Any, ::Any)

@doc """
    delete_adjacent2vertex!(tri::Triangulation, w)

Given the vertex `w`, deletes the set of edges `(u, v)` that define positively 
oriented triangles `(u, v, w)` from the adjacent2vertex map of the triangulation 
`tri`.
""" delete_adjacent2vertex!(::Triangulation, ::Any)

@doc """
    get_edges(tri::Triangulation)

Given a triangulation `tri`, returns the set of edges in `tri`. These edges will 
not be oriented, i.e. if `(i, j)` is in the set then `(j, i)` will not be.
""" get_edges(::Triangulation)

@doc """
    get_vertices(tri::Triangulation)

Given a triangulation `tri`, returns the set of vertices in `tri`. Note that this 
will include any ghost vertices.
""" get_vertices(::Triangulation)

@doc """
    get_neighbours(tri::Triangulation)

Given a triangulation `tri`, returns the set of neighbourhoods.
""" get_neighbours(::Triangulation)

@doc """
    get_neighbours(tri::Triangulation, u)

Given a vertex `u` and a triangulation `tri`, returns the set of vertices `v`
such that `(u, v)` is an edge of the triangulation for each `v` in the set.
""" get_neighbours(::Triangulation, ::Any)

@doc """
    add_vertex!(tri::Triangulation, u...)

Given vertices `u...` and a triangulation `tri`, adds all the vertices `u...` 
into the graph `tri.graph`. 

Note that this does not insert the vertex into the triangulation - 
see [`add_point!`](@ref) for this.
""" add_vertex!(::Triangulation, ::Vararg)

@doc """
    num_neighbours(tri::Triangulation, u)

Given a triangulation `tri` and a vertex `u`, returns the number of neighbours 
of `u`.
""" num_neighbours(::Triangulation, ::Any)

@doc """
    add_neighbour!(tri::Triangulation, u, v...)

Given a vertex `u` and other vertices `v...`, adds all the vertices `v...` 
into the neighbourhood of `u` from the triangulation `tri`, i.e. into `tri.graph`.
""" add_neighbour!(::Triangulation, ::Any, ::Vararg)

@doc """
    delete_neighbour!(tri::Triangulation, u, v...)

Given a vertex `u` and other vertices `v...`, deletes all the vertices `v...`
from the neighbourhood of `u` from the triangulation `tri`, i.e. from `tri.graph`.
""" delete_neighbour!(::Triangulation, ::Any, ::Vararg)

@doc """
    delete_vertex!(tri::Triangulation, u...)

Given vertices `u...` and a triangulation `tri`, deletes all the vertices 
`u...` from the graph `tri.graph`.

Note that this does not remove the vertex from the triangulation - 
see [`remove_point!`](@ref) for this.
""" delete_vertex!(tri::Triangulation, ::Vararg)

@doc """
    delete_boundary_vertices_from_graph!(tri::Triangulation)

Given a triangulation `tri`, removes all boundary vertices from the 
triangulation `tri`, i.e. all those with index less than $BoundaryIndex (
see [`is_boundary_index`](@ref)).
""" delete_boundary_vertices_from_graph!(::Triangulation)

@doc """
    each_vertex(tri::Triangulation)

Given a triangulation `tri`, returns an iterator over each vertex in the triangulation. Note that 
these are indices of points rather than the points themselves.

See also [`each_solid_vertex`](@ref) and [`each_ghost_vertex`](@ref).
""" each_vertex(::Triangulation)

@doc """
    num_vertices(tri::Triangulation)

Returns the number of vertices in the triangulation `tri`.
""" num_vertices(::Triangulation)

@doc """
    has_multiple_curves(tri::Triangulation)

Returns `true` if the triangulation `tri` has multiple boundary curves, 
and `false` otherwise.
""" has_multiple_curves(::Triangulation)

@doc """
    has_multiple_segments(tri::Triangulation)

Returns `true` if the boundary curves of the triangulation `tri` 
are broken into multiple segments. Note that if `has_multiple_curves(tri)`,
then `has_multiple_segments(tri) == true`.
""" has_multiple_segments(::Triangulation)

@doc """
    num_curves(tri::Triangulation)

Returns the number of boundary curves of the triangulation `tri`.
""" num_curves(::Triangulation)

@doc """
    num_segments(tri::Triangulation)

Returns the number of segments of the boundary curve of the 
triangulation `tri`. This is only defined if `!has_multiple_curves(tri)`
and `has_multiple_segments(tri)`.
""" num_segments(::Triangulation)

@doc """
    get_boundary_nodes(tri::Triangulation, mnℓ...) 

Given a triangulation `tri`, returns the nodes corresponding to the 
indices in `mnℓ...`. For example, if `tri` has multiple boundary curves, 
then `get_boundary_nodes(tri, m)` returns the set of nodes for the 
`m`th curve, and `get_boundary_nodes(tri, m, n)` the set of nodes for the 
`n`th segment of the `m`th curve.
""" get_boundary_nodes(::Triangulation, ::Vararg)

@doc """
    map_boundary_index(tri::Triangulation, i)

Given a triangulation `tri`, returns the position in `tri.boundary_nodes`
corresponding to the boundary index `i`, i.e. returns 
`tri.boundary_map[i]`.
""" map_boundary_index(::Triangulation, ::Any)

@doc """
    get_curve_index(tri::Triangulation, i)

Given a triangulation `tri` and a boundary index `i`, returns the curve number 
that `i` corresponds to. For example, if `map_boundary_index(tri, i) = (m, n)`,
the curve index would be `m`. If `!has_multiple_curves(tri)`, this just returns 
`1`.
""" get_curve_index(::Triangulation, ::Any)

@doc """
    get_segment_index(tri::Triangulation, i)

Given a triangulation `tri` and a boundary index `i`, returns the segment number 
that `i` corresponds to. For example, if `map_boundary_index(tri, i) = (m, n)`, 
the segment index would be `m`. If `!has_multiple_segments(tri)`, this just returns 
`1`.
""" get_segment_index(::Triangulation, ::Any)

@doc """
    get_right_boundary_node(tri::Triangulation, k, boundary_index)

Given a triangulation `tri`, a boundary node `k`, and a `boundary_index` for the associated 
boundary, returns the index of the point to the right of `k` on the outer boundary.

See also [`get_left_boundary_node`](@ref).
""" get_right_boundary_node(::Triangulation, ::Any, ::Any)

@doc """
    get_left_boundary_node(tri::Triangulation, k)

Given a triangulation `tri` and an outer boundary node `k`, returns the index 
of the point to the left of `k` on the outer boundary.

See also [`get_right_boundary_node`](@ref).
""" get_left_boundary_node(::Triangulation, ::Any, ::Any)

@doc """
    get_convex_hull_indices(tri::Triangulation)

Returns indices for the [`ConvexHull`](@ref) of the triangulation `tri`.
""" get_convex_hull_indices(::Triangulation)

@doc """
    convex_hull!(tri::Triangulation; reconstruct = is_constrained(tri))

Computess the convex hull for the points in the triangulation `tri`, updating 
the field `tri.convex_hull`. Note that this will only be needed if you are e.g. 
constructing a convex hull from an existing set of points, triangles, and edges. 
If you already have an unconstrained Delaunay triangulation, this field will 
have also been updated. 

If `!reconstruct`, the convex hull is reconstructed using [`convex_hull`](@ref), 
otherwise the ghost triangles are used to construct it (if there are no ghost 
triangles, they are added and then removed afterwards).

See also [`convex_hull`](@ref) and [`ConvexHull`](@ref).
""" convex_hull!(::Triangulation)

@doc """
    triangle_type(tri::Triangulation{P,Ts}) where {P,Ts} 

Given a triangulation `tri`, returns the type used for 
representing the triangles.
""" triangle_type(::Triangulation)

@doc """
    num_triangles(tri::Triangulation)

Given a triangulation `tri`, returns the number of triangles.
""" num_triangles(::Triangulation)

@doc """
    each_triangle(tri::Triangulation)

Given a triangulation `tri`, returns an iterator over the triangles. This iterator 
will return both the ghost and solid triangles (if there are no ghost triangles, 
then those of course are not included and `each_triangle(tri) == each_solid_triangle(tri)`).

See also [`each_ghost_triangle`](@ref) and [`each_solid_triangle`](@ref).
""" each_triangle(::Triangulation)

@doc """
    each_solid_triangle(tri::Triangulation)

Given a triangulation `tri`, returns an iterator over the solid (i.e. non-ghost) triangles.

See also [`each_ghost_triangle`](@ref) and [`each_triangle`](@ref).
""" each_solid_triangle(::Triangulation)

@doc """
    each_ghost_triangle(tri::Triangulation)

Given a triangulation `tri`, returns an iterator over the ghost triangles.

See also [`each_triangle`](@ref) and [`each_solid_triangle`](@ref).
""" each_ghost_triangle(::Triangulation)

@doc """
    contains_triangle(tri::Triangulation, T) = contains_triangle(T, get_triangles(tri))
    contains_triangle(tri::Triangulation, i, j, k) = contains_triangle(i, j, k, get_triangles(tri))

Given a triangulation `tri`, tests if the triangle `T = (i, j, k)` is in the triangulation (at least, in `triangles`). 
The returned value is a tuple, with the first result being the form of `T` that is actually in 
the triangulation (e.g., if `T = (i, j, k)` but the form of `T` in `tri` is `(k, i, j)`, `(k, i, j)` is 
returned). The second argument is a `Bool`, with `true` if `T` is contained in `tri` 
and `false` otherwise. If `false`, the first output is `T`.
""" contains_triangle(::Triangulation, ::Any)

@doc """
    construct_positively_oriented_triangle(tri::Triangulation, i, j, k)

Given a triangulation `tri` and indices `i, j, k` for a proposed triangle, returns a triangle 
with the indices ordered such that `(i, j, k)` is not negatively oriented.
""" construct_positively_oriented_triangle(::Triangulation, ::Any, ::Any, ::Any)

@doc """
    edge_type(tri::Triangulation{P,Ts,I,E}) where {P,Ts,I,E}

Given a triangulation `tri`, returns the type used for 
representing individual edges.
""" edge_type(::Triangulation)

@doc """
    num_edges(tri::Triangulation)

Given a triangulation `tri`, returns the number of edges. Note 
that this does not double count edges, e.g. `(i, j)` and `(j, i)`
do not count as two edges.
""" num_edges(::Triangulation)

@doc """
    each_edge(tri::Triangulation) 

Given a triangulation `tri`, returns an iterator over the edges 
in the triangulation. Note that e.g. `(i, j)` and `(j, i)` will not 
be iterated over twice, only one of them would be.

See also [`each_solid_edge`](@ref) and [`each_ghost_edge`](@ref).
""" each_edge(::Triangulation)

@doc """
    each_solid_edge(tri::Triangulation)

Given a triangulation `tri`, returns an iterator over the solid (i.e. non-ghost) edges.

See also [`each_ghost_edge`](@ref) and [`each_edge`](@ref).
""" each_solid_edge(::Triangulation)

@doc """
    each_ghost_edge(tri::Triangulation)

Given a triangulation `tri`, returns an iterator over the ghost edges.

See also [`each_edge`](@ref) and [`each_solid_edge`](@ref).
""" each_ghost_edge(::Triangulation)

@doc """
    get_point(tri::Triangulation, i...)

Given a triangulation `tri` and some indices `i...`, returns 
the `i`th point from the triangulation. If `i` is just a single integer, 
the result is a `Tuple` of the coordinates, but if there are multiple 
integers provided then a `Tuple` of `Tuples` of these coordinates 
is returned, with the `j`th `Tuple` the coordinates for the `i[j]`th 
point. 

It is assumed that whenever `i` is not an integer, `i` is meant to be 
a point, so `(getx(i), gety(i))` would be returned in that case. This 
makes it easier to use some predicates without having to know the index 
of the point, simply passing the point directly.

If `is_boundary_index(i)`, then instead of returning the `i`th point of the
triangulation, the centroid for the boundary curve corresponding to the 
boundary index `i` is returned. See also [`get_representative_point_coordinates`](@ref).
""" get_point(::Triangulation, ::Any)

@doc """
    each_point_index(tri::Triangulation)

Given a triangulation `tri`, returns an iterator over the indices 
of the points in `tri`. This iterator does not include the boundary 
indices. Note also that this may include points already in the triangulation.

See also [`each_vertex`](@ref), [`each_solid_vertex`](@ref), and [`each_ghost_vertex`](@ref).
""" each_point_index(::Triangulation)

@doc """
    each_point(tri::Triangulation)

Given a triangulation `tri`, returns an iterator over the coordinates 
of the points in `tri`, whether they are present in `tri` already or not. 
This iterator does not include the centroids of the boundary curves.

See also [`each_vertex`](@ref), [`each_solid_vertex`](@ref), and [`each_ghost_vertex`](@ref).
""" each_point(::Triangulation)

@doc """
    num_points(tri::Triangulation)

Given a triangulation `tri`, returns the number of points. Note that this may include 
points not already in the triangulation.

See also [`num_vertices`](@ref).
""" num_points(::Triangulation)

@doc """
    each_solid_vertex(tri::Triangulation)

Given a triangulation `tri`, returns an iterator over the solid (i.e. non-ghost) vertices.

See also [`each_ghost_vertex`](@ref) and [`each_vertex`](@ref).
""" each_solid_vertex(::Triangulation)

@doc """
    each_ghost_vertex(tri::Triangulation)

Given a triangulation `tri`, returns an iterator over the ghost vertices.

See also [`each_vertex`](@ref) and [`each_solid_vertex`](@ref).
""" each_ghost_vertex(::Triangulation)

@doc """
    is_boundary_edge(tri::Triangulation, ij)
    is_boundary_edge(tri::Triangulation, i, j) 

Given a triangulation `tri` and an edge `(i, j)`, returns `true` 
if `(i, j)` is a boundary edge or `false` if not. 

This predicate takes care for the orientation of the edge: a boundary edge 
is one that has only two solid vretices, so that `get_adjacent(tri, ij) ≤ I(BoundaryIndex)`.
""" is_boundary_edge(::Triangulation, ::Any)

@doc """
    is_boundary_triangle(tri::Triangulation, i, j, k)
    is_boundary_triangle(tri::Triangulation, T)

Given a triangulation `tri` and a triangle `T`, returns `true` 
if `T` is a boundary triangle or `false` if not.

Note that this is different from testing if `T` is a ghost 
triangle - it is assumed that `T` is a solid triangle.   
""" is_boundary_triangle(::Triangulation, ::Any)

@doc """
    triangle_orientation(tri::Triangulation, i, j, k)
    triangle_orientation(tri::Triangulation, T)

Given a triangle `T = (i, j, k)` and a triangulation `tri`, 
computes the orientation of the triangle. In particular, returns:
    
- `Certificate.PositivelyOriented`: The triangle is positively oriented.
- `Certificate.Degenerate`: The triangle is degenerate, meaning the coordinates are collinear. 
- `Certificate.NegativelyOriented`: The triangle is negatively oriented.
    
A test is also made for the case that `is_ghost_triangle(T)`: If  `T` 
is a ghost triangle, then the index corresponding to a boundary index 
points to a centroid, in which case one of the edges has its orientation 
flipped. This case will also be handled correctly. In case the boundary 
index corresponds to an interior curve, this flip is not necessary.
""" triangle_orientation(::Triangulation, ::Any)

@doc """
    point_position_relative_to_circumcircle(tri::Triangulation, i, j, k, ℓ)
    point_position_relative_to_circumcircle(tri::Triangulation, T, ℓ)

Given a triangle `T = (i, j, k)`, a triangulation `tri`, and a point 
with index `ℓ`, `pₗ`, computes the position this point `pₗ` relative to the 
triangle's circumcircle, making use of [`point_position_relative_to_circle`](@ref).

The returned values are:

- `Certificate.Outside`: `pₗ` is outside of the circumcircle.
- `Certificate.On`: `pₗ` is on the circumcircle.
- `Certificate.Inside`: `pₗ` is inside the circumcircle.

A test is also made for the case that `is_ghost_triangle(T)`: When `T`
is a ghost triangle, one of its indices is a boundary index, say `i`. 
Since this vertex is treated as being out at infinity, the circumcircle 
degenerates into the line through the other two vertices. Thus, we test 
that `pₗ` is inside this circumcircle by seeing if it is in the oriented 
outer halfplane defined by the two other vertices. See also 
[`point_position_relative_to_oriented_outer_halfplane`](@ref).
""" point_position_relative_to_circumcircle(::Triangulation, ::Any, ::Any)

@doc """
    point_position_relative_to_line(tri::Triangulation, i, j, u) 

Given indices `i`, `j`, and `u` corresponding to indices 
of points in the triangulation `tri`, corresponding to 
coordinates say `a`, `b`, and `p`, respectively, computes the 
position of `p` relative to the oriented line `(a, b)`. 

The returned values are:

- `Certificate.Left`: `p` is to the left of the line. 
- `Certificate.Collinear`: `p` is on the line.
- `Certificate.Right`: `p` is to the right of the line.

If `is_ghost_edge(i, j)`, the oriented line `(a, b)` is flipped 
since the point corresponding to the boundary index will be a 
centroid which swaps the orientation.
""" point_position_relative_to_line(::Triangulation, ::Any, ::Any, ::Any)

@doc """
    point_position_on_line_segment(tri::Triangulation, i, j, u) 

Given indices `i`, `j`, and `u` corresponding to indices 
of points in the triangulation `tri`, corresponding to 
coordinates say `a`, `b`, and `p`, respectively, computes the 
position of `p` relative to the oriented line segment `(a, b)`,
assuming that the three points are collinear.
    
The returned values are:
    
- `Certificate.On`: `p` is on the line segment, meaning between `a` and `b`.
- `Certificate.Degenerate`: Either `p == a` or `p == b`, i.e. `p` is one of the endpoints. 
- `Certificate.Left`: `p` is off and to the left of the line segment.
- `Certificate.Right`: `p` is off and to the right of the line segment.
""" point_position_on_line_segment(::Triangulation, ::Any, ::Any, ::Any)

@doc """
    line_segment_intersection_type(tri::Triangulation, u, v, i, j) 

Given two pairs of indices `(u, v)` and `(i, j)`, with 
all indices corresponding to points in the triangulation `tri`, tests 
the number of intersections between the two line segments 
associated with these indices.
    
Letting `p`, `q`, `a`, and `b` be the points referred to by 
`u`, `v`, `i`, and `j`, respectively, we return:
    
- `Certificate.None`: The line segments do not meet at any points. 
- `Certificate.Multiple`: The closed line segments `[p, q]` and `[a, b]` meet in one or several points. 
- `Certificate.Single`: The open line segments `(p, q)` and `(a, b)` meet in a single point.
- `Certificate.On`: One of the endpoints is on `[a, b]`, but there are no other intersections.
""" line_segment_intersection_type(::Triangulation, ::Any, ::Any, ::Any, ::Any)

@doc """
    point_position_relative_to_triangle(tri::Triangulation, i, j, k, u)
    point_position_relative_to_triangle(tri::Triangulation, T, u)

Given a triangle `T = (i, j, k)` and another index `u`, 
with indices referring to points in the triangulation `tri`, computes 
the position of `u` relative to the triangle.
    
Letting `p`, `q`, `a`, and `b` be the points referred to by 
`i`, `j`, `k`, and `u`, respectively, we return:
    
- `Certificate.Outside`: `p` is outside of the triangle. 
- `Certificate.On`: `p` is on one of the edges. 
- `Certificate.Inside`: `p` is inside the triangle.
    
If `T` is a ghost triangle, `Certificate.Inside` is returned.
""" point_position_relative_to_triangle(::Triangulation, ::Any, ::Any)

@doc """
    is_outer_boundary_index(tri::Triangulation, i)

Given a triangulation `tri` and an index `i`, tests if `i` is the index of a boundary 
corresponding to the outermost boundary of the triangulation.
""" is_outer_boundary_index(tri::Triangulation, ::Any)

@doc """
    is_outer_boundary_node(tri::Triangulation, i)

Given a triangulation `tri` and an index `i`, tests if `i` is the index of a node 
on the outermost boundary of the triangulation.
""" is_outer_boundary_node(tri::Triangulation, ::Any)

@doc """
    is_boundary_node(tri::Triangulation, i)

Given a triangulation `tri` and an index `i`, tests if `i` is the index of a node 
on any part of the boundary of the triangulation. The returned value is a 
`Tuple`, whose first element is the result of this test, and whose second element 
is the corresponding boundary index if the test passed and `DefaultAdjacentValue` 
otherwise.
""" is_boundary_node(tri::Triangulation, ::Any)

@doc """
    edge_exists(tri::Triangulation, i, j)
    edge_exists(tri::Triangulation, ij)

Given a triangulation `tri` and an edge `(i, j)`, tests if the edge `(i, j)` is 
in the triangulation. 
""" edge_exists(::Triangulation, ::Any)

@doc """
    has_ghost_triangles(tri::Triangulation)

Given a triangulation `tri`, tests if the triangulation has ghost triangles.
""" has_ghost_triangles(::Triangulation)

@doc """
    has_boundary_nodes(tri::Triangulation)

Given a triangulation `tri`, tests if the triangulation has boundary nodes - 
these are nodes that are constrained to be there, not those on the convex hull.
""" has_boundary_nodes(::Triangulation)

@doc """

    is_legal(tri::Triangulation, i, j)

Given a triangulation `tri` and an edge `(i, j)` in `tri`, tests if the edge is legal,
returning -`Cetificate.Legal` if so and `Certificate.Illegal` otherwise.
"""

@doc """
    brute_force_search(tri::Triangulation, q)

Given a point `q`, finds the triangle in the triangulation `tri` containing it by 
searching over all triangles.
""" brute_force_search(::Triangulation, ::Any)

@doc """
    jump_and_march(tri::Triangulation, q;
        point_indices=each_point_index(tri),
        m=default_num_samples(length(point_indices)),
        try_points=(),
        k=select_initial_point(get_points(tri), q; m, point_indices, try_points),
        check_existence::C=Val(has_multiple_segments(tri)),
        rng::AbstractRNG=Random.default_rng())

Using the jump and march algorithm, finds the triangle in the triangulation `tri` containing the 
query point `q`.

# Arguments 
- `tri`: The [`Triangulation`](@ref).
- `q`: The query point.

# Keyword Arguments 
- `m=default_num_samples(num_points(pts))`: The number of samples to use when sampling an initial point from [`select_initial_point`](@ref). Only relevant if `k` is not specified. 
- `point_indices`: The indices for the points. Only relevant if `k` is not specified. 
- `try_points=()`: Extra points to try when sampling an initial point from [`select_initial_point`](@ref). Only relevant if `k` is not specified. 
- `k=select_initial_point(pts, q; m, point_indices, try_points)`: Where to start the algorithm.
- `check_existence::C=Val(has_multiple_segments(tri))`: Whether to check that the edge exists when using [`get_adjacent`](@ref), helping to correct for incorrect boundary indices in the presence of multiple boundary segments. See [`get_adjacent`](@ref).

# Output 
Returns the triangle `V` containing the query point `q`.
""" jump_and_march(::Triangulation, ::Any)

@doc """
    integer_type(tri::Triangulation{P,Ts,I}) where {P,Ts,I}

Given a triangulation `tri`, returns the type used for representing 
the integers.
""" integer_type(::Triangulation)

@doc """
    number_type(::Triangulation{P}) where {P}

Given a triangulation `tri`, returns the type used for representing 
individual coordinates.
""" number_type(::Triangulation)

@doc """
    compute_representative_points!(tri::Triangulation; use_convex_hull=!has_multiple_segments(tri) && num_boundary_edges(get_boundary_nodes(tri)) == 0))

Given a triangulation `tri`, computes representative points for each region using 
[`pole_of_inaccessibility`](@ref). If you only want to update the point for the convex hull,
set `use_convex_hull = true` - this will only update `RepresentativePointList[1]`.
        
See also [`RepresentativePointList`](@ref).

!!! warning
    While this function computes an appropriate visual center of the polygon represented
    by the curves, i.e. by joining points, the update functions like `update_centroid_after_addition` 
    and `update_centroid_after_deletion` treat the centroid as if it is a 
    mean of all the points. If you need to use the actual visual centers of the 
    polygons, you will need to reuse this function again (as is 
    done at the end of [`triangulate`](@ref)).
""" compute_representative_points!(::Triangulation)

@doc """
    clear_empty_features!(tri::Triangulation)

Removes empty keys from the [`Adjacent`](@ref) and [`Adjacent2Vertex`](@ref)
maps from the triangulation `tri`, and all points with empty neighbourhoods 
from the [`Graph`](@ref). 

See also [`clear_empty_keys!`](@ref) and [`clear_empty_points!`](@ref).
""" clear_empty_features!(::Triangulation)

@doc """
    find_edge(tri::Triangulation, T, ℓ)

Given a triangle `T`, a triangulation `tri`, and a point index 
`ℓ` which is assumed to correspond to some point on an edge of `T`, 
returns the edge `(u, v)` that the point is on.
""" find_edge(::Triangulation, ::Any, ::Any)

@doc """
    is_constrained(tri::Triangulation)

Returns `true` if `tri` has any constrained edges, and `false` otherwise.
""" is_constrained(::Triangulation)

@doc """
    all_boundary_indices(tri::Triangulation)

Returns an iterator over all boundary indices in the triangulation.
""" all_boundary_indices(::Triangulation)

@doc """
    get_surrounding_polygon(tri::Triangulation, u; skip_boundary_indices=false)

Given a triangulation `tri` and a vertex `u`, returns a vector `S` which gives a counter-clockwise 
sequence of the neighbours of `u`.

When `u` is an outer boundary index, the returned polygon is clockwise.

When `u` is a boundary vertex and you do not have ghost triangles, then this function may return an invalid polygon.

If you want to remove all boundary indices from the result at the end, set `skip_boundary_indices=true`.
""" get_surrounding_polygon(::Triangulation, ::Any)