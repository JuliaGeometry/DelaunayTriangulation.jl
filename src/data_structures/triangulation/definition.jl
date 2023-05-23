"""
    Triangulation{P,Ts,I,E,Es,BN,BNM,B,BIR,BPL}

Struct representing a Delaunay triangulation, as constructed via e.g. [`triangulate`](@ref) or [`generate_mesh`](@ref).

# Fields 
- `points::P`

The points used to construct the triangulation. Do note that this is 
not necessarly the same as the set of points appearing in the triangulation, as could happen 
e.g. if some points were deleted from the triangulation. If that is what you need, 
then see [`each_vertex`](@ref).
- `triangles::Ts`

The triangles in the triangulation, all given with positive orientation. This field can include 
ghost triangles. If you don't want ghost triangles, see [`each_solid_triangle`](@ref).
- `adjacent::Adjacent{I, E}`

The [`Adjacent`](@ref) map.
- `adjacent2vertex::Adjacent2Vertex{I,Es,E}`

The [`Adjacent2Vertex`](@ref) map.
- `graph::Graph{I}`

The [`Graph`](@ref). 
- `boundary_nodes::BN`

The boundary nodes in the triangulation that you provide, matching the specification given in the 
documentation.
- `boundary_edge_map::BNM`

This is a `Dict` that maps all of the boundary edges to their position 
in `boundary_nodes`. See also [`construct_boundary_edge_map`](@ref).
- `boundary_index_ranges::BIR`

This is an `OrderedDict` that maps a boundary index to a range of all other boundary indices 
that the corresponding boundary curve could correspond to. For example, for a curve with 
four segments, there are four possible boundary indices that a segment could correspond to. 
See also [`construct_boundary_index_ranges`](@ref).
- `constrained_edges::Es`

This is the set of extra edges added into the triangulation that you have provided.
This will not include any of the other constrained edges from `boundary_nodes`.
- `all_constrained_edges::Es`

This is the set of all constrained edges appearing in the triangulation, essentially 
given as the union of `constrained_edges` and `boundary_nodes`. This is different from 
the `constrained_edges` field as it contains edges from both of these fields, e.g. 
there might be an edge in `constrained_edges` that is not yet in the triangulation, 
but it will definitely not appear in `all_constrained_edges`.
- `convex_hull::ConvexHull{P,Vector{I}}`

The [`ConvexHull`](@ref) of `points`. See also [`convex_hull`](@ref).
- `representative_point_list::RPL`

The `Dict` of points giving representative points for each boundary curve, or for the 
convex hull if `boundary_nodes` is empty. These representative points are used for interpreting 
ghost vertices.

# Constructors 

There are several ways to construct this struct directly, although in most cases 
you should be using [`triangulate`](@ref) or [`generate_mesh`](@ref).

## Default Constructor 

The default constructor is available, i.e. 

    Triangulation(
        points,
        triangles,
        adjacent,
        adjacent2vertex,
        graph,
        boundary_nodes,
        boundary_edge_map,
        boundary_map,
        boundary_index_ranges,
        constrained_edges,
        all_constrained_edges,
        convex_hull,
        representative_point_list
    )

## Empty Triangulation 

An empty triangulation can be initalised with the following method, 

    Triangulation(points;
        IntegerType=Int,
        EdgeType=NTuple{2,IntegerType},
        TriangleType=NTuple{3,IntegerType},
        EdgesType=Set{EdgeType},
        TrianglesType=Set{TriangleType},
        boundary_nodes=IntegerType[],
        constrained_edges=initialise_edges(EdgesType),
        representative_point_list=get_empty_representative_points(IntegerType, number_type(points))
    )

## Triangulation From an Existing Mesh 

A method is available from constructing a mesh from an existing set of points, 
triangles, and boundary nodes, mainly existing for the purpose of [`generate_mesh`](@ref):

    Triangulation(points, triangles, boundary_nodes;
        IntegerType=Int,
        EdgeType=NTuple{2,IntegerType},
        TriangleType=NTuple{3,IntegerType},
        EdgesType=Set{EdgeType},
        TrianglesType=Set{TriangleType},
        add_ghost_triangles=false
    )

with `add_ghost_triangles` calling [`add_ghost_triangles!`](@ref) at the end of the constructor.

## Wrapping a Triangulation with constrained_edge_points

This method is used in `triangulate` for wrapping a triangulation with a set of constraints:

    remake_triangulation_with_constraints(triangulation, edges, boundary_nodes)

which returns `(bn_map, bn_range, tri)`, where `bn_map` is the `boundary_map` that isn't yet added,
`bn_range` is the `boundary_index_range` that isn't yet added either, and `tri` is the wrapped 
triangulation that includes the constrained `edges` and the `boundary_nodes`. Note that 
either `edges` or `boundary_nodes` can be `nothing`.

This can be used together with the method

    replace_boundary_dict_information(triangulation, bn_map, bn_range)

which returns a new triangulation with the `boundary_map` and `boundary_index_range` replaced,
completing the wrapper.
"""
struct Triangulation{P,Ts,I,E,Es,BN,BNM,B,BIR,BPL}
    points::P
    triangles::Ts
    adjacent::Adjacent{I,E}
    adjacent2vertex::Adjacent2Vertex{I,Es,E}
    graph::Graph{I}
    boundary_nodes::BN
    boundary_edge_map::BNM
    boundary_map::B
    boundary_index_ranges::BIR
    constrained_edges::Es
    all_constrained_edges::Es
    convex_hull::ConvexHull{P,Vector{I}}
    representative_point_list::BPL
end
for n in fieldnames(Triangulation)
    name = String(n)
    @eval begin
        @doc """
        get_$($(name))(tri::Triangulation)

    Returns the $($name) field from the triangulation `tri`.
    """ ($(Symbol("get_$n")))(tri::Triangulation) = tri.$n
    end
end

function Base.show(io::IO, ::MIME"text/plain", tri::Triangulation)
    println(io, "Delaunay Triangulation.")
    println(io, "    Constrained: $(is_constrained(tri))")
    println(io, "    Has ghost triangles: $(has_ghost_triangles(tri))")
    println(io, "    Number of points: $(num_points(tri))")
    println(io, "    Number of triangles: $(num_triangles(tri))")
    print(io, "    Number of edges: $(num_edges(tri))")
end

@inline integer_type(::Triangulation{P,Ts,I}) where {P,Ts,I} = I
@inline number_type(::Triangulation{P}) where {P} = number_type(P)
