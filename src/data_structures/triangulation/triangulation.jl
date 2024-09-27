@doc """
    Triangulation{P,T,BN,W,I,E,Es,BC,BCT,BEM,GVM,GVR,BPL,C,BE}

Struct representing a triangulation, as constructed by [`triangulate`](@ref).

!!! note "Field access"

    Accessing the fields themselves using e.g. `tri.field` is not recommended and is not intended 
    to be in the public API. You should be using the 
    accessor functions, e.g. instead of `tri.points` do `get_points(tri)`. Similarly, for the iterators,
    e.g. `tri.triangles`, `each_triangle(tri)` is recommended instead.

# Fields 
- `points::P`

The point set of the triangulation. Please note that this may not necessarily correspond to each point 
in the triangulation, e.g. some points may have been deleted - see [`each_solid_vertex`](@ref) for an iterator over 
each vertex in the triangulation. 
- `triangles::T`

The triangles in the triangulation. Each triangle is oriented counter-clockwise. If your triangulation has ghost triangles,
some of these triangles will contain ghost vertices (i.e., vertices with negative indices). Solid triangles can be iterated over using 
[`each_solid_triangle`](@ref).
- `boundary_nodes::BN`

The boundary nodes of the triangulation, if the triangulation is constrained; the assumed form of these boundary nodes is outlined 
in the docs. If your triangulation is unconstrained, then `boundary_nodes` will be empty and the boundary should instead be inspected
using the convex hull field, or alternatively you can see [`lock_convex_hull!`](@ref).
- `interior_segments::Es`

Constrained segments appearing in the triangulation. These will only be those segments appearing off of the boundary. If your triangulation is unconstrained, then `segments` will be empty.
- `all_segments::Es`

This is similar to `segments`, except this includes both the interior segments and the boundary segments. If your triangulation is unconstrained, then `all_segments` will be empty.
- `weights::W`

The weights of the triangulation. If you are not using a weighted triangulation, this will be given by `ZeroWeight()`. Otherwise, 
the weights must be such that `get_weight(weights, i)` is the weight for the `i`th vertex. The weights should have the same type as the 
coordinates in `points`.
- `adjacent::Adjacent{I,E}`

The [`Adjacent`](@ref) map of the triangulation. This maps edges `(u, v)` to vertices `w` such that `(u, v, w)` is a positively 
oriented triangle in `triangles` (up to rotation).
- `adjacent2vertex::Adjacent2Vertex{I,Es}`

The [`Adjacent2Vertex`](@ref) map of the triangulation. This maps vertices `w` to sets `S` such that `(u, v, w)` is a positively 
oriented triangle in `triangles` (up to rotation) for all `(u, v) ∈ S`.
- `graph::Graph{I}`

The [`Graph`](@ref) of the triangulation, represented as an undirected graph that defines all the neighbourhood information for the triangulation.
- `boundary_curves::BC`

Functions defining the boundary curves of the triangulation, incase you are triangulating a curve-bounded domain. By default, this will be an empty `Tuple`,
indicating that the boundary is as specified in `boundary_nodes` - a piecewise linear curve. If you are triangulating a curve-bounded domain, then these will 
be the parametric curves (see [`AbstractParametricCurve`](@ref)) you provided as a `Tuple`, where the `i`th element of the `Tuple` is associated with the ghost vertex `-i`, i.e. the `i`th section as indicated by 
`ghost_vertex_map`. If the `i`th boundary was left was a sequence of edges, then the function will be a `PiecewiseLinear()`.
- `boundary_edge_map::BEM`

This is a `Dict` from [`construct_boundary_edge_map`](@ref) that maps boundary edges `(u, v)` to their corresponding position in `boundary_nodes`.
- `ghost_vertex_map::GVM`

This is a `Dict` that maps ghost vertices to their corresponding section in `boundary_nodes`, constructed by [`construct_ghost_vertex_map`](@ref).
- `ghost_vertex_ranges::GVR`

This is a `Dict` that maps ghost vertices to a range of all other ghost vertices that appear on the curve corresponding to the given ghost vertex, 
constructed by [`construct_ghost_vertex_ranges`](@ref).
- `convex_hull::ConvexHull{P,I}`

The [`ConvexHull`](@ref) of the triangulation, which is the convex hull of the point set `points`.

- `representative_point_list::BPL`

The `Dict` of points giving [`RepresentativeCoordinates`](@ref) for each boundary curve, or for the 
convex hull if `boundary_nodes` is empty. These representative points are used for interpreting 
ghost vertices.
- `polygon_hierarchy::PolygonHierarchy{I}`

The [`PolygonHierarchy`](@ref) of the boundary, defining the hierarchy of the boundary curves, giving information about which curves are contained in which other curves.
- `boundary_enricher::BE`

The [`BoundaryEnricher`](@ref) used for triangulating a curve-bounded domain. If the domain is not curve-bounded, this is `nothing`.
- `cache::C`

A [`TriangulationCache`](@ref) used as a cache for [`add_segment!`](@ref) which requires a separate `Triangulation` structure for use.
This will not contain any segments or boundary nodes. Also stores segments useful for [`lock_convex_hull!`](@ref) and [`unlock_convex_hull!`](@ref).
"""
struct Triangulation{P,T,BN,W,I,E,Es,BC,BEM,GVM,GVR,BPL,C,BE}
    # Geometry
    points::P
    triangles::T
    boundary_nodes::BN
    interior_segments::Es
    all_segments::Es
    weights::W
    # Topology
    adjacent::Adjacent{I,E}
    adjacent2vertex::Adjacent2Vertex{I,Es}
    graph::Graph{I}
    # Boundary handling
    boundary_curves::BC
    boundary_edge_map::BEM
    ghost_vertex_map::GVM
    ghost_vertex_ranges::GVR
    # Other 
    convex_hull::ConvexHull{P,I}
    representative_point_list::BPL
    polygon_hierarchy::PolygonHierarchy{I}
    boundary_enricher::BE
    cache::C
end
function Base.show(io::IO, ::MIME"text/plain", tri::Triangulation)
    println(io, "Delaunay Triangulation.")
    println(io, "   Number of vertices: ", num_solid_vertices(tri))
    println(io, "   Number of triangles: ", num_solid_triangles(tri))
    println(io, "   Number of edges: ", num_solid_edges(tri))
    println(io, "   Has boundary nodes: ", has_boundary_nodes(tri))
    println(io, "   Has ghost triangles: ", has_ghost_triangles(tri))
    println(io, "   Curve-bounded: ", is_curve_bounded(tri))
    println(io, "   Weighted: ", is_weighted(tri))
    print(io, "   Constrained: ", is_constrained(tri))
end

function Base.:(==)(tri1::Triangulation, tri2::Triangulation)
    !has_ghost_triangles(tri1) && has_ghost_triangles(tri2) && return false
    !has_ghost_triangles(tri2) && has_ghost_triangles(tri1) && return false
    get_points(tri1) ≠ get_points(tri2) && return false
    !compare_triangle_collections(get_triangles(tri1), get_triangles(tri2)) && return false
    !compare_unoriented_edge_collections(get_interior_segments(tri1), get_interior_segments(tri2)) && return false
    !compare_unoriented_edge_collections(get_all_segments(tri1), get_all_segments(tri2)) && return false
    get_adjacent(tri1) ≠ get_adjacent(tri2) && return false
    get_adjacent2vertex(tri1) ≠ get_adjacent2vertex(tri2) && return false
    get_graph(tri1) ≠ get_graph(tri2) && return false
    get_boundary_edge_map(tri1) ≠ get_boundary_edge_map(tri2) && return false
    get_ghost_vertex_map(tri1) ≠ get_ghost_vertex_map(tri2) && return false
    get_ghost_vertex_ranges(tri1) ≠ get_ghost_vertex_ranges(tri2) && return false
    rep1 = get_representative_point_list(tri1)
    rep2 = get_representative_point_list(tri2)
    length(rep1) ≠ length(rep2) && return false
    for i in 1:length(rep1)
        p1 = get_representative_point_coordinates(tri1, i)
        p2 = get_representative_point_coordinates(tri2, i)
        !isapprox([getx(p1), gety(p1)], [getx(p2), gety(p2)], atol=1e-10) && return false
    end
    get_polygon_hierarchy(tri1) ≠ get_polygon_hierarchy(tri2) && return false
    get_boundary_nodes(tri1) ≠ get_boundary_nodes(tri2) && return false
    return true
end

"""
    get_points(tri::Triangulation) -> Points 

Return the points of the triangulation. Note that this may include points not yet in `tri`.
"""
get_points(tri::Triangulation) = tri.points

"""
    get_triangles(tri::Triangulation) -> Triangles

Return the triangles of the triangulation. These triangles are all given in counter-clockwise order, 
and may include ghost triangles.
"""
get_triangles(tri::Triangulation) = tri.triangles

"""
    get_boundary_nodes(tri::Triangulation) -> BoundaryNodes 

Return the boundary nodes of the triangulation. This is only for triangulations with a constrained boundary. If the triangulation 
has no constrained boundary, then the boundary is instead given by its convex hull and this function returns an empty vector. See 
[`get_convex_hull`](@ref).
"""
get_boundary_nodes(tri::Triangulation) = tri.boundary_nodes

"""
    get_interior_segments(tri::Triangulation) -> Edges 

Return the interior segments of the triangulation. These are segments that are forced to be in the triangulation - 
they are not the same as edges.
"""
get_interior_segments(tri::Triangulation) = tri.interior_segments

"""
    get_all_segments(tri::Triangulation) -> Edges

Return all segments of the triangulation. This includes interior segments and boundary segments. Segments are 
edges that are forced to be in the triangulation.
"""
get_all_segments(tri::Triangulation) = tri.all_segments

"""
    get_weights(tri::Triangulation) -> Weights

Return the weights of the triangulation. These are the weights of the vertices of the triangulation.
"""
get_weights(tri::Triangulation) = tri.weights

"""
    get_adjacent(tri::Triangulation) -> Adjacent

Returns the adjacency map of the triangulation. This is a map from each edge `(u, v)` to a vertex `w` such that `(u, v, w)` 
is a positively oriented triangle in `tri`. 

See also [`Adjacent`](@ref).
"""
get_adjacent(tri::Triangulation) = tri.adjacent

"""
    get_adjacent2vertex(tri::Triangulation) -> Adjacent2Vertex

Returns the [`Adjacent2Vertex`](@ref) map of the triangulation `tri`. This is a map from a vertex `w` to a set of 
all edges `(u, v)` such that `(u, v, w)` is a positively oriented triangle in `tri`.
"""
get_adjacent2vertex(tri::Triangulation) = tri.adjacent2vertex

"""
    get_graph(tri::Triangulation) -> Graph

Returns the [`Graph`](@ref) of the triangulation `tri`. This is an undirected graph.
"""
get_graph(tri::Triangulation) = tri.graph

"""
    get_boundary_curves(tri::Triangulation) -> NTuple{N, Function}

Returns the functions defining the boundaries of `tri`. If `!is_curve_bounded(tri)`, 
then this returns an empty `Tuple`. Otherwise, this returns a `Tuple` of functions, one for each section of the boundary, 
where the `i`th element of the `Tuple` corresponds to the `i`th section of the boundary, which corresponds to the ghost vertex `-i`.
For curves that are defined by boundary nodes rather than by a function, the function is [`PiecewiseLinear`](@ref). For the other functions, these 
are all defined by `t -> NTuple{2, Number}`, where `t ∈ [0, 1]` and the `NTuple{2, Number}` is the coordinate on the curve at that `t`.
"""
get_boundary_curves(tri::Triangulation) = tri.boundary_curves

"""
    get_boundary_edge_map(tri::Triangulation) -> Dict

Returns the boundary edge map of the triangulation `tri`. This is a `Dict` that maps 
a boundary edge `(u, v)` to its position in `get_boundary_nodes(tri)`. In particular, 
the returned value is a `Tuple` `(position, index)` so that `boundary_nodes = get_boundary_nodes(tri, position)` are the boundary nodes associated 
with the section that `(u, v)` resides on, and `u = get_boundary_nodes(boundary_nodes, index)` and 
`v = get_boundary_nodes(boundary_nodes, index + 1)`.

See also [`construct_boundary_edge_map`](@ref).
"""
get_boundary_edge_map(tri::Triangulation) = tri.boundary_edge_map

"""
    get_ghost_vertex_map(tri::Triangulation) -> Dict

Returns the ghost vertex map of the triangulation `tri`. This is a `Dict` that maps ghost vertices 
to their associated section in `boundary_nodes`. There are three cases; below, `I` is `integer_type(tri)`:

- `has_multiple_curves(tri)`

Returns `dict::Dict{I, NTuple{2, I}}`, mapping ghost vertices `i` to `Tuple`s `(m, n)`
so that `get_boundary_nodes(tri, m, n)` are the boundary nodes associated with `i`, 
i.e. the `n`th section of the `m`th curve is associated with the ghost vertex `i`.
- `has_multiple_sections(tri)`

Returns `dict::Dict{I, I}`, mapping ghost vertices `i` to `n` so that
`get_boundary_nodes(tri, n)` are the boundary nodes associated with `i`,
i.e. the `n`th section of the boundary is associated with the ghost vertex `i`.
- `otherwise`

Returns `dict::Dict{I, A}`, mapping the ghost vertex `i` to `get_boundary_nodes(tri)`, where `A = typeof(get_boundary_nodes(tri))`.

See also [`construct_ghost_vertex_map`](@ref).
"""
get_ghost_vertex_map(tri::Triangulation) = tri.ghost_vertex_map

"""
    get_ghost_vertex_ranges(tri::Triangulation) -> Dict

Returns the ghost vertex ranges map of the triangulation `tri`. This is a `Dict` that maps ghost vertices `i` 
to the range of all other ghost vertices associated with the curve that `i` is associated with. 

See also [`construct_ghost_vertex_ranges`](@ref).
"""
get_ghost_vertex_ranges(tri::Triangulation) = tri.ghost_vertex_ranges

"""
    get_convex_hull(tri::Triangulation) -> ConvexHull 

Returns the convex hull of the points in `tri`. This is given as a [`ConvexHull`](@ref) object, where the vertices 
are sorted counter-clockwise and defined so that the first and last vertices are equal.
"""
get_convex_hull(tri::Triangulation) = tri.convex_hull

"""
    get_representative_point_list(tri::Triangulation) -> Dict 

Returns the `Dict` of [`RepresentativeCoordinates`](@ref) of `tri`, mapping curve indices `i` to the representative point for that 
curve. These representative points are how we interpret ghost triangles relative to that curve.
"""
get_representative_point_list(tri::Triangulation) = tri.representative_point_list

"""
    get_polygon_hierarchy(tri::Triangulation) -> PolygonHierarchy

Returns the [`PolygonHierarchy`](@ref) of the boundary of `tri`. This defines the hierarchy of the boundary curves, giving information about which curves are contained in which other curves.
"""
get_polygon_hierarchy(tri::Triangulation) = tri.polygon_hierarchy

"""
    get_exterior_curve_indices(tri::Triangulation) -> KeySet{Integer}

Returns the set of all curve indices that correspond to exterior curves of `tri`.
"""
get_exterior_curve_indices(tri::Triangulation) = get_exterior_curve_indices(get_polygon_hierarchy(tri))

"""
    get_boundary_enricher(tri::Triangulation) -> BoundaryEnricher

Returns the [`BoundaryEnricher`](@ref) of `tri`. If the domain is not curve-bounded, this is `nothing`.
"""
get_boundary_enricher(tri::Triangulation) = tri.boundary_enricher

"""
    get_cache(tri::Triangulation) -> TriangulationCache

Returns the cache of `tri`. This is a [`TriangulationCache`](@ref) used as a cache for [`add_segment!`](@ref) which requires a separate `Triangulation` structure for use.
"""
get_cache(tri::Triangulation) = tri.cache

"""
    get_incircle_cache(tri::Triangulation) -> Tuple 

Returns the incircle cache stored in `tri`.
"""
get_incircle_cache(tri::Triangulation) = get_incircle_cache(get_cache(tri))

"""
    get_orient3_cache(tri::Triangulation) -> Tuple

Returns the orient3 cache stored in `tri`.
"""
get_orient3_cache(tri::Triangulation) = get_orient3_cache(get_cache(tri))

"""
    get_insphere_cache(tri::Triangulation) -> Tuple

Returns the insphere cache stored in `tri`.
"""
get_insphere_cache(tri::Triangulation) = get_insphere_cache(get_cache(tri))

"""
    number_type(tri::Triangulation) -> DataType 

Returns the type used for representing individual coordinates in `tri`.
"""
number_type(::Triangulation{P}) where {P} = number_type(P)

"""
    integer_type(tri::Triangulation) -> DataType

Returns the type used for representing vertices in `tri`.
"""
integer_type(::Triangulation{P,T,BN,W,I}) where {P,T,BN,W,I} = I

"""
    edges_type(tri::Triangulation) -> DataType

Returns the type used for representing collections of edges in `tri`.
"""
edges_type(::Triangulation{P,T,BN,W,I,E,Es}) where {P,T,BN,W,I,E,Es} = Es

"""
    edge_type(tri::Triangulation) -> DataType

Returns the type used for representing individual edges in `tri`.
"""
edge_type(::Triangulation{P,T,BN,W,I,E}) where {P,T,BN,W,I,E} = E

"""
    triangle_type(tri::Triangulation) -> DataType

Returns the type used for representing collections of triangles in `tri`.
"""
triangles_type(::Triangulation{P,T}) where {P,T} = T

"""
    triangle_type(tri::Triangulation) -> DataType

Returns the type used for representing individual triangles in `tri`.
"""
triangle_type(tri::Triangulation) = triangle_type(triangles_type(tri))

@inline function __Triangulation(
    points::P, boundary_nodes, IntegerType::Type{I}, EdgeType::Type{E},
    TriangleType::Type{V}, EdgesType::Type{Es}, TrianglesType::Type{Ts},
) where {P,I,E,V,Es,Ts}
    T = TrianglesType()
    adj = Adjacent{IntegerType,EdgeType}()
    adj2v = Adjacent2Vertex{IntegerType,EdgesType}()
    graph = Graph{IntegerType}()
    all_segments = EdgesType()
    boundary_edge_map = construct_boundary_edge_map(boundary_nodes, IntegerType, EdgeType)
    ghost_vertex_map = construct_ghost_vertex_map(boundary_nodes, IntegerType)
    ghost_vertex_ranges = construct_ghost_vertex_ranges(boundary_nodes, IntegerType)
    ch = ConvexHull(points, IntegerType[])
    n = num_points(points)
    sizehint!(T, 2n - 5) # maximum number of triangles
    sizehint!(adj, 3n - 6) # maximum number of edges 
    sizehint!(adj2v, n)
    sizehint!(graph, 3n - 6, n, n)
    sizehint!(ch, n)
    polygon_hierarchy = construct_polygon_hierarchy(points; IntegerType)
    return T, adj, adj2v, graph, all_segments, boundary_edge_map, ghost_vertex_map, ghost_vertex_ranges, ch, polygon_hierarchy
end

@inline function _build_cache(
    points::P, IntegerType::Type{I}, EdgeType::Type{E},
    TriangleType::Type{V}, EdgesType::Type{Es}, TrianglesType::Type{T}, weights::W, build_cache::Val{B},
) where {P,I,E,V,Es,T,W,B}
    if B
        NT = number_type(points)
        BN = Vector{IntegerType}
        BC = Tuple{}
        BEM = Dict{EdgeType,Tuple{Vector{IntegerType},IntegerType}}
        GVM = Dict{IntegerType,Vector{IntegerType}}
        GVR = Dict{IntegerType,UnitRange{IntegerType}}
        BPL = Dict{IntegerType,RepresentativeCoordinates{IntegerType,NT}}
        incircle_cache = AP.incircleadapt_cache(NT)
        orient3_cache = AP.orient3adapt_cache(NT)
        insphere_cache = AP.insphereexact_cache(NT)
        IC = typeof(incircle_cache)
        OC = typeof(orient3_cache)
        IS = typeof(insphere_cache)
        C = EmptyTriangulationCache
        cache = TriangulationCache(
            Triangulation(points; weights, IntegerType, EdgeType, TriangleType, EdgesType, TrianglesType, build_cache=Val(false)),
            Triangulation(points; weights, IntegerType, EdgeType, TriangleType, EdgesType, TrianglesType, build_cache=Val(false)),
            I[], Es(), I[], T(), incircle_cache, orient3_cache, insphere_cache
        )::TriangulationCache{
            Triangulation{P,T,BN,W,I,E,Es,BC,BEM,GVM,GVR,BPL,C,Nothing},Vector{I},Es,Vector{I},T,IC,OC,IS
        }
    else
        cache = TriangulationCache()::EmptyTriangulationCache
    end
    return cache
end

"""
    Triangulation(points; kwargs...) -> Triangulation 

Initialises an empty `Triangulation` for triangulating `points`. The keyword arguments 
`kwargs...` match those of [`triangulate`](@ref).
"""
@inline function Triangulation(
    points::P;
    IntegerType::Type{I}=Int,
    EdgeType::Type{E}=NTuple{2,IntegerType},
    TriangleType::Type{T}=NTuple{3,IntegerType},
    EdgesType::Type{Es}=Set{EdgeType},
    TrianglesType::Type{Vs}=Set{TriangleType},
    boundary_nodes=IntegerType[],
    segments=EdgesType(),
    weights=ZeroWeight(),
    representative_point_list=Dict{IntegerType,RepresentativeCoordinates{IntegerType,number_type(points)}}(),
    boundary_curves=(),
    build_cache=true,
    boundary_enricher=nothing,
) where {P,I,E,T,Es,Vs}
    return _Triangulation(
        points,
        IntegerType, EdgeType, TriangleType, EdgesType, TrianglesType,
        boundary_nodes, segments, weights, representative_point_list,
        boundary_curves, _to_val(build_cache), boundary_enricher,
    )
end
@inline function _Triangulation(
    points, ::Type{I}, ::Type{E}, ::Type{V}, ::Type{Es}, ::Type{Ts},
    boundary_nodes, segments, weights, representative_point_list,
    boundary_curves, build_cache::Val{B}, boundary_enricher,
) where {I,E,V,Es,Ts,B}
    T, adj, adj2v, graph, all_segments, boundary_edge_map, ghost_vertex_map, ghost_vertex_ranges, ch, polygon_hierarchy = __Triangulation(points, boundary_nodes, I, E, V, Es, Ts)
    cache = _build_cache(points, I, E, V, Es, Ts, weights, build_cache)
    tri = Triangulation(
        points, T, boundary_nodes, segments, all_segments, weights,
        adj, adj2v, graph, boundary_curves, boundary_edge_map,
        ghost_vertex_map, ghost_vertex_ranges, ch, representative_point_list, polygon_hierarchy, boundary_enricher, cache,
    )
    return tri
end

"""
    Triangulation(points, triangles, boundary_nodes; kwargs...) -> Triangulation

Returns the `Triangulation` corresponding to the triangulation of `points` with `triangles` and `boundary_nodes`.

# Arguments 
- `points`: The points that the triangulation is of. 
- `triangles`: The triangles of the triangulation. These should be given in counter-clockwise order, with vertices corresponding to `points`. These should not include any ghost triangles.
- `boundary_nodes`: The boundary nodes of the triangulation. These should match the specification given in the documentation or in [`check_args`](@ref).

# Keyword Arguments 
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.
- `IntegerType=Int`: The integer type to use for the triangulation. This is used for representing vertices.
- `EdgeType=isnothing(segments) ? NTuple{2,IntegerType} : (edge_type ∘ typeof)(segments)`: The edge type to use for the triangulation. 
- `TriangleType=NTuple{3,IntegerType}`: The triangle type to use for the triangulation.
- `EdgesType=isnothing(segments) ? Set{EdgeType} : typeof(segments)`: The type to use for storing the edges of the triangulation.
- `TrianglesType=Set{TriangleType}`: The type to use for storing the triangles of the triangulation.
- `weights=ZeroWeight()`: The weights associated with the triangulation. 
- `delete_ghosts=false`: Whether to delete the ghost triangles after the triangulation is computed. This is done using [`delete_ghost_triangles!`](@ref).

# Output 
- `tri`: The [`Triangulation`](@ref).
"""
@inline function Triangulation(
    points::P, triangles::T, boundary_nodes::BN;
    IntegerType::Type{I}=Int,
    EdgeType::Type{E}=NTuple{2,IntegerType},
    TriangleType::Type{V}=NTuple{3,IntegerType},
    EdgesType::Type{Es}=Set{EdgeType},
    TrianglesType::Type{Ts}=Set{TriangleType},
    weights=ZeroWeight(),
    delete_ghosts=false,
    predicates::AbstractPredicateKernel=AdaptiveKernel(),
) where {P,T,BN,I,E,V,Es,Ts}
    _bn = copy(boundary_nodes)
    tri = Triangulation(points; boundary_nodes=_bn, weights, IntegerType, EdgeType, TriangleType, EdgesType, TrianglesType)
    return build_triangulation_from_data!(tri, triangles, _bn, delete_ghosts, predicates)
end

"""
    Triangulation(points, triangles; kwargs...) -> Triangulation

Returns the `Triangulation` corresponding to the triangulation of `points` with `triangles`.

# Arguments 
- `points`: The points that the triangulation is of. 
- `triangles`: The triangles of the triangulation. These should be given in counter-clockwise order, with vertices corresponding to `points`. These should not include any ghost triangles.

# Keyword Arguments 
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.
- `IntegerType=Int`: The integer type to use for the triangulation. This is used for representing vertices.
- `EdgeType=isnothing(segments) ? NTuple{2,IntegerType} : (edge_type ∘ typeof)(segments)`: The edge type to use for the triangulation. 
- `TriangleType=NTuple{3,IntegerType}`: The triangle type to use for the triangulation.
- `EdgesType=isnothing(segments) ? Set{EdgeType} : typeof(segments)`: The type to use for storing the edges of the triangulation.
- `TrianglesType=Set{TriangleType}`: The type to use for storing the triangles of the triangulation.
- `weights=ZeroWeight()`: The weights associated with the triangulation. 
- `delete_ghosts=false`: Whether to delete the ghost triangles after the triangulation is computed. This is done using [`delete_ghost_triangles!`](@ref).

# Output 
- `tri`: The [`Triangulation`](@ref).
"""
@inline function Triangulation(
    points::P, triangles::T;
    IntegerType::Type{I}=Int,
    EdgeType::Type{E}=NTuple{2,IntegerType},
    TriangleType::Type{V}=NTuple{3,IntegerType},
    EdgesType::Type{Es}=Set{EdgeType},
    TrianglesType::Type{Ts}=Set{TriangleType},
    weights=ZeroWeight(),
    delete_ghosts=false,
    predicates::AbstractPredicateKernel=AdaptiveKernel(),
) where {P,T,I,E,V,Es,Ts}
    _bn = [nothing]
    tri = Triangulation(points; boundary_nodes=_bn, weights, IntegerType, EdgeType, TriangleType, EdgesType, TrianglesType)
    return build_triangulation_from_data!(tri, triangles, _bn, delete_ghosts, predicates)
end

"""
    build_triangulation_from_data!(tri::Triangulation, triangles, boundary_nodes, delete_ghosts, predicates::AbstractPredicateKernel=AdaptiveKernel())

Given an empty `triangulation`, `tri`, adds all the `triangles` and `boundary_nodes` into it. Use 
`delete_ghosts=true` if you want to have all ghost triangles deleted afterwards.

The `kernel` argument determines how predicates are computed, and should be 
one of [`ExactKernel`](@ref), [`FastKernel`](@ref), and [`AdaptiveKernel`](@ref) (the default).
See the documentation for more information about these choices.
"""
@inline function build_triangulation_from_data!(tri::Triangulation, triangles, boundary_nodes, delete_ghosts, predicates::AbstractPredicateKernel=AdaptiveKernel())
    Es = edges_type(tri)
    points = get_points(tri)
    polygon_hierarchy = get_polygon_hierarchy(tri)
    construct_polygon_hierarchy!(polygon_hierarchy, points, boundary_nodes)
    adj = get_adjacent(tri)
    adj2v = get_adjacent2vertex(tri)
    graph = get_graph(tri)
    tris = get_triangles(tri)
    for τ in each_triangle(triangles)
        add_triangle!(adj, τ)
        add_triangle!(adj2v, τ)
        add_triangle!(graph, τ)
        add_triangle!(tris, τ)
    end
    add_boundary_information!(tri)
    add_ghost_triangles!(tri)
    convex_hull!(tri; predicates, reconstruct=true)
    segments = get_all_segments(tri)
    ghost_vertex_map = get_ghost_vertex_map(tri)
    all_segments = merge_segments(ghost_vertex_map, boundary_nodes, Es())
    for edge in each_edge(all_segments)
        add_edge!(segments, edge)
    end
    postprocess_triangulate!(tri; delete_ghosts)
    return tri
end
