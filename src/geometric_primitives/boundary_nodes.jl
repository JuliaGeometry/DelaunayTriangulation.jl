const AV = AbstractVector
@doc """
    has_multiple_curves(boundary_nodes) -> Bool 

Check if `boundary_nodes` has multiple curves.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> DelaunayTriangulation.has_multiple_curves([1, 2, 3, 1])
false

julia> DelaunayTriangulation.has_multiple_curves([[1, 2, 3], [3, 4, 1]])
false

julia> DelaunayTriangulation.has_multiple_curves([[[1, 2, 3], [3, 4, 1]], [[5, 6, 7, 8, 5]]])
true
```
"""
has_multiple_curves
has_multiple_curves(::AAA) where {F <: Number, A <: AV{F}, AA <: AV{A}, AAA <: AV{AA}} = true
has_multiple_curves(::AA) where {F <: Number, A <: AV{F}, AA <: AV{A}} = false
has_multiple_curves(::A) where {F <: Number, A <: AV{F}} = false
function has_multiple_curves(boundary_nodes::AV)
    #=
    Need this since, e.g.
    [
        [
            [1, 2, 3, 4, 5, 6, 7, 1]
        ],
        [
            [CircularArc((0.5, 0.0), (0.5, 0.0), (0.0, 0.0), positive=false)]
        ],
    ]
    has type Vector{Vector}, but it has multiple curves.
    =#
    if !isempty(boundary_nodes) && all(x -> typeof(x) <: AbstractVector && eltype(x) <: AbstractVector, boundary_nodes)
        return true
    else
        return false
    end
end
has_multiple_curves(::NTuple{N, <:Integer}) where {N} = false

@doc """
    has_multiple_sections(boundary_nodes) -> Bool

Check if `boundary_nodes` has multiple sections.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> DelaunayTriangulation.has_multiple_sections([1, 2, 3, 1])
false

julia> DelaunayTriangulation.has_multiple_sections([[1, 2, 3], [3, 4, 1]])
true

julia> DelaunayTriangulation.has_multiple_sections([[[1, 2, 3], [3, 4, 1]], [[5, 6, 7, 8, 5]]])
true
```    
"""
has_multiple_sections
has_multiple_sections(::AAA) where {F <: Number, A <: AV{F}, AA <: AV{A}, AAA <: AV{AA}} = true
has_multiple_sections(::AA) where {F <: Number, A <: AV{F}, AA <: AV{A}} = true
has_multiple_sections(::A) where {F <: Number, A <: AV{F}} = false
function has_multiple_sections(boundary_nodes::AV)
    #=
    Need this since, e.g.
        2-element Vector{Any}:
        [1, 2, 3]
        EllipticalArc[EllipticalArc((0.0, 0.0), 2.0, 0.5, (0.0, 1.0), 0.0, 3.141592653589793, (2.0, 0.0), (-2.0, 0.0))]
    =#
    if !isempty(boundary_nodes) && all(x -> typeof(x) <: AbstractVector, boundary_nodes)
        return true
    else
        return false
    end
end
has_multiple_sections(::NTuple{N, <:Integer}) where {N} = false

"""
    num_curves(boundary_nodes) -> Integer

Get the number of curves in `boundary_nodes`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> DelaunayTriangulation.num_curves([1, 2, 3, 1])
1

julia> DelaunayTriangulation.num_curves([[1, 2, 3], [3, 4, 1]])
1

julia> DelaunayTriangulation.num_curves([[[1, 2, 3], [3, 4, 1]], [[5, 6, 7, 8, 5]]])
2
``` 
"""
function num_curves(boundary_nodes)
    if has_multiple_curves(boundary_nodes)
        return length(boundary_nodes)
    else
        return 1
    end
end

"""
    num_sections(boundary_nodes) -> Integer

Assuming `boundary_nodes` has only one curve, get the number of sections in `boundary_nodes`.

# Examples
```jldoctest
julia> using DelaunayTriangulation

julia> DelaunayTriangulation.num_sections([1, 2, 3, 4, 5, 1])
1

julia> DelaunayTriangulation.num_sections([[1, 2, 3, 4], [4, 5, 1]])
2

julia> DelaunayTriangulation.num_sections([[1, 2, 3], [3, 4, 5, 6, 7, 8], [8, 9], [9, 1]])
4
```
"""
function num_sections(boundary_nodes)
    if has_multiple_sections(boundary_nodes)
        return length(boundary_nodes)
    else
        return 1
    end
end

"""
    num_boundary_edges(boundary_nodes) -> Integer

Get the number of boundary edges in `boundary_nodes`, assuming that `boundary_nodes` defines a 
boundary with only one curve and a single section.
"""
num_boundary_edges(boundary_nodes::AV) = max(0, length(boundary_nodes) - 1)
num_boundary_edges(boundary_nodes::NTuple{N, <:Integer}) where {N} = max(0, N - 1)

@doc """
    get_boundary_nodes(boundary_nodes, mnℓ...) 

Given a collection of `boundary_nodes`, returns the specified component of the collection.
There are several forms for the methods:

1. `get_boundary_nodes(boundary_nodes, m)`: If `boundary_nodes` has multiple curves, this returns the `m`th curve. If `boundary_nodes` has multiple sections, this returns the `m`th section. Otherwise, this returns the `m`th boundary node.
2. `get_boundary_nodes(boundary_nodes, m, n)`: If `boundary_nodes` has multiple curves, this returns the `n`th section of the `m`th curve. Otherwise, if `boundary_nodes` has multiple sections, this returns the `n`th boundary node of the `m`th section.
3. `get_boundary_nodes(boundary_nodes, (m, n))`: This is equivalent to `get_boundary_nodes(boundary_nodes, m, n)`.
4. `get_boundary_nodes(boundary_nodes::A, ::A)`: This just returns `boundary_nodes`.  

# Examples
```jldoctest
julia> using DelaunayTriangulation

julia> get_boundary_nodes([[[1, 2, 3, 4], [4, 5, 1]], [[6, 7, 8, 9], [9, 10, 6]]], 2)
2-element Vector{Vector{Int64}}:
 [6, 7, 8, 9]
 [9, 10, 6]

julia> get_boundary_nodes([[1, 2, 3, 4], [4, 5, 1]], 1)
4-element Vector{Int64}:
 1
 2
 3
 4

julia> get_boundary_nodes([1, 2, 3, 4, 5, 6, 1], 4)
4

julia> get_boundary_nodes([[[1, 2, 3, 4], [4, 5, 1]], [[6, 7, 8, 9], [9, 10, 6]]], 1, 2)
3-element Vector{Int64}:
 4
 5
 1

julia> get_boundary_nodes([[1, 2, 3, 4], [4, 5, 6, 1]], 2, 3)
6

julia> get_boundary_nodes([1, 2, 3, 4, 5, 1], [1, 2, 3, 4, 5, 1])
6-element Vector{Int64}:
 1
 2
 3
 4
 5
 1
```
"""
get_boundary_nodes
@inline get_boundary_nodes(boundary_nodes, m::Integer) = boundary_nodes[m]
@inline get_boundary_nodes(boundary_nodes, m::Integer, n::Integer) = get_boundary_nodes(get_boundary_nodes(boundary_nodes, m), n)
@inline get_boundary_nodes(boundary_nodes, (m, n)::NTuple{2, Integer}) = get_boundary_nodes(boundary_nodes, m, n) # for indexing from a boundary map 
@inline get_boundary_nodes(boundary_nodes::A, ::A) where {A} = boundary_nodes # for indexing from a boundary map
@inline get_boundary_nodes(boundary_nodes::A, ::A) where {A <: Tuple{Integer, Integer}} = boundary_nodes # ambiguity
@inline get_boundary_nodes(boundary_nodes::A, ::A) where {A <: Integer} = boundary_nodes # ambiguity

"""
    each_boundary_node(boundary_nodes) -> Iterator 

Returns an iterator that goves over each node in `boundary_nodes`.
It is assumed that `boundary_nodes` represents a single section or a single
contiguous boundary; if you do want to loop over every boundary nodes for a boundary 
with multiple sections, you should to see the result from [`construct_ghost_vertex_map`](@ref).

# Examples 
```jldoctest 
julia> using DelaunayTriangulation

julia> DelaunayTriangulation.each_boundary_node([7, 8, 19, 2, 17])
5-element Vector{Int64}:
  7
  8
 19
  2
 17

julia> DelaunayTriangulation.each_boundary_node([7, 8, 19, 2, 17, 7])
6-element Vector{Int64}:
  7
  8
 19
  2
 17
  7
```
"""
each_boundary_node(boundary_nodes) = boundary_nodes

"""
    construct_ghost_vertex_map(boundary_nodes::A, IntegerType::Type{I}=number_type(boundary_nodes)) where {A,I} -> Dict

Given a set of `boundary_nodes`, returns a `Dict` that maps ghost vertices
to their associated section in `boundary_nodes`. There are three cases:

- `has_multiple_curves(boundary_nodes)`

Returns `dict::Dict{I, NTuple{2, I}}`, mapping ghost vertices `i` to `Tuple`s `(m, n)`
so that `get_boundary_nodes(boundary_nodes, m, n)` are the boundary nodes associated with `i`, 
i.e. the `n`th section of the `m`th curve is associated with the ghost vertex `i`.
- `has_multiple_sections(boundary_nodes)`

Returns `dict::Dict{I, I}`, mapping ghost vertices `i` to `n` so that
`get_boundary_nodes(boundary_nodes, n)` are the boundary nodes associated with `i`,
i.e. the `n`th section of the boundary is associated with the ghost vertex `i`.
- `otherwise`

Returns `dict::Dict{I, A}`, mapping the ghost vertex `i` to `boundary_nodes`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> gv_map = DelaunayTriangulation.construct_ghost_vertex_map([1, 2, 3, 4, 5, 1])
Dict{Int64, Vector{Int64}} with 1 entry:
  -1 => [1, 2, 3, 4, 5, 1]

julia> gv_map = DelaunayTriangulation.construct_ghost_vertex_map([[17, 29, 23, 5, 2, 1], [1, 50, 51, 52], [52, 1]])
Dict{Int64, Int64} with 3 entries:
  -1 => 1
  -3 => 3
  -2 => 2

julia> gv_map = DelaunayTriangulation.construct_ghost_vertex_map([[[1, 5, 17, 18, 1]], [[23, 29, 31, 33], [33, 107, 101], [101, 99, 85, 23]]])
Dict{Int64, Tuple{Int64, Int64}} with 4 entries:
  -1 => (1, 1)
  -3 => (2, 2)
  -2 => (2, 1)
  -4 => (2, 3)
```

# Extended help
This map can be useful for iterating over all boundary nodes. For example, you can iterate 
over all sections of a boundary using:

```julia
gv_map = construct_ghost_vertex_map(boundary_nodes)
for (ghost_vertex, section) in gv_map 
    nodes = get_boundary_nodes(boundary_nodes, section)
    # do something with nodes
end
```

This works for any form of `boundary_nodes`.
"""
@inline function construct_ghost_vertex_map(boundary_nodes, IntegerType::Type{I} = number_type(boundary_nodes)) where {I}
    if has_multiple_curves(boundary_nodes)
        return _construct_ghost_vertex_map_multiple_curves(boundary_nodes, IntegerType)
    elseif has_multiple_sections(boundary_nodes)
        return _construct_ghost_vertex_map_multiple_sections(boundary_nodes, IntegerType)
    else
        return _construct_ghost_vertex_map_contiguous(boundary_nodes, IntegerType)
    end
end
function _construct_ghost_vertex_map_multiple_curves(boundary_nodes, IntegerType::Type{I}) where {I}
    dict = Dict{I, NTuple{2, I}}()
    nc = num_curves(boundary_nodes)
    current_idx = I(𝒢)
    for m in 1:nc
        bn_m = get_boundary_nodes(boundary_nodes, m)
        ns = num_sections(bn_m)
        for n in 1:ns
            dict[current_idx] = (m, n)
            current_idx -= 1
        end
    end
    return dict
end
function _construct_ghost_vertex_map_multiple_sections(boundary_nodes, IntegerType::Type{I}) where {I}
    dict = Dict{I, I}()
    ns = num_sections(boundary_nodes)
    current_idx = I(𝒢)
    for n in 1:ns
        dict[current_idx] = n
        current_idx -= 1
    end
    return dict
end
function _construct_ghost_vertex_map_contiguous(boundary_nodes, IntegerType::Type{I}) where {I}
    dict = Dict(I(𝒢) => boundary_nodes)
    return dict
end

"""
    construct_boundary_edge_map(boundary_nodes::A, IntegerType::Type{I}=number_type(boundary_nodes), EdgeType::Type{E}=NTuple{2,IntegerType}) where {A,I,E} -> Dict

Constructs a map that takes boundary edges in `boundary_nodes` to a `Tuple` giving the edge's position 
in `boundary_nodes`. In particular, if `dict = construct_boundary_edge_map(boundary_nodes)`, then
`dict[e] = (pos, ℓ)` so that `bn = get_boundary_nodes(boundary_nodes, pos)` gives the boundary nodes 
associated with the section that `e` lives on, and `get_boundary_nodes(bn, ℓ)` is the first vertex of `e`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> DelaunayTriangulation.construct_boundary_edge_map([17, 18, 15, 4, 3, 17])
Dict{Tuple{Int64, Int64}, Tuple{Vector{Int64}, Int64}} with 5 entries:
  (18, 15) => ([17, 18, 15, 4, 3, 17], 2)
  (3, 17)  => ([17, 18, 15, 4, 3, 17], 5)
  (17, 18) => ([17, 18, 15, 4, 3, 17], 1)
  (4, 3)   => ([17, 18, 15, 4, 3, 17], 4)
  (15, 4)  => ([17, 18, 15, 4, 3, 17], 3)

julia> DelaunayTriangulation.construct_boundary_edge_map([[5, 17, 3, 9], [9, 18, 13, 1], [1, 93, 57, 5]])
Dict{Tuple{Int64, Int64}, Tuple{Int64, Int64}} with 9 entries:
  (18, 13) => (2, 2)
  (17, 3)  => (1, 2)
  (9, 18)  => (2, 1)
  (13, 1)  => (2, 3)
  (3, 9)   => (1, 3)
  (93, 57) => (3, 2)
  (5, 17)  => (1, 1)
  (57, 5)  => (3, 3)
  (1, 93)  => (3, 1)

julia> DelaunayTriangulation.construct_boundary_edge_map([[[2, 5, 10], [10, 11, 2]], [[27, 28, 29, 30], [30, 31, 85, 91], [91, 92, 27]]])
Dict{Tuple{Int64, Int64}, Tuple{Tuple{Int64, Int64}, Int64}} with 12 entries:
  (92, 27) => ((2, 3), 2)
  (2, 5)   => ((1, 1), 1)
  (11, 2)  => ((1, 2), 2)
  (10, 11) => ((1, 2), 1)
  (30, 31) => ((2, 2), 1)
  (91, 92) => ((2, 3), 1)
  (29, 30) => ((2, 1), 3)
  (31, 85) => ((2, 2), 2)
  (27, 28) => ((2, 1), 1)
  (5, 10)  => ((1, 1), 2)
  (28, 29) => ((2, 1), 2)
  (85, 91) => ((2, 2), 3)
```
"""
@inline function construct_boundary_edge_map(
        boundary_nodes::A,
        IntegerType::Type{I} = number_type(boundary_nodes),
        EdgeType::Type{E} = NTuple{2, IntegerType},
    ) where {A, I, E}
    if has_multiple_curves(boundary_nodes)
        return _construct_boundary_edge_map_multiple_curves(boundary_nodes, IntegerType, EdgeType)
    elseif has_multiple_sections(boundary_nodes)
        return _construct_boundary_edge_map_multiple_sections(boundary_nodes, IntegerType, EdgeType)
    else
        return _construct_boundary_edge_map_contiguous(boundary_nodes, IntegerType, EdgeType)
    end
end
function _construct_boundary_edge_map_multiple_curves(boundary_nodes, IntegerType::Type{I}, EdgeType::Type{E}) where {I, E}
    dict = Dict{E, Tuple{NTuple{2, I}, I}}()
    return _construct_boundary_edge_map!(dict, boundary_nodes)
end
function _construct_boundary_edge_map!(dict::Dict{E, Tuple{NTuple{2, I}, I}}, boundary_nodes) where {E, I}
    nc = num_curves(boundary_nodes)
    for m in 1:nc
        bn_m = get_boundary_nodes(boundary_nodes, m)
        ns = num_sections(bn_m)
        for n in 1:ns
            bn_n = get_boundary_nodes(bn_m, n)
            ne = num_boundary_edges(bn_n)
            for ℓ in 1:ne
                u = get_boundary_nodes(bn_n, ℓ)
                v = get_boundary_nodes(bn_n, ℓ + 1)
                e = construct_edge(E, u, v)
                dict[e] = ((m, n), ℓ)
            end
        end
    end
    return dict
end
function _construct_boundary_edge_map_multiple_sections(boundary_nodes, IntegerType::Type{I}, EdgeType::Type{E}) where {I, E}
    dict = Dict{E, NTuple{2, I}}()
    return _construct_boundary_edge_map!(dict, boundary_nodes)
end
function _construct_boundary_edge_map!(dict::Dict{E, NTuple{2, I}}, boundary_nodes) where {E, I} 
    ns = num_sections(boundary_nodes)
    for n in 1:ns
        bn_n = get_boundary_nodes(boundary_nodes, n)
        ne = num_boundary_edges(bn_n)
        for ℓ in 1:ne
            u = get_boundary_nodes(bn_n, ℓ)
            v = get_boundary_nodes(bn_n, ℓ + 1)
            e = construct_edge(E, u, v)
            dict[e] = (n, ℓ)
        end
    end
    return dict
end
function _construct_boundary_edge_map_contiguous(boundary_nodes::A, IntegerType::Type{I}, EdgeType::Type{E}) where {A, I, E}
    dict = Dict{E, Tuple{A, I}}()
    return _construct_boundary_edge_map!(dict, boundary_nodes)
end
function _construct_boundary_edge_map!(dict::Dict{E, Tuple{A, I}}, boundary_nodes) where {A, I, E}
    ne = num_boundary_edges(boundary_nodes)
    for ℓ in 1:ne
        u = get_boundary_nodes(boundary_nodes, ℓ)
        v = get_boundary_nodes(boundary_nodes, ℓ + 1)
        e = construct_edge(E, u, v)
        dict[e] = (boundary_nodes, ℓ)
    end
    return dict
end

function _bemcopy(dict::Dict{E, Tuple{A, I}}; boundary_nodes) where {A, I, E} # so that boundary_nodes remains aliased with this dict when copying in Triangulation
    new_dict = Dict{E, Tuple{A, I}}()
    for (e, (pos, ℓ)) in dict
        new_dict[e] = (boundary_nodes, ℓ)
    end
    return new_dict
end
_bemcopy(dict::Dict{E, NTuple{2, I}}; boundary_nodes) where {I, E} = copy(dict)
_bemcopy(dict::Dict{E, Tuple{NTuple{2, I}, I}}; boundary_nodes) where {E, I} = copy(dict)

"""
    insert_boundary_node!(boundary_nodes, pos, node)

Inserts a boundary node `node` into `boundary_nodes` at the position 
`pos`. Here, `pos[1]` is such that `get_boundary_nodes(boundary_nodes, pos[1])`
is the section that the node will be inserted onto, and `pos[2]` gives the position
of the array to insert `node` into. In particular, 

    insert_boundary_node!(boundary_nodes, pos, node)

is the same as 

    insert!(get_boundary_nodes(boundary_nodes, pos[1]), pos[2], node)

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> boundary_nodes = [1, 2, 3, 4, 5, 1];

julia> DelaunayTriangulation.insert_boundary_node!(boundary_nodes, (boundary_nodes, 4), 23)
7-element Vector{Int64}:
  1
  2
  3
 23
  4
  5
  1

julia> boundary_nodes = [[7, 13, 9, 25], [25, 26, 29, 7]];

julia> DelaunayTriangulation.insert_boundary_node!(boundary_nodes, (2, 1), 57)
2-element Vector{Vector{Int64}}:
 [7, 13, 9, 25]
 [57, 25, 26, 29, 7]

julia> boundary_nodes = [[[17, 23, 18, 25], [25, 26, 81, 91], [91, 101, 17]], [[1, 5, 9, 13], [13, 15, 1]]];

julia> DelaunayTriangulation.insert_boundary_node!(boundary_nodes, ((1, 3), 3), 1001)
2-element Vector{Vector{Vector{Int64}}}:
 [[17, 23, 18, 25], [25, 26, 81, 91], [91, 101, 1001, 17]]
 [[1, 5, 9, 13], [13, 15, 1]]
```
"""
function insert_boundary_node!(boundary_nodes, pos, node)
    nodes = get_boundary_nodes(boundary_nodes, pos[1])
    insert!(nodes, pos[2], node)
    return boundary_nodes
end

"""
    delete_boundary_node!(boundary_nodes, pos)

Deletes the boundary node at the position `pos` from `boundary_nodes`. Here, `pos[1]` is such that
`get_boundary_nodes(boundary_nodes, pos[1])` is the section that the node will be deleted from, and
`pos[2]` gives the position of the array to delete from. In particular, 

    delete_boundary_node!(boundary_nodes, pos)

is the same as

    deleteat!(get_boundary_nodes(boundary_nodes, pos[1]), pos[2])

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> boundary_nodes = [71, 25, 33, 44, 55, 10];

julia> DelaunayTriangulation.delete_boundary_node!(boundary_nodes, (boundary_nodes, 4))
5-element Vector{Int64}:
 71
 25
 33
 55
 10

julia> boundary_nodes = [[7, 13, 9, 25], [25, 26, 29, 7]];

julia> DelaunayTriangulation.delete_boundary_node!(boundary_nodes, (2, 3))
2-element Vector{Vector{Int64}}:
 [7, 13, 9, 25]
 [25, 26, 7]

julia> boundary_nodes = [[[17, 23, 18, 25], [25, 26, 81, 91], [91, 101, 17]], [[1, 5, 9, 13], [13, 15, 1]]];

julia> DelaunayTriangulation.delete_boundary_node!(boundary_nodes, ((2, 2), 2))
2-element Vector{Vector{Vector{Int64}}}:
 [[17, 23, 18, 25], [25, 26, 81, 91], [91, 101, 17]]
 [[1, 5, 9, 13], [13, 1]]
```
"""
function delete_boundary_node!(boundary_nodes, pos)
    nodes = get_boundary_nodes(boundary_nodes, pos[1])
    deleteat!(nodes, pos[2])
    return boundary_nodes
end

@doc """
    get_curve_index(dict, ghost_vertex) -> Int
    get_curve_index(ghost_vertex) -> Int

Given a `Dict` from [`construct_ghost_vertex_map`](@ref) and a `ghost_vertex`,
returns the index of the curve corresponding to that ghost vertex. The second method 
maps `ghost_vertex` to `1` if it is an `Integer` or a `Vector`, and `ghost_vertex[1]` if it is a `Tuple`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> DelaunayTriangulation.get_curve_index(-1)
1

julia> DelaunayTriangulation.get_curve_index((5, 3))
5

julia> gv_map = DelaunayTriangulation.construct_ghost_vertex_map([[[1, 5, 17, 18, 1]], [[23, 29, 31, 33], [33, 107, 101], [101, 99, 85, 23]]])
Dict{Int64, Tuple{Int64, Int64}} with 4 entries:
  -1 => (1, 1)
  -3 => (2, 2)
  -2 => (2, 1)
  -4 => (2, 3)

julia> DelaunayTriangulation.get_curve_index(gv_map, -1)
1

julia> DelaunayTriangulation.get_curve_index(gv_map, -2)
2

julia> DelaunayTriangulation.get_curve_index(gv_map, -3)
2

julia> DelaunayTriangulation.get_curve_index(gv_map, -4)
2
```
"""
get_curve_index
get_curve_index(ghost_vertex::Tuple) = ghost_vertex[1]
get_curve_index(ghost_vertex) = 1
get_curve_index(dict, ghost_vertex) = get_curve_index(dict[ghost_vertex])

@doc """
    get_section_index(dict, ghost_vertex) -> Int
    get_section_index(ghost_vertex) -> Int

Given a `Dict` from [`construct_ghost_vertex_map`](@ref) and a `ghost_vertex`,
returns the index of the section corresponding to that ghost vertex. The second method
maps `ghost_vertex` to itself if it is an `Integer`, `1` if it is a `Vector`, and 
`ghost_vertex[2]` if it is a `Tuple`. 

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> DelaunayTriangulation.get_section_index((2, 3)) # 3rd section of the 2nd curve
3

julia> DelaunayTriangulation.get_section_index(4)
4

julia> DelaunayTriangulation.get_section_index([1, 2, 3, 4, 5, 1])
1

julia> gv_map = DelaunayTriangulation.construct_ghost_vertex_map([[[1, 5, 17, 18, 1]], [[23, 29, 31, 33], [33, 107, 101], [101, 99, 85, 23]]])     
Dict{Int64, Tuple{Int64, Int64}} with 4 entries:
  -1 => (1, 1)
  -3 => (2, 2)
  -2 => (2, 1)
  -4 => (2, 3)

julia> DelaunayTriangulation.get_section_index(gv_map, -1)
1

julia> DelaunayTriangulation.get_section_index(gv_map, -2)
1

julia> DelaunayTriangulation.get_section_index(gv_map, -3)
2

julia> DelaunayTriangulation.get_section_index(gv_map, -4)
3
```
"""
get_section_index
get_section_index(ghost_vertex::Tuple) = ghost_vertex[2]
get_section_index(ghost_vertex::Integer) = ghost_vertex
get_section_index(ghost_vertex) = 1
get_section_index(dict, ghost_vertex) = get_section_index(dict[ghost_vertex])

"""
    construct_ghost_vertex_ranges(boundary_nodes::A, IntegerType::Type{I}=number_type(boundary_nodes)) where {A,I} -> Dict

Given a set of `boundary_nodes`, returns a `Dict` that maps ghost vertices to 
the range of all ghost vertices that the corresponding boundary curve could correspond to.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> boundary_nodes = [
                  [
                      [1, 2, 3, 4], [4, 5, 6, 1]
                  ],
                  [
                      [18, 19, 20, 25, 26, 30]
                  ],
                  [
                      [50, 51, 52, 53, 54, 55], [55, 56, 57, 58], [58, 101, 103, 105, 107, 120], [120, 121, 122, 50]
                  ]
              ]
3-element Vector{Vector{Vector{Int64}}}:
 [[1, 2, 3, 4], [4, 5, 6, 1]]
 [[18, 19, 20, 25, 26, 30]]
 [[50, 51, 52, 53, 54, 55], [55, 56, 57, 58], [58, 101, 103, 105, 107, 120], [120, 121, 122, 50]]

julia> DelaunayTriangulation.construct_ghost_vertex_ranges(boundary_nodes)
Dict{Int64, UnitRange{Int64}} with 7 entries:
  -5 => -7:-4
  -1 => -2:-1
  -7 => -7:-4
  -3 => -3:-3
  -2 => -2:-1
  -4 => -7:-4
  -6 => -7:-4
```
"""
@inline function construct_ghost_vertex_ranges(boundary_nodes, IntegerType::Type{I} = number_type(boundary_nodes)) where {I}
    if has_multiple_curves(boundary_nodes)
        return _construct_ghost_vertex_ranges_multiple_curves(boundary_nodes, IntegerType)
    elseif has_multiple_sections(boundary_nodes)
        return _construct_ghost_vertex_ranges_multiple_sections(boundary_nodes, IntegerType)
    else
        return _construct_ghost_vertex_ranges_contiguous(IntegerType)
    end
end
function _construct_ghost_vertex_ranges_multiple_curves(boundary_nodes, IntegerType::Type{I}) where {I}
    start = I(𝒢)
    current_ghost_vertex = I(𝒢)
    dict = Dict{I, UnitRange{I}}()
    nc = num_curves(boundary_nodes)
    for i in 1:nc
        bn = get_boundary_nodes(boundary_nodes, i)
        ns = num_sections(bn)
        stop = start - ns + 1
        for _ in 1:ns
            dict[current_ghost_vertex] = stop:start
            current_ghost_vertex -= 1
        end
        start -= ns
    end
    return dict
end
function _construct_ghost_vertex_ranges_multiple_sections(boundary_nodes, IntegerType::Type{I}) where {I}
    current_ghost_vertex = I(𝒢)
    dict = Dict{I, UnitRange{I}}()
    ns = num_sections(boundary_nodes)
    range = (current_ghost_vertex - ns + 1):current_ghost_vertex
    for _ in 1:ns
        dict[current_ghost_vertex] = range
        current_ghost_vertex -= 1
    end
    return dict
end
function _construct_ghost_vertex_ranges_contiguous(IntegerType::Type{I}) where {I}
    current_ghost_vertex = I(𝒢)
    dict = Dict{I, UnitRange{I}}()
    dict[current_ghost_vertex] = current_ghost_vertex:current_ghost_vertex
    return dict
end

"""
    get_skeleton(boundary_nodes, IntegerType) -> empty(boundary_nodes)

Given a set of boundary nodes, returns the empty skeleton of the boundary nodes. This is essentially `empty` applied to `boundary_nodes`,
but with vertices of type `IntegerType`. This is mainly needed for [`convert_boundary_curves!`](@ref). You will need to implement a new method for this 
if you want your custom boundary node interface to be supported for curve-bounded domains.
"""
@inline function get_skeleton(boundary_nodes, ::Type{I}) where {I}
    if has_multiple_curves(boundary_nodes)
        return _get_skeleton_multiple_curves(boundary_nodes, I)
    elseif has_multiple_sections(boundary_nodes)
        return _get_skeleton_multiple_sections(boundary_nodes, I)
    else
        return _get_skeleton_contiguous(boundary_nodes, I)
    end
end
@inline function _get_skeleton_multiple_curves(boundary_nodes, ::Type{I}) where {I}
    boundary_nodes′ = Vector{Vector{Vector{I}}}(undef, num_curves(boundary_nodes))
    for curve_index in 1:num_curves(boundary_nodes)
        boundary_nodes_curve = get_boundary_nodes(boundary_nodes, curve_index)
        boundary_nodes′[curve_index] = _get_skeleton_multiple_sections(boundary_nodes_curve, I)
    end
    return boundary_nodes′
end
@inline function _get_skeleton_multiple_sections(boundary_nodes, ::Type{I}) where {I}
    boundary_nodes′ = Vector{Vector{I}}(undef, num_sections(boundary_nodes))
    for section_index in 1:num_sections(boundary_nodes)
        boundary_nodes_section = get_boundary_nodes(boundary_nodes, section_index)
        boundary_nodes′[section_index] = _get_skeleton_contiguous(boundary_nodes_section, I)
    end
    return boundary_nodes′
end
@inline function _get_skeleton_contiguous(boundary_nodes, ::Type{I}) where {I}
    return I[]
end

"""
    set_boundary_node!(boundary_nodes, pos, node)

Given a set of `boundary_nodes`, sets the boundary node at position `pos` to `node`.
Here, `pos[1]` is such that `get_boundary_nodes(boundary_nodes, pos[1])`
is the section that the node will be set onto, and `pos[2]` gives the position
of the array to set `node` into. In particular, 

    set_boundary_node!(boundary_nodes, pos, node)

is the same as

    get_boundary_nodes(boundary_nodes, pos[1])[pos[2]] = node

assuming `setindex!` is defined for the type of `boundary_nodes`.
"""
function set_boundary_node!(boundary_nodes, pos, node)
    nodes = get_boundary_nodes(boundary_nodes, pos[1])
    nodes[pos[2]] = node
    return boundary_nodes
end