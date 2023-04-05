const AV = AbstractVector

"""
    has_multiple_curves(bn::A) where {A}

Returns `true` if the given set of boundary nodes `bn` defines multiple curves, 
meaning disjoint boundary curves. We currently define the methods 

    has_multiple_curves(::AAA) where {F<:Number,A<:AV{F},AA<:AV{A},AAA<:AV{AA}} 
    has_multiple_curves(::AA) where {F<:Number,A<:AV{F},AA<:AV{A}} 
    has_multiple_curves(::A) where {F<:Number,A<:AV{F}}

with the first method returning `true`, while the last two methods return `false`,
and `AV = AbstractVector`. You can extend this function as you need.
"""
function has_multiple_curves end
function has_multiple_curves(::F) where {F}
    return error("The has_multiple_curves function has not been defined for the type $F.")
end
function has_multiple_curves(::AAA) where {F<:Number,A<:AV{F},AA<:AV{A},
    AAA<:AV{AA}}
    return true
end
has_multiple_curves(::AA) where {F<:Number,A<:AV{F},AA<:AV{A}} = false
has_multiple_curves(::A) where {F<:Number,A<:AV{F}} = false

"""
    has_multiple_segments(bn::A) where {A}

Returns `true` if the given set of boundary nodes `bn` contains multiple segments, 
meaning disjoint boundary curves. We currently define the methods 

    has_multiple_segments(::AAA) where {F<:Number,A<:AV{F},AA<:AV{A},AAA<:AV{AA}} 
    has_multiple_segments(::AA) where {F<:Number,A<:AV{F},AA<:AV{A}} 
    has_multiple_segments(::A) where {F<:Number,A<:AV{F}}

with the first and second methods returning `true`, while the last method returns `false`,
and `AV = AbstractVector`. You can extend this function as you need.
"""
function has_multiple_segments end
function has_multiple_segments(::F) where {F}
    return error("The has_multiple_segments function has not been defined for the type $F.")
end
function has_multiple_segments(::AAA) where {F<:Number,A<:AV{F},AA<:AV{A},
    AAA<:AV{AA}}
    return true
end
has_multiple_segments(::AA) where {F<:Number,A<:AV{F},AA<:AV{A}} = true
has_multiple_segments(::A) where {F<:Number,A<:AV{F}} = false

"""
    num_curves(bn::A)

Returns the number of `curves` defined by the boundary nodes `bn`. We currently 
define the methods

    num_curves(bn::AAA) where {F<:Number,A<:AV{F},AA<:AV{A},AAA<:AV{AA}} 

which simply returns `bn`.
"""
function num_curves end
function num_curves(::F) where {F}
    return error("The num_curves function has not been defined for the type $F.")
end
num_curves(bn::AAA) where {F<:Number,A<:AV{F},AA<:AV{A},AAA<:AV{AA}} = length(bn)
num_curves(bn::AA) where {F<:Number,A<:AV{F},AA<:AV{A}} = 1
num_curves(bn::A) where {F<:Number,A<:AV{F}} = 1

"""
    num_segments(bn::A)

Returns the number of `segments` defined by the boundary nodes `bn`. We currently 
define the method

    num_segments(bn::AA) where {F<:Number,A<:AV{F},AA<:AV{A}}

which simply returns `length(bn)`.
"""
function num_segments end
function num_segments(::F) where {F}
    return error("The num_segments function has not been defined for the type $F.")
end
num_segments(bn::AA) where {F<:Number,A<:AV{F},AA<:AV{A}} = length(bn)

"""
    num_boundary_edges(bn)

Given a collection of boundary nodes `bn`, returns the number of edges. This only 
needs to be defined for individual segments. We define the method 

    num_boundary_edges(bn::A) where {A<:AbstractVector}

which returns `length(bn) -1` (`-1` because it is assumed that `bn[begin] == bn[end]`). 
This is the only method that needs to be extended.

See also [`getboundarynodes`](@ref).
"""
function num_boundary_edges end
function num_boundary_edges(::F) where {F}
    return error("The num_boundary_edges function has not been defined for the type $F.")
end
num_boundary_edges(bn::A) where {A<:AbstractVector} = max(length(bn) - 1, 0)

"""
    getboundarynodes(bn::A, mnℓ...)

Given a collection of boundary nodes `bn`, returns the specified component of the 
collection. There are several forms for the methods. In these methods, it is assumed 
that one-based indexing is used for accessing all the boundary nodes. If you want to 
use offsets, for example, then define `getboundarynodes` appropriately (e.g. maybe 
`getboundarynodes(bn, m)` could map `m` to `m-4` if `4` is your offset).

The methods that you need to define are those that go down a level, i.e. from a set of curves 
to a curve, from a set of segments to a set of nodes, and from a set of nodes to a node. Of course, 
if you only ever use e.g. a set of nodes, then you need only define that method. The methods that 
we define for this are

    getboundarynodes(bn::AAA, m::Integer) where {F<:Number,A<:AV{F},AA<:AV{A},AAA<:AV{AA}} 
    getboundarynodes(bn::AA, n::Integer) where {F<:Number,A<:AV{F},AA<:AV{A}} 
    getboundarynodes(bn::A, ℓ::Integer) where {F<:Number,A<:AV{F}} 

The first method takes a set of curves to the `m`th curve, the second takes a set of segments to the 
`n`th segment, and the third takes a set of nodes to the `ℓ`th node. These are the only methods 
that need to be extended. For the set of curves case, we also define

    getboundarynodes(bn, m::Integer, n::Integer)
    getboundarynodes(bn, (m, n)::NTuple{2,Integer})

which calls `getboundarynodes(getboundarynodes(bn, m), n)`. This does not need to be extended. Lastly, 
we also define 

    getboundarynodes(bn::A, ::A) where {A}

which simply returns `bn`. This is useful when using the result of [`construct_boundary_map`](@ref).
"""
function getboundarynodes end
function getboundarynodes(::F, ::Any) where {F}
    return error("The getboundarynodes function has not been defined for the type $F.")
end
function getboundarynodes(bn::AAA,
    m::Integer) where {F<:Number,A<:AV{F},AA<:AV{A},
    AAA<:AV{AA}}
    return bn[m]
end
getboundarynodes(bn::AA, n::Integer) where {F<:Number,A<:AV{F},AA<:AV{A}} = bn[n]
getboundarynodes(bn::A, ℓ::Integer) where {F<:Number,A<:AV{F}} = bn[ℓ]
getboundarynodes(bn, m::Integer, n::Integer) = getboundarynodes(getboundarynodes(bn, m), n)
getboundarynodes(bn, (m, n)::NTuple{2,Integer}) = getboundarynodes(bn, m, n)
getboundarynodes(bn::A, ::A) where {A} = bn # This is for the case where construct_boundary_map's result just returns bn 

"""
    get_boundary_nodes(bn, mnℓ...)

Get the boundary nodes from `bn` corresponding to the specified indices. 
See [`getboundarynodes`](@ref).
"""
get_boundary_nodes(bn, mnℓ...) = getboundarynodes(bn, mnℓ...)

"""
    each_boundary_node(bn::A)

Returns an iterator that goes over each node in `bn`. Only defined for 
single segments so that `bn` acts like a vector of numbers. The only method 
currently defined is 

    each_boundary_node(bn::A) where {F<:Number,A<:AbstractVector{F}}
    
which just returns `bn`. You can extend this function as you need. If you really 
want to loop over every boundary node, you can make use of the result from 
[`construct_boundary_map`](@ref).
"""
function each_boundary_node end
function each_boundary_node(::F) where {F}
    return error("The each_boundary_node function has not been defined for the type $F.")
end
each_boundary_node(bn::A) where {F<:Number,A<:AV{F}} = bn

"""
    construct_boundary_map(bn; IntegerType::Type{I} = Int64) where {I}

Given a set of boundary nodes `bn`, returns a `OrderedDict` that maps boundary indices 
to their position in `bn`. In particular:

- `has_multiple_curves(bn)`

In this case, the result is a `dict = OrderedDict{I, NTuple{2, I}}`. The results will be of the form 
`dict[i] = (m, n)`, so that boundary indices with value `i` correspond to nodes at 
`get_boundary_nodes(bn, m, n)`, i.e. the `n`th segment of the `m`th curve.
- `has_multiple_segments(bn)`

In this case, the result is a `dict = OrderedDict{I, I}`. The results will be of the form `dict[i] = n`,
so that boundary indices with value `i` correspond to nodes at `get_boundary_nodes(bn, n)`, i.e. 
the `n`th segment.

- `else`

Here, the result is a `dict = OrderedDict{I, F}`, mapping `$BoundaryIndex` back to `bn` and `F = typeof(bn)`.

# Iteration Tips 
This dict can be useful for iterating over all boundary nodes. For example, you could do

```julia 
bn_map = construct_boundary_map(bn)
for segment_index in values(bn_map)
    nodes = get_boundary_nodes(bn, segment_index)
    ## Do something with the nodes 
end 
```

The above will work for any form of `bn` also.
"""
Base.@constprop :aggressive function construct_boundary_map(bn;
    IntegerType::Type{I}=Int64) where {I}
    if has_multiple_curves(bn)
        dict = OrderedDict{I,NTuple{2,I}}()
        nc = num_curves(bn)
        current_idx = I(BoundaryIndex)
        for m in 1:nc
            bn_m = get_boundary_nodes(bn, m)
            ns = num_segments(bn_m)
            for n in 1:ns
                dict[current_idx] = (m, n)
                current_idx -= 1
            end
        end
    elseif has_multiple_segments(bn)
        dict = OrderedDict{I,I}()
        ns = num_segments(bn)
        current_idx = I(BoundaryIndex)
        for n in 1:ns
            dict[current_idx] = n
            current_idx -= 1
        end
    else
        dict = OrderedDict(I(BoundaryIndex) => bn)
    end
    return dict
end

function construct_boundary_edge_map(bn::A; IntegerType::Type{I}=Int64, EdgeType::Type{E}=NTuple{2,IntegerType}) where {A,I,E}
    if has_multiple_curves(bn)
        dict = Dict{E,Tuple{NTuple{2,I},I}}()
        nc = num_curves(bn)
        for m in 1:nc
            bn_m = get_boundary_nodes(bn, m)
            ns = num_segments(bn_m)
            for n in 1:ns
                bn_n = get_boundary_nodes(bn_m, n)
                ne = num_boundary_edges(bn_n)
                for ℓ in 1:ne
                    u = get_boundary_nodes(bn_n, ℓ)
                    v = get_boundary_nodes(bn_n, ℓ + 1)
                    dict[(u, v)] = ((m, n), ℓ)
                end
            end
        end
    elseif has_multiple_segments(bn)
        dict = Dict{E,Tuple{I,I}}()
        ns = num_segments(bn)
        for n in 1:ns
            bn_n = get_boundary_nodes(bn, n)
            ne = num_boundary_edges(bn_n)
            for ℓ in 1:ne
                u = get_boundary_nodes(bn_n, ℓ)
                v = get_boundary_nodes(bn_n, ℓ + 1)
                dict[(u, v)] = (n, ℓ)
            end
        end
    else
        dict = Dict{E,Tuple{A,I}}()
        ne = num_boundary_edges(bn)
        for ℓ in 1:ne
            u = get_boundary_nodes(bn, ℓ)
            v = get_boundary_nodes(bn, ℓ + 1)
            dict[(u, v)] = (bn, ℓ)
        end
    end
    return dict
end

"""
    has_multiple_segments(bm::OrderedDict)

Given a map from [`construct_boundary_map`](@ref), tests if `length(bm) > 1`, i.e 
if there are multiple segments in the corresponding boundary.
"""
has_multiple_segments(bm::OrderedDict) = length(bm) > 1

"""
    map_boundary_index(dict, i)

Given a `dict` from [`construct_boundary_map`](@ref), returns `dict[i]`. Also works 
for a `dict` from [`construct_boundary_index_ranges`](@ref).
"""
map_boundary_index(dict, i) = dict[i]

"""
    get_curve_index(dict, i)
    get_curve_index(i)

Given a `dict` from [`construct_boundary_map`](@ref) and a boundary index `i`, 
returns the index of the curve corresponding to that boundary index. The 
second method maps `i` to `1` if it is an integer, and `i[1]` if it is a `Tuple`.
"""
function get_curve_index end
get_curve_index(i::Tuple) = i[1]
get_curve_index(i) = 1
get_curve_index(dict, i) = get_curve_index(map_boundary_index(dict, i))

"""
    get_segment_index(dict, i)
    get_segment_index(i)

Given a `dict` from [`construct_boundary_map`](@ref) and a boundary index `i`, 
returns the index of the segment corresponding to that boundary index. The 
second method maps `i` to `i` if it is an integer, `1` if it is a vector, 
and `i[2]` if it is a `Tuple`.
"""
function get_segment_index end
get_segment_index(i::Tuple) = i[2]
get_segment_index(i::Integer) = i
get_segment_index(i) = 1
get_segment_index(dict, i) = get_segment_index(map_boundary_index(dict, i))

"""
    num_outer_boundary_segments(boundary_nodes)

Given a set of `boundary_nodes`, returns the number of segments 
that correspond to the outer boundary. Note that this also gives 
the range of outer boundary indices, i.e. 
`$BoundaryIndex:-1:-num_outer_boundary_segments(boundary_nodes)`.
"""
function num_outer_boundary_segments end
function num_outer_boundary_segments(boundary_nodes)
    if has_multiple_curves(boundary_nodes)
        bn = get_boundary_nodes(boundary_nodes, 1)
        ns = num_segments(bn)
        return ns
    elseif has_multiple_segments(boundary_nodes)
        ns = num_segments(boundary_nodes)
        return ns
    else
        return 1
    end
end

"""
    construct_boundary_index_ranges(boundary_nodes; IntegerType::Type{I}=Int64) where {I}

Given a set of `boundary_nodes`, creates an `OrderedDict` that maps boundary indices 
to the range of all boundary indices that the corresponding boundary curve could 
correspond to. For example, suppose we have 

```julia-repl 
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
```

Then the first curve, `[[1, 2, 3, 4], [4, 5, 6, 1]]` has boundary indices `$BoundaryIndex` and `$(BoundaryIndex-1)`, 
so the range would be `$(BoundaryIndex-1):$BoundaryIndex`. The full `Dict` we obtain will be 

```julia-repl
julia> construct_boundary_index_ranges(boundary_nodes)
OrderedDict{Int64, UnitRange{Int64}} with 7 entries:
  -1 => -2:-1
  -2 => -2:-1
  -3 => -3:-3
  -4 => -7:-4
  -5 => -7:-4
  -6 => -7:-4
  -7 => -7:-4
```
"""
function construct_boundary_index_ranges(boundary_nodes;
    IntegerType::Type{I}=Int64) where {I}
    start = I(BoundaryIndex)
    current_boundary_index = I(BoundaryIndex)
    dict = OrderedDict{I,UnitRange{I}}()
    if has_multiple_curves(boundary_nodes)
        nc = num_curves(boundary_nodes)
        for i in 1:nc
            bn = get_boundary_nodes(boundary_nodes, i)
            ns = num_segments(bn)
            stop = start - ns + 1
            for _ in 1:ns
                dict[current_boundary_index] = stop:start
                current_boundary_index -= 1
            end
            start -= ns
        end
    elseif has_multiple_segments(boundary_nodes)
        ns = num_segments(boundary_nodes)
        range = (current_boundary_index-ns+1):current_boundary_index
        for _ in 1:ns
            dict[current_boundary_index] = range
            current_boundary_index -= 1
        end
    else
        dict[current_boundary_index] = current_boundary_index:current_boundary_index
    end
    return dict
end