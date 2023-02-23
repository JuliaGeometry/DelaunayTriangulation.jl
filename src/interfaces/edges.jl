"""
    construct_edge(::Type{E}, i, j) where {E}

Constructs an edge with indices `(i, j)` with the type `E`. The 
following methods are currently defined:

    construct_edge(::Type{NTuple{2, I}}, i, j) where {I}
    construct_edge(::Type{A}, i, j) where {I,A<:AbstractVector} 

You can extend this function as you need, ensuring you define it 
on the type rather than on an instance of the type. 
"""
function construct_edge end
function construct_edge(::Type{E}, i, j) where {E}
    return error("The construct_edge function has not been defined for the type $E.")
end
construct_edge(::Type{NTuple{2,I}}, i, j) where {I} = (I(i), I(j))
construct_edge(::Type{A}, i, j) where {I,A<:AbstractVector{I}} = A([i, j])

"""
    initial(e::E)

Given an edge `e`, returns the index for the initial point. The following 
methods are currently defined:

    initial(e::NTuple{2,I}) where {I}
    initial(e::A) where {I,A<:AbstractVector{I}} 

You can extend this function as you need.
"""
function initial end
function initial(::E) where {E}
    return error("The initial function has not been defined for the type $E.")
end
initial(e::NTuple{2,I}) where {I} = e[1]
initial(e::A) where {I,A<:AbstractVector{I}} = e[begin]

"""
    terminal(e::E)

Given an edge `e`, returns the index for the terminal point. The following 
methods are currently defined:

    terminal(e::NTuple{2,I}) where {I}
    terminal(e::A) where {I,A<:AbstractVector{I}} 

You can extend this function as you need.
"""
function terminal end
function terminal(::E) where {E}
    return error("The terminal function has not been defined for the type $E.")
end
terminal(e::NTuple{2,I}) where {I} = e[2]
terminal(e::A) where {I,A<:AbstractVector{I}} = e[end]

"""
    edge_indices(e)

Given an edge `e`, returns `(initial(e), terminal(e))`. 

See also [`initial`](@ref) and [`terminal`](@ref).
"""
function edge_indices end
edge_indices(e) = (initial(e), terminal(e))

"""
    initialise_edges(::Type{S})

For a given type `S` for some collection (e.g. a `Set`), returns an
empty instance of that collection. The only method defined is

    initialise_edges(::Type{S}) where {E, S <: Set{E}}
    initialise_edges(::Type{A}) where {E, A <: AbstractVector{E}}

which returns a `Set{E}()` or a `A()`, respectively. You can extend this 
function as you need, making sure you extend it for the type rather than 
for instances of that type.
"""
function initialise_edges end
function initialise_edges(::Type{S}) where {S}
    return error("The initialise_edges function has not been defined for the type $S.")
end
initialise_edges(::Type{S}) where {E,S<:Set{E}} = S()
initialise_edges(::Type{A}) where {E,A<:AbstractVector{E}} = A()

"""
    edge_type(::Type{S}) where {S}

For a given type `S` representing a collection of edges, 
returns the type of triangle used inside `S`, e.g. `NTuple{2, Int64}`
if `S = Set{NTuple{2, Int64}}`. The only methods defined are 

    edge_type(::Type{S}) where {E,S<:Set{E}}
    edge_type(::Type{A}) where {E,A<:AbstractVector{E}} 

which return `E`. You can extend this function as you need, making sure 
you extend it for the type rather than for instances of that type.
"""
function edge_type end
function edge_type(::Type{F}) where {F}
    return error("The edge_type function has not been defined for the type $F.")
end
edge_type(::Type{S}) where {E,S<:Set{E}} = E
edge_type(::Type{A}) where {E,A<:AbstractVector{E}} = E

"""
    num_edges(E::S) where {S}

Given a collection of edges `E`, returns the number of edges
in `E`. The only method currently defined is 

    num_edges(E::Set)

which returns `length(E)`. You can extend this function as you need.
"""
function num_edges end
function num_edges(::S) where {S}
    return error("The num_edges function has not been defined for the type $S.")
end
num_edges(E::Set) = length(E)

"""
    contains_edge(e::E, Es::S) where {E, S}

Given a collection of edges `Es` of type `S`, containing edges
of type `E`, checks if `Es` includes the edge `e`, returning `true` 
if so. The methods currently defined are

    contains_edge(e::E, Es::Set{E}) where {E} 
    contains_edge(e::E, Es::A) where {E,A<:AbstractVector{E}}
    contains_edge(i, j, Es::E) 

The first two methods simply return `e ∈ E`, while the latter
constructs the edge `e = (i, j)` of type `edge_type(E)` and call the first two 
methods. Only the method `contains_edge(::E, ::Es)` needs to be extended if you need, 
the last method makes use of this definition.
"""
function contains_edge end
function contains_edge(::Any, ::F) where {F}
    return error("The contains_edge function has not been defined for the type $F.")
end
contains_edge(e::E, Es::Set{E}) where {E} = e ∈ Es
contains_edge(e::E, Es::A) where {E,A<:AbstractVector{E}} = e ∈ Es
contains_edge(i, j, Es::E) where {E} = contains_edge(construct_edge(edge_type(E), i, j), Es)

"""
    add_to_edges!(E::S, e) where {S}

Given a collection of edges `E`, pushes `e` into it. The only 
methods currently defined are

    add_to_edges!(E::Set, e)
    add_to_edges!(E::Vector, e)

which simply call `push!(E, e)`. You can extend this function  
as you need. 
"""
function add_to_edges! end
function add_to_edges!(::S, e) where {S}
    return error("The add_to_edges! function has not been defined for the type $S.")
end
function add_to_edges!(E::Set, e)
    push!(E, e)
    return nothing
end
function add_to_edges!(E::Vector, e)
    push!(E, e)
    return nothing 
end

"""
    add_edge!(E, e...)

Given a collection of edges `E`, adds all the triangles `e...` into it. 
To extend this method to other collections, see [`add_to_edges!`](@ref).
"""
function add_edge!(E, e::Vararg{F,N}) where {F,N}
    for i in 1:N
        add_to_edges!(E, e[i])
    end
    return nothing
end

"""
    delete_from_edges!(E::S, e::F) where {S}

Given a collection of edges `E` of type `S`, containing 
edges of type `F`, deletes the edge `e` from `E`.  The 
methods currently defined are

    delete_from_edges!(E::Set{F}, T::F) where {F}
    delete_from_edges!(Es::A, e::E) where {E, A<:AbstractVector{E}}
    
which just calls `delete!` on `E` in the first case, or `filter!` 
in the second case. This is the form of the function that needs 
to be extended. We also define 

    delete_from_edges!(Es::E, i::Integer, j::Integer) where {E}

which constructs the edge `(i, j)` and then deletes it from `Es`, calling 
the methods above. You do not need to extend this last method.
"""
function delete_from_edges! end
function delete_from_edges!(::S, e) where {S}
    return error("The delete_from_edges! function has not been defined for the type $S.")
end
function delete_from_edges!(E::Set{F}, e::F) where {F}
    has_e = contains_edge(e, E)
    if has_e
        delete!(E, e)
    end
    return nothing
end
function delete_from_edges!(Es::A, e::E) where {E,A<:AbstractVector{E}}
    has_e = contains_edge(e, Es)
    if has_e
        filter!(E -> e ≠ E, Es)
    end
    return nothing
end
function delete_from_edges!(Es::E, i::Integer, j::Integer) where {E}
    V = edge_type(E)
    e = construct_edge(V, i, j)
    delete_from_edges!(Es, e)
    return nothing
end

"""
    delete_edge!(E, e...)

Given a collection of edges `E`, deletes all the edges `e...` from it. 
To extend this method to other collections, see [`delete_from_edges!`](@ref).
"""
function delete_edge!(E, e::Vararg{F,N}) where {F,N}
    for i in 1:N
        delete_from_edges!(E, e[i])
    end
    return nothing
end

"""
    each_edge(E::F) where {F}

For a given collection of edges `E`, returns an iterator that 
goes over each edge in the collection. The methods currently 
defined are 

    each_edge(E::Set)
    each_edge(E::AbstractMatrix)

with the first method simply returning `E`, and the second returning 
`eachcol(E)`. You can extend this function as you need.
"""
function each_edge end
function each_edge(::F) where {F}
    return error("The each_edge function has not been defined for the type $F.")
end
each_edge(E::Set) = E
each_edge(E::AbstractMatrix) = eachcol(E)

"""
    random_edge(E)

Returns a random edge from edge set `E`.
"""
function random_edge end
random_edge(E, rng::AbstractRNG=Random.default_rng()) = rand(rng, each_edge(E))

"""
    is_empty(E)

Tests if the edge set `E` is empty.
"""
function is_empty end
is_empty(E) = isempty(E)
