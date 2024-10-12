geti(T) = T[1]

getj(T) = T[2]

getk(T) = T[3] 

triangle_vertices(T) = (geti(T), getj(T), getk(T))

triangle_edges(i, j, k) = ((i, j), (j, k), (k, i))
triangle_edges(T) = triangle_edges(geti(T), getj(T), getk(T))

function sort_triangle(T) 
    i, j, k = triangle_vertices(T)
    minijk = min(i, j, k)
    return minijk == i ? (j, k, i) : minijk == j ? (k, i, j) : (i, j, k)
end
sort_triangle(i::Integer, j::Integer, k::Integer) = sort_triangle((i, j, k))
@inline function sort_triangle(i, j, k) # in case one of the vertices is a point rather than an integer 
    # Note that the return type is a 3-Union 
    return is_ghost_vertex(i) ? (j, k, i) : is_ghost_vertex(j) ? (k, i, j) : (i, j, k)
end

function compare_triangles(T, V)
    i, j, k = triangle_vertices(T)
    u, v, w = triangle_vertices(V)
    return (i, j, k) == (u, v, w) ||
        (i, j, k) == (v, w, u) ||
        (i, j, k) == (w, u, v)
end

function construct_positively_oriented_triangle(i, j, k, points, predicates::AbstractPredicateKernel = AdaptiveKernel())
    p, q, r = get_point(points, i, j, k)
    orientation = triangle_orientation(predicates, p, q, r)
    if is_negatively_oriented(orientation)
        return (j, i, k)
    else
        return (i, j, k)
    end
end

function rotate_triangle(T, ::Val{N}) where {N}
    i, j, k = triangle_vertices(T)
    N < 0 && throw(ArgumentError(lazy"Cannot rotate triangle $T by a negative amount."))
    return N == 0 ? T : N == 1 ? (j, k, i) : N == 2 ? (k, i, j) : rotate_triangle(T, Val(N % 3))
end
rotate_triangle(T, N::Integer) = rotate_triangle(Val(N))

struct Triangles{I,E} <: AbstractSet{NTuple{3,I}}
    adjacent::Adjacent{I,E}
end
Triangles{I,E}() where {I,E} = Triangle{I,E}(Adjacent{I,E}())
get_adjacent(triangles::Triangles) = triangles.adjacent
get_dict(triangles::Triangles) = get_adjacent(get_adjacent(triangles))
num_triangles(triangles) = length(triangles)
contains_triangle(T, triangles) = in(sort_triangle(T), triangles) # not doing ::Triangles because of an ambiguity
add_triangle!(triangles::Triangles, T...) = add_triangle!(get_adjacent(triangles), T...)
delete_triangle!(triangles::Triangles, T) = delete_triangle!(get_adjacent(triangles), T)
each_triangle(triangles) = triangles 
function compare_triangle_collections(T::Triangles, V::Triangles)
    num_triangles(T) == num_triangles(V) || return false
    for τ in each_triangle(T)
        contains_triangle(τ, V) || return false 
    end
    return true
end
integer_type(::Triangles{I}) where {I} = I
edge_type(::Triangles{I,E}) where {I,E} = E

# AbstractSet
Base.length(triangles::Triangles) = length(get_dict(triangles)) ÷ 3
function Base.iterate(triangles::Triangles, state...)
    ew_state = Base.iterate(get_dict(triangles), state...)
    ew_state === nothing && return nothing
    (e, w), state = ew_state
    u, v = edge_vertices(e)
    T = (u, v, w)
    while T != sort_triangle(T) # unique representation
        ew_state = Base.iterate(get_dict(triangles), state...)
        ew_state === nothing && return nothing
        (e, w), state = ew_state
        u, v = edge_vertices(e)
        T = (u, v, w)
    end
    return T, state
end
Base.empty(triangles::Triangles{I,E}) where {I, E} = Triangles{I,E}()
Base.isempty(triangles::Triangles{I,E}) where {I, E} = isempty(get_dict(triangles))
function Base.in(T, triangles::Triangles)
    u, v, w = triangle_vertices(T)
    adj = get_adjacent(triangles)
    return edge_exists(adj, u, v) && edge_exists(adj, v, w) && edge_exists(adj, w, u)
end
Base.copy(triangles::Triangles) = Triangles(copy(get_adjacent(triangles)))
Base.empty!(triangles::Triangles) = empty!(get_adjacent(triangles))

concretize_triangles(triangles) = Set(each_triangle(triangles))