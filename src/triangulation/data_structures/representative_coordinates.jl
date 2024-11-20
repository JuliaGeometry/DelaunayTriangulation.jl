mutable struct RepresentativeCoordinates{T}
    x::T 
    y::T 
    n::Int 
end
RepresentativeCoordinates{T}() where {T} = RepresentativeCoordinates{T}(zero(T), zero(T), 0)

getx(c::RepresentativeCoordinates) = c.x
gety(c::RepresentativeCoordinates) = c.y
getn(c::RepresentativeCoordinates) = c.n

function Base.:(==)(p::RepresentativeCoordinates, q::RepresentativeCoordinates)
    return getxy(p) == getxy(q) 
end 

function Base.convert(::Type{RepresentativeCoordinates{T}}, c::RepresentativeCoordinates) where T
    x, y = getxy(c)
    n = getn(c)
    return RepresentativeCoordinates(T(x), T(y), n)
end
function Base.convert(::Type{RepresentativeCoordinates{T}}, c::RepresentativeCoordinates{T}) where T
    return c
end

Base.copy(c::RepresentativeCoordinates) = RepresentativeCoordinates(getx(c), gety(c), getn(c))

function reset!(c::RepresentativeCoordinates{T}) where T
    c.x = zero(T)
    c.y = zero(T)
    c.n = 0
    return c
end

function add_point!(c::RepresentativeCoordinates, p)
    n = getn(c)
    q = getxy(c)
    x, y = 1 / (n + 1) .* (n .* q .+ p)
    c.x = x
    c.y = y
    c.n += 1
    return c
end

function delete_point!(c::RepresentativeCoordinates, p)
    n = getn(c)
    q = getxy(c)
    x, y = 1 / (n - 1) .* (n .* q .- p)
    c.x = x
    c.y = y
    c.n -= 1
    return c
end

function compute_centroid!(c::RepresentativeCoordinates, points)
    reset!(c)
    foreach(each_point(points)) do p
        add_point!(c, p)
    end
end

struct RepresentativeCoordinatesList{T}
    coordinates::Vector{RepresentativeCoordinates{T}}
end
@inline get_coordinates(list::RepresentativeCoordinatesList) = list.coordinates
@inline Base.getindex(list::RepresentativeCoordinatesList, i::Int) = get_coordinates(list)[i]
@inline Base.copy(list::RepresentativeCoordinatesList) = RepresentativeCoordinatesList(copy(get_coordinates(list)))
@inline Base.length(list::RepresentativeCoordinatesList) = length(get_coordinates(list))
@inline function Base.resize!(list::RepresentativeCoordinatesList{T}, n::Int) where {T}
    coordinates = get_coordinates(list)
    if n < length(coordinates)
        resize!(coordinates, n)
    else
        for _ in length(coordinates):n
            push!(coordinates, RepresentativeCoordinates{T}())
        end
    end
    return list
end
    
function reset!(list::RepresentativeCoordinatesList)
    foreach(get_coordinates(list)) do c 
        reset!(c)
    end
    return list
end

function Base.empty!(list::RepresentativeCoordinatesList)
    coordinates = get_coordinates(list)
    empty!(coordinates)
    return list
end

function update_after_addition!(list::RepresentativeCoordinatesList{T}, curve_index, p) where {T}
    curve_index > length(list) && resize!(list, curve_index)
    q = list[curve_index]
    add_point!(q, p)
    return list
end

function update_after_deletion!(list::RepresentativeCoordinatesList{T}, curve_index, p) where {T}
    curve_index > length(list) && resize!(list, curve_index)
    q = list[curve_index]
    delete_point!(q, p)
    return list
end