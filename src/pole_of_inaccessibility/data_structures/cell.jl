struct Cell{T}
    x::T 
    y::T 
    half_width::T 
    dist::T 
    max_dist::T 
    function Cell(x::T, y::T, half_width::T, points, boundary_nodes) where {T}
        q = (x, y)
        dist = distance_to_polygon(q, points, boundary_nodes)
        max_dist = dist + half_width * sqrt(2)
        return new{T}(x, y, half_width, dist, max_dist)
    end
end

Base.:(<)(p::Cell, q::Cell) = p.max_dist < q.max_dist
Base.:(==)(p::Cell, q::Cell) = p.max_dist == q.max_dist
function Base.hash(cell::Cell, h::UInt)
    h = hash(cell.max_dist, h)
    return hash(Cell, h)
end
 
getx(c::Cell) = c.x

gety(c::Cell) = c.y
