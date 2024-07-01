"""
    Cell{T}

A cell in a grid. The cell is a square with side length `2half_width`. The cell is centered at `(x, y)`. The cell is 
assumed to live in a polygon.

# Fields 
- `x::T`

The x-coordinate of the center of the cell.
- `y::T`

The y-coordinate of the center of the cell.
- `half_width::T`

The half-width of the cell.
- `dist::T`

The distance from the center of the cell to the polygon.
- `max_dist::T`

The maximum distance from the center of the cell to the polygon. This is `dist + half_width * sqrt(2)`.

# Constructors
    Cell(x::T, y::T, half_width::T, points, boundary_nodes)

Constructs a cell with center `(x, y)` and half-width `half_width`. The cell is assumed to live in the polygon defined by `points` and `boundary_nodes`.
"""
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

"""
    getx(c::Cell) -> Number

Returns the x-coordinate of `c`.
"""
getx(c::Cell) = c.x

"""
    gety(c::Cell) -> Number

Returns the y-coordinate of `c`.
"""
gety(c::Cell) = c.y