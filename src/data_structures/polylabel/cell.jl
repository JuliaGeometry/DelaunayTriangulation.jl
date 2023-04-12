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
    `Cell(x::T, y::T, half_width::T, pts, boundary_nodes)`

Constructs a cell with center `(x, y)` and half-width `half_width`. The cell is assumed to live in the polygon defined by `pts` and `boundary_nodes`.
"""
struct Cell{T}
    x::T
    y::T
    half_width::T
    dist::T
    max_dist::T
    function Cell(x::T, y::T, half_width::T, pts, boundary_nodes) where {T}
        dist = distance_to_polygon((x, y), pts, boundary_nodes)
        max_dist = dist + half_width * sqrt(2)
        return new{T}(x, y, half_width, dist, max_dist)
    end
end
getx(c::Cell) = c.x
gety(c::Cell) = c.y
Base.:(<)(p::Cell, q::Cell) = p.max_dist < q.max_dist
Base.:(==)(p::Cell, q::Cell) = p.max_dist == q.max_dist
function Base.hash(cell::Cell, h::UInt)
    #= 
    If you remove this definition and run the test labelled "A previously broken example" 
    in test/polygon_utils.jl, you get an error where two cells that would typically be == 
    are added into the CellQueue, and so we get a BoundsError since they both map to an index 
    "9", but you end up with only eight keys. This seems to fix it. 
    =#
    h = hash(cell.max_dist, h)
    return hash(Cell, h)
end