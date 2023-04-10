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

"""
    CellQueue{T}

A struct representing the priority queue of [`Cell`](@ref)s, using for sorting the cells in a grid
according to their maximum distance.

# Fields
- `queue::PriorityQueue{Cell{T},T,typeof(Base.Order.Reverse)}`

The priority queue of cells.

# Constructors
    CellQueue{T}()

Constructs a new `CellQueue` with elements of type `Cell{T}`.
"""
struct CellQueue{T} # Could a heap be used for this? Duplicate keys could show up...
    queue::PriorityQueue{Cell{T},T,typeof(Base.Order.Reverse)}
    function CellQueue{T}() where {T}
        return new{T}(PriorityQueue{Cell{T},T,typeof(Base.Order.Reverse)}(Base.Order.Reverse))
    end
end
@inline function insert_cell!(queue::CellQueue, cell::Cell)
    return cell ∉ keys(queue.queue) && enqueue!(queue.queue, cell, cell.max_dist)
end
@inline get_next_cell!(queue::CellQueue) = dequeue!(queue.queue)
@inline Base.isempty(queue::CellQueue) = Base.isempty(queue.queue)

"""
    polygon_features(pts, boundary_nodes)

Returns features of the polygon represented by the points `pts` with `boundary_nodes` defining the polygon 
connections. The features returned are `(a, c)`, where `a` is the area of the polygon and 
`c = (cx, cy)` is the centroid. 

!!! notes 

    - The polygon is assumed to be simple, i.e. no self-intersections.
    - The function works with holes, provided `boundary_nodes` represents these as described in the documentation.
    - The polygon is assumed to have a consistent orientation for each boundary. If the orientation is positive, `a > 0`, and `a < 0` otherwise.
"""
function polygon_features(pts, boundary_nodes)
    if has_multiple_curves(boundary_nodes)
        return polygon_features_multiple_curves(pts, boundary_nodes)
    elseif has_multiple_segments(boundary_nodes)
        return polygon_features_multiple_segments(pts, boundary_nodes)
    else
        return polygon_features_single_segment(pts, boundary_nodes)
    end
end
function polygon_features_single_segment(pts, boundary_nodes; scale=Val(true))
    F = number_type(pts)
    cx = zero(F)
    cy = zero(F)
    a = zero(F)
    n_edge = num_boundary_edges(boundary_nodes)
    vᵢ = get_boundary_nodes(boundary_nodes, 1)
    pᵢ = get_point(pts, vᵢ)
    xᵢ, yᵢ = getxy(pᵢ)
    for j in 2:(n_edge+1)
        vᵢ₊₁ = get_boundary_nodes(boundary_nodes, j)
        pᵢ₊₁ = get_point(pts, vᵢ₊₁)
        xᵢ₊₁, yᵢ₊₁ = getxy(pᵢ₊₁)
        area_contribution = xᵢ * yᵢ₊₁ - xᵢ₊₁ * yᵢ
        cx += (xᵢ + xᵢ₊₁) * area_contribution
        cy += (yᵢ + yᵢ₊₁) * area_contribution
        a += area_contribution
        vᵢ, pᵢ, xᵢ, yᵢ = vᵢ₊₁, pᵢ₊₁, xᵢ₊₁, yᵢ₊₁
    end
    if is_true(scale)
        return a / 2, (cx / (3a), cy / (3a))
    else
        return a, (cx, cy)
    end
end
function polygon_features_multiple_segments(pts, boundary_nodes)
    F = number_type(pts)
    cx = zero(F)
    cy = zero(F)
    a = zero(F)
    ns = num_segments(boundary_nodes)
    for i in 1:ns
        bn = get_boundary_nodes(boundary_nodes, i)
        sub_a, (sub_cx, sub_cy) = polygon_features_single_segment(pts, bn;
            scale=Val(false))
        a += sub_a
        cx += sub_cx
        cy += sub_cy
    end
    return a / 2, (cx / (3a), cy / (3a))
end
function polygon_features_multiple_curves(pts, boundary_nodes)
    F = number_type(pts)
    cx = zero(F)
    cy = zero(F)
    a = zero(F)
    nc = num_curves(boundary_nodes)
    for i in 1:nc
        bn = get_boundary_nodes(boundary_nodes, i)
        sub_a, (sub_cx, sub_cy) = polygon_features_multiple_segments(pts, bn)
        cx += sub_a * sub_cx
        cy += sub_a * sub_cy
        a += sub_a
    end
    return a, (cx / a, cy / a)
end

"""
    squared_distance_to_segment(x₁, y₁, x₂, y₂, x, y)

Given a line segment `(x₁, y₁) → (x₂, y₂)` and a query point `(x, y)`, returns the 
squared distance from `(x, y)` to the line segment.
"""
function squared_distance_to_segment(x₁, y₁, x₂, y₂, x, y)
    qp₁_x = x - x₁
    qp₁_y = y - y₁
    p₁p₂_x = x₂ - x₁
    p₁p₂_y = y₂ - y₁
    t = (qp₁_x * p₁p₂_x + qp₁_y * p₁p₂_y) / (p₁p₂_x^2 + p₁p₂_y^2)
    ts = min(max(t, zero(t)), one(t)) # https://math.stackexchange.com/a/330329/861404
    intersect_x = x₁ + ts * p₁p₂_x
    intersect_y = y₁ + ts * p₁p₂_y
    dx = x - intersect_x
    dy = y - intersect_y
    return dx^2 + dy^2
end


"""
    distance_to_polygon(q, pts, boundary_nodes)

Given a polygon represented by the points `pts` with `boundary_nodes` defining the polygon 
connections, and a query point `q`, returns the signed distance from `q` to the polygon. If 
`q` is outside of the polygon, then the returned distance is negative, and if it is inside 
then the distance is positive. Works with holes, provided `boundary_nodes` matches the 
specification of a boundary given in the documentation.
"""
function distance_to_polygon(q, pts, boundary_nodes)
    if has_multiple_curves(boundary_nodes)
        return distance_to_polygon_multiple_curves(q, pts, boundary_nodes)
    elseif has_multiple_segments(boundary_nodes)
        return distance_to_polygon_multiple_segments(q, pts, boundary_nodes)
    else
        return distance_to_polygon_single_segment(q, pts, boundary_nodes)
    end
end
function distance_to_polygon_single_segment(q, pts, boundary_nodes, is_in_outer=false; return_sqrt=Val(true))
    x, y = getxy(q)
    F = number_type(pts)
    dist = typemax(F)
    n_edge = num_boundary_edges(boundary_nodes)
    vᵢ = get_boundary_nodes(boundary_nodes, 1)
    pᵢ = get_point(pts, vᵢ)
    xᵢ, yᵢ = getxy(pᵢ)
    for j in 2:(n_edge+1)
        vᵢ₊₁ = get_boundary_nodes(boundary_nodes, j)
        pᵢ₊₁ = get_point(pts, vᵢ₊₁)
        xᵢ₊₁, yᵢ₊₁ = getxy(pᵢ₊₁)
        if (yᵢ₊₁ > y) ≠ (yᵢ > y)
            intersect_x = (xᵢ - xᵢ₊₁) * (y - yᵢ₊₁) / (yᵢ - yᵢ₊₁) + xᵢ₊₁
            if x < intersect_x
                is_in_outer = !is_in_outer
            end
        end
        new_dist = squared_distance_to_segment(xᵢ, yᵢ, xᵢ₊₁, yᵢ₊₁, x, y)
        dist = new_dist < dist ? new_dist : dist
        vᵢ, pᵢ, xᵢ, yᵢ = vᵢ₊₁, pᵢ₊₁, xᵢ₊₁, yᵢ₊₁
    end
    dist = is_true(return_sqrt) ? sqrt(dist) : dist
    return is_in_outer ? dist : -dist
end
function distance_to_polygon_multiple_segments(q, pts, boundary_nodes, is_in_outer=-one(number_type(pts)); return_sqrt=Val(true))
    F = number_type(pts)
    dist = typemax(F)
    ns = num_segments(boundary_nodes)
    for i in 1:ns
        bn = get_boundary_nodes(boundary_nodes, i)
        new_dist = distance_to_polygon_single_segment(q, pts, bn, is_in_outer == one(F); return_sqrt=Val(false))
        is_in_outer = sign(new_dist)
        new_dist = abs(new_dist)
        dist = new_dist < dist ? new_dist : dist
    end
    dist = is_true(return_sqrt) ? sqrt(dist) : dist
    return is_in_outer * dist
end
function distance_to_polygon_multiple_curves(q, pts, boundary_nodes)
    F = number_type(pts)
    is_in_outer = -one(F)
    dist = typemax(F)
    nc = num_curves(boundary_nodes)
    for i in 1:nc
        bn = get_boundary_nodes(boundary_nodes, i)
        new_dist = distance_to_polygon_multiple_segments(q, pts, bn, is_in_outer == one(F);
            return_sqrt=Val(false))
        is_in_outer = sign(new_dist)
        new_dist = abs(new_dist)
        dist = new_dist < dist ? new_dist : dist
    end
    return is_in_outer * sqrt(dist)
end

"""
    polygon_bounds(pts, boundary_nodes)

Given a polygon represented by the points `pts` with `boundary_nodes` defining the polygon 
connections, returns a bounding box of the polygon. The bounding box is returned 
in the order `(xmin, xmax, ymin, ymax)`.
"""
function polygon_bounds(pts, boundary_nodes)
    if has_multiple_curves(boundary_nodes)
        return polygon_bounds_multiple_segments(pts, get_boundary_nodes(boundary_nodes, 1)) # 1 is the outermost boundary
    elseif has_multiple_segments(boundary_nodes)
        return polygon_bounds_multiple_segments(pts, boundary_nodes)
    else
        return polygon_bounds_single_segment(pts, boundary_nodes)
    end
end
function polygon_bounds_single_segment(pts, boundary_nodes)
    F = number_type(pts)
    xmin, xmax, ymin, ymax = typemax(F), typemin(F), typemax(F), typemin(F)
    n_edge = num_boundary_edges(boundary_nodes)
    for i in 1:n_edge
        vᵢ = get_boundary_nodes(boundary_nodes, i)
        pᵢ = get_point(pts, vᵢ)
        xᵢ, yᵢ = getxy(pᵢ)
        xmin = min(xᵢ, xmin)
        xmax = max(xᵢ, xmax)
        ymin = min(yᵢ, ymin)
        ymax = max(yᵢ, ymax)
    end
    return xmin, xmax, ymin, ymax
end
function polygon_bounds_multiple_segments(pts, boundary_nodes)
    F = number_type(pts)
    xmin, xmax, ymin, ymax = typemax(F), typemin(F), typemax(F), typemin(F)
    ns = num_segments(boundary_nodes)
    for i in 1:ns
        bn = get_boundary_nodes(boundary_nodes, i)
        _xmin, _xmax, _ymin, _ymax = polygon_bounds_single_segment(pts, bn)
        xmin = xmin > _xmin ? _xmin : xmin
        xmax = xmax < _xmax ? _xmax : xmax
        ymin = ymin > _ymin ? _ymin : ymin
        ymax = ymax < _ymax ? _ymax : ymax
    end
    return xmin, xmax, ymin, ymax
end

"""
    pole_of_inaccessibility(pts, boundary_nodes; precision = one(number_type(pts)))

Given a collection of points `pts` and a set of `boundary_nodes` defining the polygon
connections, finds the pole of inaccessibility. This works for multiply-connected polygons, 
provided `boundary_nodes` matches the specification given in the documentation. You can 
control the tolerance of the returned pole using `precision`.

This function is also commonly called `polylabel`.

!!! notes 

    The pole of inaccessibility is a point within a polygon that is furthest from an 
    edge. It is useful for our purposes since it is a representative point that is 
    guaranteed to be inside the polygon, in contrast to for example a centroid which 
    is not always inside the polygon.

    For more information about this, see e.g. [this blog post](https://blog.mapbox.com/a-new-algorithm-for-finding-a-visual-center-of-a-polygon-7c77e6492fbc)
    or [the original repo](https://github.com/mapbox/polylabel). This implementation was partially based 
    on [the python implementation](https://github.com/Twista/python-polylabel) and [this other Julia implementation](https://github.com/asinghvi17/Polylabel.jl).
"""
function pole_of_inaccessibility(pts, boundary_nodes; precision=one(number_type(pts)))
    F = number_type(pts)
    ## Initiate
    xmin, xmax, ymin, ymax = polygon_bounds(pts, boundary_nodes)
    width = xmax - xmin
    height = ymax - ymin
    min_extent = min(width, height)
    half_width = min_extent / 2

    ## Initialise the priority queue and decide if the polygon centroid of bounding box centroid is the best initial guess
    _, centroid = polygon_features(pts, boundary_nodes)
    centroid_cell = Cell(getx(centroid), gety(centroid), zero(half_width), pts, boundary_nodes)
    bounding_box_cell = Cell(xmin + width / 2, ymin + height / 2, zero(half_width), pts, boundary_nodes)
    best_cell = centroid_cell.dist > bounding_box_cell.dist ? centroid_cell : bounding_box_cell
    queue = CellQueue{F}()
    insert_cell!(queue, best_cell)

    ## Now fill the bounding box with more cells
    x = xmin
    while x < xmax
        y = ymin
        while y < ymax
            cell = Cell(x + half_width, y + half_width, half_width, pts, boundary_nodes)
            y += min_extent
            insert_cell!(queue, cell)
        end
        x += min_extent
    end

    ## Now let us process all the current cells. The idea is to try and find the best cell, and any bad cells are split into four.
    while !isempty(queue)
        best_cell = process_cell!(queue, best_cell, pts, boundary_nodes, precision)
    end

    ## We are done, and the last best_cell is the solution 
    return best_cell.x, best_cell.y
end
function process_cell!(queue::CellQueue{F}, best_cell::Cell{F}, pts, boundary_nodes,
    precision) where {F}
    next_cell = get_next_cell!(queue)
    if next_cell.dist > best_cell.dist
        best_cell = next_cell # This cell will have a large circle, so let's choose it 
    end
    if next_cell.max_dist - best_cell.dist ≤ precision
        return best_cell # In the current set of cells, we have now found the largest radius, so we can continue 
    end
    # We now break the cell into four, since we have not yet found the largest radius. This is the "quadtree" part of the approach.
    h = next_cell.half_width / 2
    x = getx(next_cell)
    y = gety(next_cell)
    cell_1 = Cell(x - h, y - h, h, pts, boundary_nodes)
    cell_2 = Cell(x + h, y - h, h, pts, boundary_nodes)
    cell_3 = Cell(x - h, y + h, h, pts, boundary_nodes)
    cell_4 = Cell(x + h, y + h, h, pts, boundary_nodes)
    insert_cell!(queue, cell_1)
    insert_cell!(queue, cell_2)
    insert_cell!(queue, cell_3)
    insert_cell!(queue, cell_4)
    return best_cell
end
