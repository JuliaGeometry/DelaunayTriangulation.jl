const MAX_INITIAL_POLE_CELLS = 10_000

"""
    pole_of_inaccessibility(points, boundary_nodes; precision = nothing)

Finds the pole of inaccessibility for the polygon defined by `points` and `boundary_nodes`. The `boundary_nodes` 
must match the specification given in the documentation. See also [`check_args`](@ref) for this specification.

# Arguments 
- `points`: The points defining the polygon.
- `boundary_nodes`: The boundary nodes defining the polygon.

# Keyword Arguments
- `precision=nothing`: The precision of the returned pole, as a distance: the distance from the returned pole to the boundary is within `precision` of 
   the largest possible distance. If `precision === nothing`, then `max(width, height) / 1000` is used, where `width` and `height` are the dimensions 
   of the bounding box of the polygon, so that the result does not depend on the scale of the polygon.

# Outputs 
- `x`: The x-coordinate of the pole.
- `y`: The y-coordinate of the pole.

# Extended help 
The pole of inaccessibility is the point within a polygon that is furthest from an edge. For DelaunayTriangulation.jl, 
this is useful as it is a representative point for ghost edges that is guaranteed to be inside the polygon, in contrast to for example
a centroid which is not always inside the polygon. Some useful links are [this blog post](https://blog.mapbox.com/a-new-algorithm-for-finding-a-visual-center-of-a-polygon-7c77e6492fbc)
and the [the original repo](https://github.com/mapbox/polylabel). Our implementation is partially based on 
on [the python implementation](https://github.com/Twista/python-polylabel) and [this other Julia implementation](https://github.com/asinghvi17/Polylabel.jl).
"""
function pole_of_inaccessibility(points, boundary_nodes; precision = nothing)
    ## Initiate
    xmin, xmax, ymin, ymax = polygon_bounds(points, boundary_nodes)
    width = xmax - xmin
    height = ymax - ymin
    min_extent = min(width, height)
    if iszero(min_extent)
        if xmin == xmax
            return xmin, midpoint(ymin, ymax)
        else
            return midpoint(xmin, xmax), ymin
        end
    end
    # Covering the bounding box with squares of side min_extent needs about max_extent / min_extent squares, which can exhaust
    # memory for very thin polygons (e.g. the convex hull of nearly collinear points), so we use larger squares in that case.
    max_extent = max(width, height)
    cell_size = max(min_extent, max_extent / MAX_INITIAL_POLE_CELLS)
    half_width = cell_size / 2
    # Cells are refined until their size is about tol, so a tol relative to max_extent also limits the number of cells 
    # needed to cover a ridge of near-optimal points, e.g. along the middle of a long rectangle.
    tol = isnothing(precision) ? max_extent / 1000 : precision

    ## Initialise the priority queue and decide if the polygon centroid of bounding box centroid is the best initial guess
    _, centroid = polygon_features(points, boundary_nodes)
    centroid_cell = Cell(getx(centroid), gety(centroid), zero(half_width), points, boundary_nodes)
    bounding_box_cell = Cell(xmin + width / 2, ymin + height / 2, zero(half_width), points, boundary_nodes)
    best_cell = centroid_cell.dist > bounding_box_cell.dist ? centroid_cell : bounding_box_cell
    queue = CellQueue{number_type(points)}()
    insert_cell!(queue, best_cell)

    ## Now fill the bounding box with more cells
    x = xmin
    while x < xmax
        y = ymin
        while y < ymax
            cell = Cell(x + half_width, y + half_width, half_width, points, boundary_nodes)
            y += cell_size
            insert_cell!(queue, cell)
        end
        x += cell_size
    end

    ## Now let us process all the current cells. The idea is to try and find the best cell, and any bad cells are split into four.
    while !isempty(queue)
        best_cell = process_cell!(queue, best_cell, points, boundary_nodes, tol)
    end

    ## We are done, and the last best_cell is the solution 
    return best_cell.x, best_cell.y
end
function process_cell!(
        queue::CellQueue{F}, best_cell::Cell{F}, points, boundary_nodes,
        precision,
    ) where {F}
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
    cell_1 = Cell(x - h, y - h, h, points, boundary_nodes)
    cell_2 = Cell(x + h, y - h, h, points, boundary_nodes)
    cell_3 = Cell(x - h, y + h, h, points, boundary_nodes)
    cell_4 = Cell(x + h, y + h, h, points, boundary_nodes)
    insert_cell!(queue, cell_1)
    insert_cell!(queue, cell_2)
    insert_cell!(queue, cell_3)
    insert_cell!(queue, cell_4)
    return best_cell
end
