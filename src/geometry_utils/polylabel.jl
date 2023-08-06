"""
    pole_of_inaccessibility(pts, boundary_nodes; precision = one(number_type(pts)))

Given a collection of points `pts` and a set of `boundary_nodes` defining the polygon
connections, finds the pole of inaccessibility. This works for multiply-connected polygons, 
provided `boundary_nodes` matches the specification given in the documentation. You can 
control the tolerance of the returned pole using `precision`.

This function is also commonly called `polylabel`.

!!! note 

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
    if iszero(min_extent)
        if xmin == xmax 
            return xmin, (ymin + ymax)/2 
        else 
            return (xmin + xmax)/2, ymin
        end
    end
    half_width = min_extent / 2

    ## Initialise the priority queue and decide if the polygon centroid of bounding box centroid is the best initial guess
    _, centroid = polygon_features(pts, boundary_nodes)
    centroid_cell = Cell(_getx(centroid), _gety(centroid), zero(half_width), pts, boundary_nodes)
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
    x = _getx(next_cell)
    y = _gety(next_cell)
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