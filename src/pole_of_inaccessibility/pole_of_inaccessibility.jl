function pole_of_inaccessibility(points, boundary_nodes; precision = one(number_type(points)))
    # Initiate 
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
    half_width = min_extent / 2

    # Initialise the priority queue and decide if the polygon centroid of bounding box centroid is the best initial guess
    _, centroid = polygon_features(points, boundary_nodes)
    cx, cy = getxy(centroid)
    centroid_cell = Cell(cx, cy, zero(half_width), points, boundary_nodes)
    bounding_box_cell = Cell(xmin + width / 2, ymin + height / 2, zero(half_width), points, boundary_nodes)
    best_cell = centroid_cell.dist > bounding_box_cell.dist ? centroid_cell : bounding_box_cell
    queue = CellQueue{number_type(points)}()
    insert_cell!(queue, best_cell)

    # Now fill the bounding box with more cells
    x = xmin
    while x < xmax
        y = ymin
        while y < ymax
            cell = Cell(x + half_width, y + half_width, half_width, points, boundary_nodes)
            y += min_extent
            insert_cell!(queue, cell)
        end
        x += min_extent
    end

    # Now let us process all the current cells. The idea is to try and find the best cell, and any bad cells are split into four.
    itr = 1
    while !isempty(queue)
        itr += 1
        best_cell = process_cell!(queue, best_cell, points, boundary_nodes, precision)
    end

    # We are done, and the last best_cell is the solution
    return getxy(best_cell)
end

function process_cell!(queue::CellQueue{F}, best_cell::Cell{F}, points, boundary_nodes, precision) where {F}
    next_cell = get_next_cell!(queue)
    if next_cell.dist > best_cell.dist
        best_cell = next_cell # This cell will have a large circle, so let's choose it
    end
    if next_cell.max_dist - best_cell.dist ≤ precision
        return best_cell # In the current set of cells, we have now found the largest radius, so we can continue
    end
    # We now break the cell into four, since we have not yet found the largest radius. This is the "quadtree" part of the approach.
    h = next_cell.half_width / 2    
    x, y = getxy(next_cell)
    corners = ((x-h, y-h), (x+h, y-h), (x-h, y+h), (x+h, y+h))
    foreach(corners) do corner 
        x, y = getxy(corner)
        cell = Cell(x, y, h, points, boundary_nodes)
        insert_cell!(queue, cell)
    end 
    return best_cell
end