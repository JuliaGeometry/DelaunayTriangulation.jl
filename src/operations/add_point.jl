"""
    add_point!(tri::Triangulation, new_point; 
        point_indices=Iterators.filter(!is_boundary_index, get_vertices(tri)),
        m = default_num_samples(length(point_indices)),
        try_points = (),
        initial_search_point = integer_type(tri)(select_initial_point(get_points(tri), new_point; point_indices, m, try_points)))
       
Adds the point `new_point` into the triangulation.

This function will not update the convex hull - if you need it to 
be corrected, you could use e.g. [`convex_hull!`](@ref). You should 
also be careful with using this if you have deleted ghost triangles.

# Arguments 
- `tri::Triangulation`: The triangulation. 
- `new_point`: The index of the point in `get_points(tri)` to add. 

# Keyword Arguments 
- `point_indices=Iterators.filter(!is_boundary_index, get_vertices(tri))`: The currently inserted points in `tri`. 
- `m=default_num_samples(length(point_indices))`: How many points to sample. 
- `try_points()`: Points to consider when sampling for point location. 
- `initial_search_point=integer_type(tri)(select_initial_point(get_points(tri), new_point; point_indices, m, try_points)))`: Where to start the point location. 

# Outputs 
There are no outputs, but `tri` is updated in-place.
"""
function add_point!(tri::Triangulation, new_point;
                    point_indices=get_vertices(tri),
                    m=default_num_samples(length(point_indices)),
                    try_points=(),
                    initial_search_point=integer_type(tri)(select_initial_point(get_points(tri),
                                                                                new_point;
                                                                                point_indices,
                                                                                m,
                                                                                try_points)))
    add_point_bowyer_watson!(tri, new_point, initial_search_point)
    return nothing
end
