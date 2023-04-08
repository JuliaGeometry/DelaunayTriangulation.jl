"""
    add_point!(tri::Triangulation, new_point;
    point_indices=get_vertices(tri),
    m=default_num_samples(length(point_indices)),
    try_points=(),
    initial_search_point=integer_type(tri)(select_initial_point(get_points(tri),
                                                                new_point;
                                                                point_indices,
                                                                m,
                                                                try_points)))
    add_point!(tri::Triangulation, new_point_x, new_point_y;
    point_indices=get_vertices(tri),
    m=default_num_samples(length(point_indices)),
    try_points=(),
    rng::AbstractRNG=Random.default_rng(),
    initial_search_point=integer_type(tri)(select_initial_point(get_points(tri),
        (new_point_x, new_point_y),
        point_indices,
        m,
        try_points,
        rng)))                                                                
                                                                
       
Adds the point `new_point` into the triangulation.

This function will not update the convex hull - if you need it to 
be corrected, you could use e.g. [`convex_hull!`](@ref). You should 
also be careful with using this if you have deleted ghost triangles.

# Arguments 
- `tri::Triangulation`: The triangulation. 
- `new_point`: The index of the point in `get_points(tri)` to add. If instead this is simply a point rather than an index, it is assumed it is not present in `get_points(tri)` and instead gets inserted into `get_points(tri)` via `push_point!(tri, new_point)`, meaning `new_point` has index `num_points(tri) + 1`.

# Keyword Arguments 
- `point_indices=get_vertices(tri)`: The currently inserted points in `tri`. 
- `m=default_num_samples(length(point_indices))`: How many points to sample. 
- `try_points()`: Points to consider when sampling for point location. 
- `rng::AbstractRNG=Random.default_rng()`: The RNG.
- `initial_search_point=integer_type(tri)(select_initial_point(get_points(tri), new_point; point_indices, m, try_points)))`: Where to start the point location. 

# Outputs 
There are no outputs, but `tri` is updated in-place.
"""
function add_point!(tri::Triangulation, new_point;
    point_indices=get_vertices(tri),
    m=default_num_samples(length(point_indices)),
    try_points=(),
    rng::AbstractRNG=Random.default_rng(),
    initial_search_point=integer_type(tri)(select_initial_point(get_points(tri),new_point;point_indices,m,try_points,rng)),
    update_representative_point=false)
    if !(new_point isa Integer)
        push_point!(tri, new_point)
        new_point = num_points(tri)
    end
    add_point_bowyer_watson!(tri, new_point, initial_search_point, rng, update_representative_point)
    return nothing
end

function add_point!(tri::Triangulation, new_point_x, new_point_y;
    point_indices=get_vertices(tri),
    m=default_num_samples(length(point_indices)),
    try_points=(),
    rng::AbstractRNG=Random.default_rng(),
    initial_search_point=integer_type(tri)(select_initial_point(get_points(tri),(new_point_x, new_point_y); point_indices, m, try_points, rng)),
    update_representative_point=false)
    push_point!(tri, new_point_x, new_point_y)
    return add_point!(tri, num_points(tri); point_indices=point_indices, m=m, try_points=try_points, rng=rng, initial_search_point=initial_search_point, update_representative_point=update_representative_point)
end