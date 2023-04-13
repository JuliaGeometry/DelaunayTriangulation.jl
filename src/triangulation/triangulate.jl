"""
    triangulate(points::P; edges=nothing, boundary_nodes=nothing,
        IntegerType::Type{I}=Int64,
        EdgeType::Type{E}=NTuple{2,IntegerType},
        TriangleType::Type{V}=NTuple{3,IntegerType},
        EdgesType::Type{Es}=Set{EdgeType},
        TrianglesType::Type{Ts}=Set{TriangleType},
        randomise=true,
        delete_ghosts=true,
        delete_empty_features=true,
        try_last_inserted_point=true,
        skip_points=Set{IntegerType}(),
        num_sample_rule::M=default_num_samples,
        rng::AbstractRNG=Random.default_rng(),
        point_order=get_point_order(points, randomise, skip_points, IntegerType, rng),
        recompute_representative_point=true,
        delete_holes=true,
        check_arguments=true
    ) where {P,I,E,V,Es,Ts,M}

Computes the unconstrained Delaunay triangulation of a set of `points`. If `edges` is provided, 
they will be inserted. If `boundary_nodes` is provided, a boundary will also be made from these 
nodes, with all triangles inside the boundaries deleted.

# Arguments 
- `points::P`: The set of points to compute the triangulation of. 

# Keyword Arguments 
- `edges=nothing`: Any constrained edges to insert. If `nothing`, an unconstrained triangulation is built. The constrained edges should not intersect each other, and they should not cross over boundary edges.
- `boundary_nodes=nothing`: Any boundaries to define. The specification of these boundary nodes is outlined in the boundary handling section of the documentation. All triangles away from a defined boundary are deleted if `delete_holes`.
- `IntegerType::Type{I}=Int64`: The integer type to use for indexing. 
- `EdgeType::Type{E}=NTuple{2,IntegerType}`: The type to use for representing edges. 
- `TriangleType::Type{V}=NTuple{3,IntegerType}`: The type to use for representing triangles. 
- `EdgesType::Type{Es}=Set{EdgeType}`: The type to use for representing collections of edges. 
- `TrianglesType::Type{Ts}=Set{TriangleType}`: The type to use for representing collections of triangles. 
- `randomise=true`: Whether to randomise the insertion order. 
- `delete_ghosts=true`: Whether to remove the ghost triangles at the end of the triangulation. 
- `delete_empty_features=true`: Whether to delete any empty neighbourhoods and adjacencies at the end of the triangulation. 
- `try_last_inserted_point=true`: When finding the next point, this decides if the previously inserted point should also be attempted. 
- `skip_points=Set{IntegerType}()`: Points to skip over when triangulationg, i.e. points to not include in the triangulation. 
- `num_sample_rule::M=default_num_samples`: A function of the form `n -> Number`, with `n` the number of points currently in the triangulation, that returns the number of points to sample during the point location steps. 
- `rng::AbstractRNG=Random.default_rng()`: The RNG to use.
- `point_order=get_point_order(points, randomise, skip_points, IntegerType, rng)`: The insertion order. 
- `recompute_representative_point=true`: At the end of the triangulation, will recompute the `RepresentativePointList` if `true`.
- `delete_holes=true`: Whether to delete the exterior faces of all boundaries.
- `check_arguments=true`: Whether to check the arguments for validity.
- `check_ghost_triangles=!isnothing(boundary_nodes) && num_curves(boundary_nodes) > 1 && is_true(has_interiors_within_interiors(points, boundary_nodes))`: When you have disjoint domains, the method for deleting holes doesn't always clean up outermost boundaries with index corresponding to BoundaryIndex ($BoundaryIndex). If this argument is set to `true`, then we iterate over every triangle and check that it is valid.

# Outputs 
Returns a [`Triangulation`](@ref).
"""
function triangulate(points::P; edges=nothing, boundary_nodes=nothing,
    IntegerType::Type{I}=Int64,
    EdgeType::Type{E}=!isnothing(edges) ? edge_type(typeof(edges)) : NTuple{2,IntegerType},
    TriangleType::Type{V}=NTuple{3,IntegerType},
    EdgesType::Type{Es}=!isnothing(edges) ? typeof(edges) : Set{EdgeType},
    TrianglesType::Type{Ts}=Set{TriangleType},
    randomise=true,
    delete_ghosts=true,
    delete_empty_features=true,
    try_last_inserted_point=true,
    skip_points=Set{IntegerType}(),
    num_sample_rule::M=default_num_samples,
    rng::AbstractRNG=Random.default_rng(),
    point_order=get_point_order(points, randomise, skip_points,
        IntegerType, rng),
    recompute_representative_point=true,
    delete_holes=true,
    check_arguments=true,
    check_ghost_triangles=!isnothing(boundary_nodes) && num_curves(boundary_nodes) > 1 && is_true(has_interiors_within_interiors(points, boundary_nodes))) where {P,I,E,V,Es,Ts,M}
    check_arguments && check_args(points, boundary_nodes)
    tri = triangulate_bowyer_watson(points;
        IntegerType,
        EdgeType,
        TriangleType,
        EdgesType,
        TrianglesType,
        randomise,
        delete_ghosts=false,
        delete_empty_features=false,
        try_last_inserted_point,
        skip_points,
        num_sample_rule,
        rng,
        point_order,
        recompute_representative_point=false)
    if !isnothing(edges) || !isnothing(boundary_nodes)
        bn_map, bn_range, tri = remake_triangulation_with_constraints(tri, edges, boundary_nodes)
        triangulate_constrained!(tri, bn_map; rng)
        tri = replace_boundary_dict_information(tri, bn_map, bn_range)
    end
    if !isnothing(boundary_nodes) && delete_holes
        delete_holes!(tri, check_ghost_triangles)
        add_boundary_information!(tri)
        !delete_ghosts && add_ghost_triangles!(tri)
    end
    recompute_representative_point && compute_representative_points!(tri)
    delete_ghosts && delete_ghost_triangles!(tri)
    delete_empty_features && clear_empty_features!(tri)
    return tri
end
