"""
    triangulate(points::P;
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
        recompute_representative_point=true
    ) where {P,I,E,V,Es,Ts,M}

Computes the unconstrained Delaunay triangulation of a set of `points`.

# Arguments 
- `points::P`: The set of points to compute the triangulation of. 

# Keyword Arguments 
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

# Outputs 
Returns a [`Triangulation`](@ref).
"""
function triangulate(points::P;
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
                     point_order=get_point_order(points, randomise, skip_points,
                                                 IntegerType, rng),
                     recompute_representative_point=true) where {P,I,E,V,Es,Ts,M}
    points_are_unique(points) || throw("Duplicate points are not allowed.")
    tri = triangulate_bowyer_watson(points;
                                    IntegerType,
                                    EdgeType,
                                    TriangleType,
                                    EdgesType,
                                    TrianglesType,
                                    randomise,
                                    delete_ghosts,
                                    delete_empty_features,
                                    try_last_inserted_point,
                                    skip_points,
                                    num_sample_rule,
                                    rng,
                                    point_order,
                                    recompute_representative_point)
    return tri
end
