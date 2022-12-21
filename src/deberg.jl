function add_point_berg!(T, adj::Adjacent{I,E}, adj2v, DG, HG::HistoryGraph{Tri}, pts, r) where {I,E,Tri}
    r = I(r)
    Tᵢⱼₖ, flag = locate_triangle(HG, pts, r, construct_triangle(Tri, I(LowerRightBoundingIndex), I(UpperBoundingIndex), I(LowerLeftBoundingIndex)))
    if flag == 1
        i, j, k = indices(Tᵢⱼₖ)
        split_triangle!(i, j, k, r, T, adj, adj2v, DG, HG)
        legalise_edge!(i, j, r, T, HG, adj, adj2v, DG, pts)
        legalise_edge!(j, k, r, T, HG, adj, adj2v, DG, pts)
        legalise_edge!(k, i, r, T, HG, adj, adj2v, DG, pts)
    elseif flag == 0
        i, j = find_edge(Tᵢⱼₖ, pts, r)
        k = get_edge(adj, i, j)
        ℓ = get_edge(adj, j, i)
        split_edge!(i, j, r, T, adj, adj2v, DG, HG)
        if !is_boundary_edge(j, i, adj)
            split_edge!(j, i, r, T, adj, adj2v, DG, HG)
        end
        legalise_edge!(i, ℓ, r, T, HG, adj, adj2v, DG, pts)
        legalise_edge!(ℓ, j, r, T, HG, adj, adj2v, DG, pts)
        legalise_edge!(j, k, r, T, HG, adj, adj2v, DG, pts)
        legalise_edge!(k, i, r, T, HG, adj, adj2v, DG, pts)
    end
    return nothing
end

"""
    triangulate_berg(pts;
        IntegerType::Type{I}=Int64,
        EdgeType::Type{E}=NTuple{2,IntegerType},
        TriangleType::Type{V}=NTuple{3,IntegerType},
        EdgesType::Type{Es}=Set{EdgeType},
        TrianglesType::Type{Ts}=Set{TriangleType},
        randomise=true,
        trim=true,
        trim_empty_features=true,
        skip_pts=Set{Int64}()) where {I,E,V,Es,Ts}

Triangulates the set of points `pts` using the de Berg's algorithm.

# Inputs 
- `pts`: The point set.

# Keyword Arguments 
- `IntegerType::Type{I}=Int64`: Type used to represent integers. 
- `EdgeType::Type{E}=NTuple{2,IntegerType}`: Type used to represent edges. 
- `TriangleType::Type{V}=NTuple{3,IntegerType}`: Type used to represent triangles. 
- `EdgesType::Type{Es}=Set{EdgeType}`: Type used to represent a collection of edges.
- `TrianglesType::Type{Ts}=Set{TriangleType}`: Type used to represent a collection of triangles. 
- `randomise=true`: Whether to randomise the insertion order.
- `trim=true`: Whether to remove the ghost triangles at the end.
- `trim_empty_features = true`: Whether to remove keys from the `Adjacent` map that refer to deleted edges, and edges from the `DelaunayGraph` that refer to deleted edges.
- `skip_pts`: There may be points you want to avoid adding into the triangulation, but they still exist in `pts`. If this is the case, add the indices for the points into a `Set` and set this keyword to this set.

# Outputs 
The output is a `Triangulation` struct, containing

- `T`: The set of triangles. 
- `adj`: The [`Adjacent`](@ref) map.
- `adj2v`: The [`Adjacent2Vertex`](@ref) map.
- `DG`: The `[DelaunayGraph`](@ref).
- `pts`: The provided point set.

The second output is `HG`, a [`HistoryGraph`](@ref).
"""
function triangulate_berg(pts;
    IntegerType::Type{I}=Int64,
    EdgeType::Type{E}=NTuple{2,IntegerType},
    TriangleType::Type{V}=NTuple{3,IntegerType},
    EdgesType::Type{Es}=Set{EdgeType},
    TrianglesType::Type{Ts}=Set{TriangleType},
    randomise=true,
    trim=true,
    trim_empty_features=true,
    skip_pts=Set{Int64}()) where {I,E,V,Es,Ts}
    !all_points_are_unique(pts) && throw("The points must all be unique.")
    pt_order = randomise ? shuffle(_eachindex(pts)) : collect(_eachindex(pts))
    setdiff!(pt_order, skip_pts)
    T = construct_triangles(Ts)
    adj = Adjacent{I,E}()
    adj2v = Adjacent2Vertex{I,Es,E}()
    DG = DelaunayGraph{I}()
    HG = HistoryGraph{V}()
    add_triangle!(I(LowerRightBoundingIndex), I(UpperBoundingIndex), I(LowerLeftBoundingIndex), T, adj, adj2v, DG, HG)
    for r in pt_order
        add_point_berg!(T, adj, adj2v, DG, HG, pts, r)
    end
    trim_empty_features && clear_empty_keys!(adj)
    trim_empty_features && clear_empty_points!(DG)
    trim && remove_bounding_triangle!(T, adj, adj2v, DG)
    return Triangulation(T, adj, adj2v, DG, pts), HG
end

function remove_bounding_triangle!(T::Ts, adj::Adjacent{I,E}, adj2v, DG) where {I,E,Ts}
    V = triangle_type(Ts)
    for w ∈ (I(LowerRightBoundingIndex), I(UpperBoundingIndex), I(LowerLeftBoundingIndex))
        for (u, v) in get_edge(adj2v, w) # (u, v, w) is a triangle..
            delete_edge!(adj, w, u)
            delete_edge!(adj, w, v)
            delete_edge!(adj, u, w)
            delete_edge!(adj, v, w)
            delete_edge!(adj2v, u, v, w)
            delete_edge!(adj2v, v, w, u)
            if u ≥ I(FirstPointIndex) && v ≥ I(FirstPointIndex) # This can only be a boundary edge
                add_edge!(adj2v, I(BoundaryIndex), u, v)
                add_edge!(adj, u, v, I(BoundaryIndex))
                add_neighbour!(DG, I(BoundaryIndex), u, v)
            end
            delete_triangle!(T, construct_triangle(V, u, v, w))
        end
        delete_point!(DG, w)
        delete_point!(adj2v, w)
    end
    delete_edge!(adj2v, I(BoundaryIndex), I(LowerRightBoundingIndex), I(LowerLeftBoundingIndex))
    delete_edge!(adj2v, I(BoundaryIndex), I(LowerLeftBoundingIndex), I(UpperBoundingIndex))
    delete_edge!(adj2v, I(BoundaryIndex), I(UpperBoundingIndex), I(LowerRightBoundingIndex))
    return nothing
end
