#function add_point_bowyer!(T::Ts, adj, adj2v, DG, pts, r;

function add_point_berg!(T, adj::Adjacent{I,E}, adj2v, DG, HG, pts, r) where {I,E}
    r = I(r)
    Tᵢⱼₖ, flag = locate_triangle(HG, pts, r, (I(LowerRightBoundingIndex), I(UpperBoundingIndex), I(LowerLeftBoundingIndex)))
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

function triangulate_berg(pts;
    IntegerType::Type{I}=Int64,
    EdgeType::Type{E}=NTuple{2,IntegerType},
    TriangleType::Type{V}=NTuple{3,IntegerType},
    EdgesType::Type{Es}=Set{EdgeType},
    TrianglesType::Type{Ts}=Set{TriangleType},
    randomise=true,
    trim=true) where {I,E,V,Es,Ts}
    pt_order = randomise ? shuffle(eachindex(pts)) : eachindex(pts)
    T = Ts()
    adj = Adjacent{I,E}()
    adj2v = Adjacent2Vertex{I,Es,E}()
    DG = DelaunayGraph{I}()
    HG = HistoryGraph{V}()
    add_triangle!(I(-1), I(-2), I(-3), T, adj, adj2v, DG, HG)
    for r in pt_order
        add_point_berg!(T, adj, adj2v, DG, HG, pts, r)
    end
    trim && remove_bounding_triangle!(T, adj, adj2v, DG)
    return T, adj, adj2v, DG, HG
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
