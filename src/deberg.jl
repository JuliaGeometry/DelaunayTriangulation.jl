function triangulate_berg(pts;
    IntegerType::Type{I}=Int64,
    EdgeType::Type{E}=NTuple{2,IntegerType},
    TriangleType::Type{V}=NTuple{3,IntegerType},
    EdgesType::Type{Es}=Set{EdgeType},
    TrianglesType::Type{Ts}=Set{TriangleType}) where {I,E,V,Es,Ts}
    pt_order = shuffle(eachindex(pts))
    T = Ts()
    adj = Adjacent{I,E}()
    adj2v = Adjacent2Vertex{I,Es,E}()
    DG = DelaunayGraph{I}()
    HG = HistoryGraph{V}()
    add_triangle!(I(-1), I(-2), I(-3), T, adj, adj2v, DG, HG)
    for r in pt_order
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
    end
    remove_bounding_triangle!(T, adj, adj2v, DG)
    return T, adj, adj2v, DG, HG
end

function remove_bounding_triangle!(T, adj::Adjacent{I,E}, adj2v, DG) where {I,E}
    for w ∈ (I(LowerRightBoundingIndex), I(UpperBoundingIndex), I(LowerLeftBoundingIndex))
        for (u, v) in get_edge(adj2v, w)
            delete_triangle!(u, v, w, T, adj, adj2v, DG)
        end
        delete_point!(DG, w)
        delete_point!(adj2v, w)
    end
    return nothing
end