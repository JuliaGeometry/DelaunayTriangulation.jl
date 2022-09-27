using Test
_isless(p::Point{2,T}, q::Point{2,T}) where {T} = isless(p.coords, q.coords)

p0 = Point(5, 5);
p1 = Point(4.5, 2.5);
p2 = Point(2.5, 1.5);
p3 = Point(3, 3.5);
p4 = Point(0, 2);
p5 = Point(1, 5);
p6 = Point(1, 3);
p7 = Point(4, -1);
p8 = Point(-1, 4);
pts = PointSet([p6, p1, p2, p3, p4, p5, p0, p7, p8])
partialsort!(pts.items, 1, lt=_isless, rev=true)
@test pts.items[begin] == p0

function Triangulate(P::PointSet)
    pts = P.items
    partialsort!(pts, 1, lt=_isless, rev=true) # p0 = pts[begin]
    𝒯 = Connectivity{Triangle,3}[connect((1, 0, -1), Triangle)]
    𝒟 = TriDAG()
    add!(𝒟, 𝒯[begin])
    @views shuffle!(pts[begin+1:end])
    for r = (firstindex(pts)+1):lastindex(pts)
        pᵣ = pts[r]
        𝒯ᵢⱼₖ, interior_flag = locate_triangle(pᵣ, 𝒯, 𝒟)
        pᵢ, pⱼ, pₖ = extract_vertices(pts, 𝒯ᵢⱼₖ)
        if interior_flag
            add_edges!(𝒯, 𝒟, 𝒯ᵢⱼₖ, pᵣ)
            legalise_edge!(𝒯, 𝒟, 𝒯ᵢⱼₖ, pᵣ, pᵢ, pⱼ)
            legalise_edge!(𝒯, 𝒟, 𝒯ᵢⱼₖ, pᵣ, pⱼ, pₖ)
            legalise_edge!(𝒯, 𝒟, 𝒯ᵢⱼₖ, pᵣ, pₖ, pᵢ)
        else
            eᵢⱼ, pₖ = locate_edge(pᵣ, 𝒯ᵢⱼₖ)
            𝒯ᵢₗⱼ = adjacent_triangles(𝒯, 𝒯ᵢⱼₖ, eᵢⱼ)
            pₗ = select_adjacent_vertex(𝒯, eᵢⱼ, 𝒯ᵢₗⱼ)
            add_edges!(𝒯, 𝒟, pᵣ, pₖ)
            add_edges!(𝒯, 𝒟, pᵣ, pₗ)
            legalise_edge!(𝒯, 𝒟, 𝒯ᵢₗⱼ, pᵣ, pᵢ, pₗ)
            legalise_edge!(𝒯, 𝒟, 𝒯ᵢₗⱼ, pᵣ, pₗ, pⱼ)
            legalise_edge!(𝒯, 𝒟, 𝒯ᵢⱼₖ, pᵣ, pⱼ, pₖ)
            legalise_edge!(𝒯, 𝒟, 𝒯ᵢⱼₖ, pᵣ, pₖ, pᵢ)
        end
    end
end




function triangulate(pts; shuffle_pts=true)
    # partialsort!(pts, 1, lt=_isless, rev=true) # p0 = pts[begin]
    𝒯 = TriangleType[(LargeRightIdx + 1, LargeLeftIdx, LargeRightIdx)]
    𝒟 = TriDAG()
    𝒜 = TriAdjacent()
    root = 𝒯[begin]
    add!(𝒟, root)
    𝒜[(LargeRightIdx + 1, LargeLeftIdx)] = LargeRightIdx
    𝒜[(LargeLeftIdx, LargeRightIdx)] = LargeRightIdx + 1
    𝒜[(LargeRightIdx, LargeRightIdx + 1)] = LargeLeftIdx
    𝒜[(LargeLeftIdx, LargeRightIdx + 1)] = EmptyIdx
    𝒜[(LargeRightIdx, LargeLeftIdx)] = EmptyIdx
    𝒜[(LargeRightIdx + 1, LargeRightIdx)] = EmptyIdx
    shuffle_pts && @views shuffle!(pts[begin+1:end])
    for r = (firstindex(pts)+1):lastindex(pts)
        pᵣ = pts[r]
        𝒯ᵢⱼₖ, interior_flag = locate_triangle(𝒟, pts, pᵣ, root)
        i, j, k = 𝒯ᵢⱼₖ
        if interior_flag == 1
            @show 𝒯
            add_point!(𝒯, 𝒟, 𝒜, 𝒯ᵢⱼₖ, r)
            legalise_edge!(𝒯, 𝒟, 𝒜, i, j, r, pts)
            legalise_edge!(𝒯, 𝒟, 𝒜, j, k, r, pts)
            legalise_edge!(𝒯, 𝒟, 𝒜, k, i, r, pts)
        else
            eᵢⱼ, pₖ = locate_edge(pᵣ, 𝒯ᵢⱼₖ)
            𝒯ᵢₗⱼ = adjacent_triangles(𝒯, 𝒯ᵢⱼₖ, eᵢⱼ)
            pₗ = select_adjacent_vertex(𝒯, eᵢⱼ, 𝒯ᵢₗⱼ)
            add_edges!(𝒯, 𝒟, pᵣ, pₖ)
            add_edges!(𝒯, 𝒟, pᵣ, pₗ)
            legalise_edge!(𝒯, 𝒟, 𝒯ᵢₗⱼ, pᵣ, pᵢ, pₗ)
            legalise_edge!(𝒯, 𝒟, 𝒯ᵢₗⱼ, pᵣ, pₗ, pⱼ)
            legalise_edge!(𝒯, 𝒟, 𝒯ᵢⱼₖ, pᵣ, pⱼ, pₖ)
            legalise_edge!(𝒯, 𝒟, 𝒯ᵢⱼₖ, pᵣ, pₖ, pᵢ)
        end
    end
    return 𝒯, 𝒟, 𝒜
end


function triangulate(pts; shuffle_pts=true, trim=true)
    # partialsort!(pts, 1, lt=_isless, rev=true) # p0 = pts[begin]
    𝒯, 𝒟, 𝒜, 𝒱𝒩, root = initialise_triangulation()
    shuffle_pts && @views shuffle!(pts[begin+1:end])
    for r = (firstindex(pts)+1):lastindex(pts)
        pᵣ = pts[r]
        𝒯ᵢⱼₖ, interior_flag = locate_triangle(𝒟, pts, pᵣ, root)
        i, j, k = 𝒯ᵢⱼₖ
        if interior_flag == 1
            add_point!(𝒯, 𝒟, 𝒜, 𝒱𝒩, 𝒯ᵢⱼₖ, r)
            legalise_edge!(𝒯, 𝒟, 𝒜, 𝒱𝒩, i, j, r, pts)
            legalise_edge!(𝒯, 𝒟, 𝒜, 𝒱𝒩, j, k, r, pts)
            legalise_edge!(𝒯, 𝒟, 𝒜, 𝒱𝒩, k, i, r, pts)
        else
            eᵢⱼ, pₖ = locate_edge(pᵣ, 𝒯ᵢⱼₖ)
            𝒯ᵢₗⱼ = adjacent_triangles(𝒯, 𝒯ᵢⱼₖ, eᵢⱼ)
            pₗ = select_adjacent_vertex(𝒯, eᵢⱼ, 𝒯ᵢₗⱼ)
            add_edges!(𝒯, 𝒟, pᵣ, pₖ)
            add_edges!(𝒯, 𝒟, pᵣ, pₗ)
            legalise_edge!(𝒯, 𝒟, 𝒯ᵢₗⱼ, pᵣ, pᵢ, pₗ)
            legalise_edge!(𝒯, 𝒟, 𝒯ᵢₗⱼ, pᵣ, pₗ, pⱼ)
            legalise_edge!(𝒯, 𝒟, 𝒯ᵢⱼₖ, pᵣ, pⱼ, pₖ)
            legalise_edge!(𝒯, 𝒟, 𝒯ᵢⱼₖ, pᵣ, pₖ, pᵢ)
        end
    end
    trim && remove_bounding_triangle!(𝒯, 𝒜, 𝒱𝒩)
    return 𝒯, 𝒟, 𝒜, 𝒱𝒩
end
