using .DelaunayTriangulation
using CairoMakie
using StaticArrays
using SimpleGraphs
using Test
using Random
using StatsBase
const DT = DelaunayTriangulation
import DataStructures: DefaultDict

function triplot!(ax, T, pts; markersize=11, plot_ghost_edges=false, DG=nothing, recompute_centroid=true, plot_all_pts=false, kwargs...)
    Tmat = zeros(Int64, length(T), 3)
    tri_length = 0
    for T in T
        if !any(==(DT.BoundaryIndex), (geti(T), getj(T), getk(T)))
            tri_length += 1
            Tmat[tri_length, :] .= [geti(T), getj(T), getk(T)]
        end
    end
    Tmat = @views Tmat[1:tri_length, :]
    pmat = zeros(length(pts), 2)
    for i in DT._eachindex(pts)
        p = _get_point(pts, i)
        pmat[i, :] = [getx(p), gety(p)]
    end
    poly!(ax, pmat, Tmat; kwargs...)
    if plot_all_pts || isnothing(DG)
        scatter!(ax, pmat, color=:red, markersize=markersize)
    else
        scatter!(ax, pmat[setdiff(collect(DG.graph.V), 0), :], color=:red, markersize=markersize)
    end
    if plot_ghost_edges
        if recompute_centroid
            DT.compute_centroid!(pts)
        end
        pcx, pcy = DT.get_point(pts, DT.BoundaryIndex)
        for u in DT.get_neighbour(DG, DT.BoundaryIndex)
            px, py = DT.get_point(pts, u)
            p2x, p2y = (pcx, pcy) .* (1 - 100) .+ 100 .* (px, py)
            lines!(ax, [(px, py), (p2x, p2y)], color=:blue, linestyle=:dash)
        end
    end
end

function DT._get_point(pts::AbstractMatrix, i)
    return @view pts[:, i]
end
function DT._eachindex(pts::AbstractMatrix)
    return axes(pts, 2)
end

function voronoiplot!(ax, vorn, pts, DG; markersize=11, kwargs...)
    makie_polygons = Vector{Makie.Polygon}(undef, num_points(pts))
    xmax = -Inf
    ymax = -Inf
    for i in DT._eachindex(pts)
        p = get_point(pts, i)
        x, y = p
        if x > xmax
            xmax = x
        end
        if y > ymax
            ymax = y
        end
    end
    BN = convex_hull(DG, pts)
    for i in DT._eachindex(pts)
        if DT.BoundaryIndex ∉ vorn.polygons[i]
            makie_polygons[i] = Makie.Polygon(Makie.Point2.(vorn.circumcenters[vorn.polygons[i]]))
        else
            _pts = Vector{Float64}[]
            for (j, v) in pairs(vorn.polygons[i])
                if v == DT.BoundaryIndex # BoundaryIndices come in pairs
                    mid_pt = Float64[]
                    prev_v = vorn.polygons[i][j == 1 ? length(vorn.polygons[i]) : j - 1]
                    next_v = vorn.polygons[i][j == length(vorn.polygons[i]) ? 1 : j + 1]
                    if next_v == DT.BoundaryIndex # then v is the first boundary index
                        τ = vorn.idx_to_triangle[prev_v]
                        a, b, c = indices(τ)
                        if c ∉ DT.get_neighbour(DG, DT.BoundaryIndex)
                            pa, pb = get_point(pts, a, b)
                            pax, pay = pa
                            pbx, pby = pb
                            push!(mid_pt, (pax + pbx) / 2)
                            push!(mid_pt, (pay + pby) / 2) # This is the midpoint of the edge. Now we connect the circumcenter to a point "at infinity"
                        elseif a ∉ DT.get_neighbour(DG, DT.BoundaryIndex)
                            pa, pb = get_point(pts, b, c)
                            pax, pay = pa
                            pbx, pby = pb
                            push!(mid_pt, (pax + pbx) / 2)
                            push!(mid_pt, (pay + pby) / 2) # This is the midpoint of the edge. Now we connect the circumcenter to a point "at infinity"
                        elseif b ∉ DT.get_neighbour(DG, DT.BoundaryIndex)
                            pa, pb = get_point(pts, c, a)
                            pax, pay = pa
                            pbx, pby = pb
                            push!(mid_pt, (pax + pbx) / 2)
                            push!(mid_pt, (pay + pby) / 2) # This is the midpoint of the edge. Now we connect the circumcenter to a point "at infinity"
                        else # all pts are boundary points. This means the point is a corner. So the corner point will have the least degree, so lets put the edge through there 
                            adeg, bdeg, cdeg = DT.num_neighbours(DG, a), DT.num_neighbours(DG, b), DT.num_neighbours(DG, c)
                            if minimum((adeg, bdeg, cdeg)) == adeg
                                pa = get_point(pts, a)
                                pax, pay = pa
                                push!(mid_pt, pax)
                                push!(mid_pt, pay)
                            elseif minimum((adeg, bdeg, cdeg)) == bdeg
                                pb = get_point(pts, b)
                                pbx, pby = pb
                                push!(mid_pt, pbx)
                                push!(mid_pt, pby)
                            else
                                pc = get_point(pts, c)
                                pcx, pcy = pc
                                push!(mid_pt, pcx)
                                push!(mid_pt, pcy)
                            end
                        end
                        prev_circ = vorn.circumcenters[prev_v]
                        if DT.isinconvexhull(pts, BN, prev_circ) == 1
                            t = 10000
                            next_pt = prev_circ .* (1 - t) .+ mid_pt .* t
                            while getx(next_pt)^2 / xmax^2 + gety(next_pt)^2 / ymax^2 > 12
                                t = t / 2
                                next_pt = prev_circ .* (1 - t) .+ mid_pt .* t
                            end
                            push!(_pts, next_pt)
                        else
                            t = 10000
                            next_pt = mid_pt .* (1 - t) .+ prev_circ .* t
                            while getx(next_pt)^2 / xmax^2 + gety(next_pt)^2 / ymax^2 > 12
                                t = t / 2
                                next_pt = mid_pt .* (1 - t) .+ prev_circ .* t
                            end
                            push!(_pts, next_pt)
                        end
                    elseif prev_v == DT.BoundaryIndex
                        τ = vorn.idx_to_triangle[next_v]
                        a, b, c = indices(τ)
                        if c ∉ DT.get_neighbour(DG, DT.BoundaryIndex)
                            pa, pb = get_point(pts, a, b)
                            pax, pay = pa
                            pbx, pby = pb
                            push!(mid_pt, (pax + pbx) / 2)
                            push!(mid_pt, (pay + pby) / 2) # This is the midpoint of the edge. Now we connect the circumcenter to a point "at infinity"
                        elseif a ∉ DT.get_neighbour(DG, DT.BoundaryIndex)
                            pa, pb = get_point(pts, b, c)
                            pax, pay = pa
                            pbx, pby = pb
                            push!(mid_pt, (pax + pbx) / 2)
                            push!(mid_pt, (pay + pby) / 2) # This is the midpoint of the edge. Now we connect the circumcenter to a point "at infinity"
                        elseif b ∉ DT.get_neighbour(DG, DT.BoundaryIndex)
                            pa, pb = get_point(pts, c, a)
                            pax, pay = pa
                            pbx, pby = pb
                            push!(mid_pt, (pax + pbx) / 2)
                            push!(mid_pt, (pay + pby) / 2) # This is the midpoint of the edge. Now we connect the circumcenter to a point "at infinity"
                        else # all pts are boundary points. This means the point is a corner. So the corner point will have the least degree, so lets put the edge through there 
                            adeg, bdeg, cdeg = DT.num_neighbours(DG, a), DT.num_neighbours(DG, b), DT.num_neighbours(DG, c)
                            if minimum((adeg, bdeg, cdeg)) == adeg
                                pa = get_point(pts, a)
                                pax, pay = pa
                                push!(mid_pt, pax)
                                push!(mid_pt, pay)
                            elseif minimum((adeg, bdeg, cdeg)) == bdeg
                                pb = get_point(pts, b)
                                pbx, pby = pb
                                push!(mid_pt, pbx)
                                push!(mid_pt, pby)
                            else
                                pc = get_point(pts, c)
                                pcx, pcy = pc
                                push!(mid_pt, pcx)
                                push!(mid_pt, pcy)
                            end
                        end
                        next_circ = vorn.circumcenters[next_v]
                        if DT.isinconvexhull(pts, BN, next_circ) == 1
                            t = 10000
                            next_pt = next_circ .* (1 - t) .+ mid_pt .* t
                            while getx(next_pt)^2 / xmax^2 + gety(next_pt)^2 / ymax^2 > 12
                                t = t / 2
                                next_pt = next_circ .* (1 - t) .+ mid_pt .* t
                            end
                            push!(_pts, next_pt)
                        else
                            t = 10000
                            next_pt = mid_pt .* (1 - t) .+ next_circ .* t
                            while getx(next_pt)^2 / xmax^2 + gety(next_pt)^2 / ymax^2 > 12
                                t = t / 2
                                next_pt = mid_pt .* (1 - t) .+ next_circ .* t
                            end
                            push!(_pts, next_pt)
                        end
                    end
                else
                    push!(_pts, [vorn.circumcenters[v]...])
                end
            end
            makie_polygons[i] = Makie.Polygon(Makie.Point2.(_pts))
        end
    end
    for i in BN # the boundary orientations can get a bit messed up
        mkpi = makie_polygons[i].exterior.points.parent.data
        bn_idx = collect(eachindex(mkpi))
        DT.sort_boundary!(mkpi, bn_idx, get_point(pts, i))
        fixed_pts = mkpi[bn_idx]
        makie_polygons[i] = Makie.Polygon(fixed_pts)
    end
    poly!(ax, makie_polygons, color=rand(RGBf, length(makie_polygons)); kwargs...)
    scatter!(ax, vorn.circumcenters, color=:purple, markersize=markersize)
    return makie_polygons
end

function example_triangulation()
    p1 = @SVector[0.0, 1.0]
    p2 = @SVector[3.0, -1.0]
    p3 = @SVector[2.0, 0.0]
    p4 = @SVector[-1.0, 2.0]
    p5 = @SVector[4.0, 2.0]
    p6 = @SVector[-2.0, -1.0]
    pts = [p1, p2, p3, p4, p5, p6]
    T = Set{NTuple{3,Int64}}([
        (6, 3, 1),
        (3, 2, 5),
        (4, 1, 5),
        (4, 6, 1),
        (5, 1, 3)
    ])
    A = [
        0 0 1 1 1 1
        0 0 1 0 1 0
        1 1 0 0 1 1
        1 0 0 0 1 1
        1 1 1 1 0 0
        1 0 1 1 0 0
    ]
    DG = DT.DelaunayGraph(UndirectedGraph(A))
    adj = DT.Adjacent(DefaultDict(DT.DefaultAdjacentValue,
        Dict(
            (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
            (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
            (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
            (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
            (5, 1) => 3, (1, 3) => 5, (3, 5) => 1,
            (4, 5) => DT.BoundaryIndex, (5, 2) => DT.BoundaryIndex,
            (2, 3) => DT.BoundaryIndex, (3, 6) => DT.BoundaryIndex,
            (6, 4) => DT.BoundaryIndex
        )
    ))
    adj2v = DT.Adjacent2Vertex(Dict(
        DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 2), (2, 3), (3, 6), (6, 4)]),
        1 => Set{NTuple{2,Int64}}([(5, 4), (3, 5), (6, 3), (4, 6)]),
        2 => Set{NTuple{2,Int64}}([(5, 3)]),
        3 => Set{NTuple{2,Int64}}([(1, 6), (5, 1), (2, 5)]),
        4 => Set{NTuple{2,Int64}}([(1, 5), (6, 1)]),
        5 => Set{NTuple{2,Int64}}([(4, 1), (1, 3), (3, 2)]),
        6 => Set{NTuple{2,Int64}}([(1, 4), (3, 1)])
    ))
    all(all(DT.get_edge(adj, i, j) == k for (i, j) in DT.get_edge(adj2v, k)) for k in keys(adj2v.adjacent2vertex)) &&
        return pts, T, DG, adj, adj2v
end

function example_empty_triangulation()
    p1 = @SVector[0.0, 1.0]
    p2 = @SVector[3.0, -1.0]
    p3 = @SVector[2.0, 0.0]
    pts = [p1, p2, p3]
    T = Set{NTuple{3,Int64}}([])
    A = zeros(Int64, 0, 0)
    DG = DT.DelaunayGraph(UndirectedGraph(A))
    adj = DT.Adjacent(DefaultDict(DT.DefaultAdjacentValue, Dict{NTuple{2,Int64},Int64}()))
    adj2v = DT.Adjacent2Vertex(Dict(DT.BoundaryIndex => Set{NTuple{2,Int64}}()))
    return pts, T, DG, adj, adj2v
end

function circle_three_points(a, b, c)
    ax, ay = a
    bx, by = b
    cx, cy = c
    d1x, d1y = [by - ay, ax - bx]
    d2x, d2y = [cy - ay, ax - cx]
    k = d2x * d1y - d2y * d1x
    s1x, s1y = [(ax + bx) / 2, (ay + by) / 2]
    s2x, s2y = [(ax + cx) / 2, (ay + cy) / 2]
    ℓ = d1x * (s2y - s1y) - d1y * (s2x - s1x)
    m = ℓ / k
    centx, centy = [s2x + m * d2x, s2y + m * d2y]
    dx = centx - ax
    dy = centy - ay
    r = sqrt(dx^2 + dy^2)
    return [centx, centy], r
end

function circle_data(a, b, c)
    θ = LinRange(0, 2π, 500)
    (h, k), r = circle_three_points(a, b, c)
    circx = h .+ r * sin.(θ)
    circy = k .+ r * cos.(θ)
    return circx, circy
end