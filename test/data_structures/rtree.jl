using ..DelaunayTriangulation
using SpatialIndexing
using CairoMakie
using LinearAlgebra
const DT = DelaunayTriangulation
const SI = SpatialIndexing
_functions.jl")

@testset "Expanding a BoundingBox" begin
    bbox = DT.BoundingBox(2, 5, 13, 20)
    tol = 0.10
    bbox2 = DT.expand(bbox, tol)
    @test bbox2 == DT.BoundingBox(2 - 0.3, 5 + 0.3, 13 - 0.7, 20 + 0.7)
end

@testset "BoundingBox for a set of points" begin
    for _ in 1:100
        points = rand(2, 50000)
        bbox = DT.bounding_box(points)
        a, b, c, d = extrema(points[1, :])..., extrema(points[2, :])...
        _bbox = DT.BoundingBox(a, b, c, d)
        @test bbox == _bbox
    end
end

@testset "BoundingBox for a triangle" begin
    for _ in 1:100
        points = rand(2, 50000)
        i, j, k = rand(1:50000, 3)
        p, q, r = get_point(points, i, j, k)
        px, py = p
        qx, qy = q
        rx, ry = r
        A = min(px, qx, rx)
        B = max(px, qx, rx)
        C = min(py, qy, ry)
        D = max(py, qy, ry)
        bbox = DT.bounding_box(p, q, r)
        _bbox = DT.BoundingBox(A, B, C, D)
        @test bbox == _bbox
    end
end

@testset "DiametralBoundingBox" begin
    for _ in 1:10
        p, q = rand(2), rand(2)
        center = (p .+ q) ./ 2
        radius = norm(p .- q) / 2
        center′, radius′ = DT.diametral_circle(p, q)
        @test center ⪧ center′ && radius ≈ radius

        points = rand(2, 500)
        for i in axes(points, 2)
            for j in axes(points, 2)
                bbox = DT.bounding_box(points, i, j)
                _bbox = DT.get_bounding_box(bbox)
                _edge = DT.get_edge(bbox)
                @test _edge == (i, j)
                p, q = get_point(points, i, j)
                center = (p .+ q) ./ 2
                radius = norm(p .- q) / 2
                a′, b′, c′, d′ = center[1] - radius, center[1] + radius, center[2] - radius, center[2] + radius
                a, b, c, d = _bbox.x.a, _bbox.x.b, _bbox.y.a, _bbox.y.b
                @test (a′, b′, c′, d′) ⪧ (a, b, c, d)
            end
        end
    end
end

@testset "Leaf/Branch ==" begin
    parent = nothing
    bbox = DT.BoundingBox(0.0, 1.0, 2.0, 3.0)
    bbox1 = DT.BoundingBox(2.3, 5.0, 7.7, 13.3)
    bbox2 = DT.BoundingBox(17.7, 23.3, 50.0, 75.0)
    dbx1 = DT.DiametralBoundingBox(bbox1, (1, 2))
    dbx2 = DT.DiametralBoundingBox(bbox2, (5, 10))
    leaf1 = DT.Leaf{DT.Branch}(parent, bbox, [dbx1, dbx2])
    leaf2 = DT.Leaf{DT.Branch}(parent, bbox, [dbx1, dbx2])
    @test leaf1 == leaf2
    leaf2 = DT.Leaf{DT.Branch}(parent, bbox1, [dbx1, dbx2])
    @test leaf1 ≠ leaf2
    leaf2 = DT.Leaf{DT.Branch}(parent, bbox, [dbx1])
    @test leaf1 ≠ leaf2
    leaf1 = DT.Leaf{DT.Branch}(DT.Branch(), bbox, [dbx1, dbx2])
    leaf2 = DT.Leaf{DT.Branch}(DT.Branch(), bbox, [dbx1, dbx2])
    @test leaf1 == leaf2

    leaf1 = DT.Leaf{DT.Branch}(parent, bbox, [dbx1, dbx2])
    leaf2 = DT.Leaf{DT.Branch}(DT.Branch(), bbox, [dbx1, dbx2])
    branch1 = DT.Branch(parent, bbox, [leaf1, leaf2], 1)
    branch2 = DT.Branch(parent, bbox, [leaf1, leaf2], 1)
    @test branch1 == branch2
    branch2 = DT.Branch(parent, bbox, [leaf1, leaf2], 2)
    @test branch1 ≠ branch2
    branch2 = DT.Branch(parent, bbox1, [leaf1, leaf2], 1)
    @test branch1 ≠ branch2
    branch2 = DT.Branch(parent, bbox, [leaf1], 1)
    @test branch1 ≠ branch2
    branch1 = DT.Branch(DT.Branch(), bbox, [leaf1, leaf2], 1)
    branch2 = DT.Branch(DT.Branch(), bbox, [leaf1, leaf2], 1)
    @test branch1 == branch2
end

@testset "Some random RTrees for testing" begin
    @testset "Points" begin
        for _ in 1:50
            n = rand(1:1001)
            points = 5 .* randn(2, n)
            tree_si = SI.RTree{Float64,2}(Int, NTuple{2,Int}; variant=SI.RTreeLinear)
            tree = DT.RTree()
            for idx in axes(points, 2) |> shuffle
                x, y = get_point(points, idx)
                insert!(tree_si, SI.Rect((x, y), (x, y)), idx, (idx, idx))
                @test tree ≠ tree_si
                insert!(tree, DT.DiametralBoundingBox(DT.BoundingBox((x, y)), (idx, idx)))
                @test tree == tree_si

                for _ in 1:20
                    p = Tuple(10 .* randn(2))
                    ints = DT.get_intersections(tree, p; cache_id=rand(1:2))
                    ints_si = SI.intersects_with(tree_si, SI.Rect(SI.Point(p)))
                    @test ints == ints_si
                    @inferred collect(ints)
                end
            end
            for idx in axes(points, 2) |> shuffle
                x, y = get_point(points, idx)
                rect_si = SI.Rect((x, y), (x, y))
                rect = DT.DiametralBoundingBox(DT.BoundingBox((x, y)), (idx, idx))
                delete!(tree, rect)
                @test tree ≠ tree_si
                delete!(tree_si, rect_si)
                @test tree == tree_si

                for _ in 1:20
                    p = Tuple(10 .* rand(2))
                    ints = DT.get_intersections(tree, p)
                    ints_si = SI.intersects_with(tree_si, SI.Rect(SI.Point(p)))
                    @test ints == ints_si
                end
            end
        end
    end

    @testset "Edges" begin
        for _ in 1:25
            n = rand(1:250)
            points = 5 .* randn(2, n)
            tree_si = SI.RTree{Float64,2}(Int, NTuple{2,Int}; variant=SI.RTreeLinear)
            tree = DT.BoundaryRTree(points)
            ctr = 1
            range1 = rand(axes(points, 2), max(isqrt(n), 5))
            range2 = rand(axes(points, 2), max(isqrt(n), 5))
            for idx in range1
                for idx2 in range2
                    i, j = idx, idx2
                    si_rect = si_diametral_bounding_box(points, minmax(i, j)...)
                    insert!(tree_si, si_rect, ctr, minmax(i, j))
                    @test tree.tree ≠ tree_si
                    insert!(tree, i, j)
                    @test tree.tree == tree_si
                    ctr += 1

                    for _ in 1:5
                        a, b, c, d = 10 .* randn(4)
                        a, b = minmax(a, b)
                        c, d = minmax(c, d)
                        rect = DT.BoundingBox(a, b, c, d)
                        rect_si = SI.Rect((a, c), (b, d))
                        int = DT.get_intersections(tree, rect; cache_id=rand(1:2))
                        int_si = SI.intersects_with(tree_si, rect_si)
                        @test int == int_si

                        i = rand(1:n)
                        int = DT.get_intersections(tree, i)
                        p = SI.Point(get_point(points, i))
                        int_si = SI.intersects_with(tree_si, SI.Rect(p))
                        @test int == int_si

                        i, j = rand(1:n, 2)
                        int = DT.get_intersections(tree, i, j)
                        bbox = DT.bounding_box(points, minmax(i, j)...).bounding_box
                        bbox_si = SI.Rect((bbox.x.a, bbox.y.a), (bbox.x.b, bbox.y.b))
                        int_si = SI.intersects_with(tree_si, bbox_si)
                        @test int == int_si

                        i, j, k = rand(1:n, 3)
                        int = DT.get_intersections(tree, i, j, k; cache_id=rand(1:2))
                        bbox = DT.bounding_box(get_point(points, i, j, k)...)
                        bbox_si = SI.Rect((bbox.x.a, bbox.y.a), (bbox.x.b, bbox.y.b))
                        int_si = SI.intersects_with(tree_si, bbox_si)
                        @test int == int_si
                    end
                end
            end
            for idx in range1
                for idx2 in range2
                    i, j = minmax(idx, idx2)
                    si_rect = si_diametral_bounding_box(points, i, j)
                    delete!(tree, i, j)
                    @test tree.tree ≠ tree_si
                    delete!(tree_si, si_rect)
                    @test tree.tree == tree_si

                    for _ in 1:5
                        a, b, c, d = 10 .* randn(4)
                        a, b = minmax(a, b)
                        c, d = minmax(c, d)
                        rect = DT.BoundingBox(a, b, c, d)
                        rect_si = SI.Rect((a, c), (b, d))
                        int = DT.get_intersections(tree, rect)
                        int_si = SI.intersects_with(tree_si, rect_si)
                        @test int == int_si

                        i = rand(1:n)
                        int = DT.get_intersections(tree, i)
                        p = SI.Point(get_point(points, i))
                        int_si = SI.intersects_with(tree_si, SI.Rect(p))
                        @test int == int_si

                        i, j = rand(1:n, 2)
                        int = DT.get_intersections(tree, i, j)
                        bbox = DT.bounding_box(points, minmax(i, j)...).bounding_box
                        bbox_si = SI.Rect((bbox.x.a, bbox.y.a), (bbox.x.b, bbox.y.b))
                        int_si = SI.intersects_with(tree_si, bbox_si)
                        @test int == int_si
                    end
                end
            end
        end
    end
end

@testset "split_edge!" begin
    n = 25
    points = 5 .* randn(2, n)
    tree = DT.BoundaryRTree(points)
    insert!(tree, 1, 2)
    rects, els = get_dt_rectangles(tree.tree)
    @test length(rects) == length(els) == 1
    @test DT.get_edge(els[1]) == (1, 2)
    DT.split_edge!(tree, 1, 2, 10)
    rects, els = get_dt_rectangles(tree.tree)
    @test length(rects) == 1 && length(els) == 2
    @test DT.get_edge(els[1]) == (1, 10) || DT.get_edge(els[2]) == (1, 10)
    @test DT.get_edge(els[1]) == (2, 10) || DT.get_edge(els[2]) == (2, 10)
    delete!(tree, 10, 1)
    rects, els = get_dt_rectangles(tree.tree)
    @test length(rects) == length(els) == 1
    @test DT.get_edge(els[1]) == (2, 10)
end

@testset "RTree ==" begin
    tree1 = DT.RTree()
    tree2 = DT.RTree()
    @test tree1 == tree2
    insert!(tree1, DT.DiametralBoundingBox(DT.BoundingBox((0.0, 0.0)), (1, 2)))
    @test tree1 ≠ tree2
    insert!(tree2, DT.DiametralBoundingBox(DT.BoundingBox((0.0, 0.0)), (1, 2)))
    @test tree1 == tree2
    insert!(tree1, DT.DiametralBoundingBox(DT.BoundingBox((0.3, 0.0)), (1, 3)))
    @test tree1 ≠ tree2
    insert!(tree2, DT.DiametralBoundingBox(DT.BoundingBox((0.3, 0.0)), (1, 3)))
    @test tree1 == tree2
end

@testset "BoundaryRTree ==" begin
    n = 25
    points = 5 .* randn(2, n)
    tree1 = DT.BoundaryRTree(points)
    tree2 = DT.BoundaryRTree(points)
    @test tree1 == tree2
    tree3 = DT.BoundaryRTree(rand(2, 50))
    @test tree1 ≠ tree3
    insert!(tree1, 1, 2)
    @test tree1 ≠ tree2
    insert!(tree2, 1, 2)
    @test tree1 == tree2
end