using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using DataStructures
using StableRNGs
const GV = DT.𝒢
using ..DelaunayTriangulation: Certificate

Base.include(@__MODULE__, "../helper_functions.jl")

x, y = complicated_geometry()
rng = StableRNG(99988)
boundary_nodes, points = convert_boundary_points_to_indices(x, y)
tri = triangulate(points; rng, boundary_nodes, delete_ghosts=false)
A = get_area(tri)
refine!(tri; max_area=1e-2A, rng, use_circumcenter=true)

tri2, label_map, index_map = simple_geometry()
add_ghost_triangles!(tri2)
DT.compute_representative_points!(tri2)
pts = get_points(tri2)
adj = get_adjacent(tri2)
adj2v = get_adjacent2vertex(tri2)
boundary_map = get_ghost_vertex_map(tri2)
tri2.representative_point_list[1].x = 10.0
tri2.representative_point_list[1].y = 10.0
graph = DT.get_graph(tri2)
boundary_nodes = get_boundary_nodes(tri2)

@testset "is_ghost_vertex" begin
    @test DT.is_ghost_vertex(GV)
    @test DT.is_ghost_vertex(GV - 1)
    @test DT.is_ghost_vertex(GV - 5)
    @inferred DT.is_ghost_vertex(GV)
end

@testset "is_boundary_edge" begin
    for tri in (tri, tri2)
        for e in each_edge(tri)
            i, j = e
            if haskey(tri.boundary_edge_map, e)
                @test DT.is_boundary_edge(tri, e)
                @test DT.is_boundary_edge(tri, i, j)
            elseif haskey(tri.boundary_edge_map, DT.reverse_edge(e))
                @test DT.is_boundary_edge(tri, DT.reverse_edge(e))
                @test DT.is_boundary_edge(tri, j, i)
            else
                @test !DT.is_boundary_edge(tri, e)
                @test !DT.is_boundary_edge(tri, i, j)
            end
            @inferred DT.is_boundary_edge(tri, e)
            @inferred DT.is_boundary_edge(tri, i, j)
        end
    end
end

@testset "is_boundary_triangle" begin
    for tri in (tri, tri2)
        for T in tri.triangles
            i, j, k = T
            if (i, j) ∈ keys(tri.boundary_edge_map) || (j, k) ∈ keys(tri.boundary_edge_map) || (k, i) ∈ keys(tri.boundary_edge_map)
                @test DT.is_boundary_triangle(tri, T)
                @test DT.is_boundary_triangle(tri, i, j, k)
            else
                @test !DT.is_boundary_triangle(tri, T)
                @test !DT.is_boundary_triangle(tri, i, j, k)
            end
            @inferred DT.is_boundary_triangle(tri, T)
        end
    end
end

@testset "is_ghost_edge" begin
    e = [(1, 2), (3, 5), (GV, 2), (5, GV)]
    results = [false, false, true, true]
    for ((i, j), result) in zip(e, results)
        @test DT.is_ghost_edge(i, j) == result
        @inferred DT.is_ghost_edge(i, j)
        @test DT.is_ghost_edge((i, j)) == result
    end
end

@testset "is_ghost_triangle" begin
    for tri in (tri, tri2)
        for T in tri.triangles
            i, j, k = T
            if any(≤(GV), (i, j, k))
                @test DT.is_ghost_triangle(T)
                @test DT.is_ghost_triangle(i, j, k)
            else
                @test !DT.is_ghost_triangle(T)
                @test !DT.is_ghost_triangle(i, j, k)
            end
            @inferred DT.is_ghost_triangle(T)
            @inferred DT.is_ghost_triangle(i, j, k)
        end
    end
end

@testset "is_interior_curve" begin
    for tri in (tri, tri2)
        for (gv, (i, j)) in tri.ghost_vertex_map
            if i ≠ 1
                @test DT.is_interior_curve(tri, i)
                @test DT.is_interior_ghost_vertex(tri, gv)
            else
                @test !DT.is_interior_curve(tri, i)
                @test !DT.is_interior_ghost_vertex(tri, gv)
            end
            @inferred DT.is_interior_curve(tri, i)
            @inferred DT.is_interior_ghost_vertex(tri, gv)
        end
    end
end

@testset "is_exterior_ghost_triangle" begin
    for tri in (tri, tri2)
        for T in tri.triangles
            i, j, k = T
            u, v, w = DT.sort_triangle(T)
            if DT.is_ghost_vertex(w)
                i, j = DT.map_ghost_vertex(tri, w)
                if i ≠ 1
                    @test !DT.is_exterior_ghost_triangle(tri, T...)
                else
                    @test DT.is_exterior_ghost_triangle(tri, T...)
                end
                @inferred DT.is_exterior_ghost_triangle(tri, T...)
            else
                @test !DT.is_exterior_ghost_triangle(tri, T...)
                @inferred DT.is_exterior_ghost_triangle(tri, T...)
            end
        end
    end
end

@testset "is_exterior_ghost_edge" begin
    for tri in (tri, tri2)
        for e in each_edge(tri)
            i, j = e
            if DT.is_ghost_edge(i, j)
                u = DT.is_ghost_vertex(i) ? i : j
                a, b = DT.map_ghost_vertex(tri, u)
                if a ≠ 1
                    @test !DT.is_exterior_ghost_edge(tri, i, j)
                    @test !DT.is_exterior_ghost_edge(tri, j, i)
                else
                    @test DT.is_exterior_ghost_edge(tri, i, j)
                    @test DT.is_exterior_ghost_edge(tri, j, i)
                end
                @inferred DT.is_exterior_ghost_edge(tri, i, j)
                @inferred DT.is_exterior_ghost_edge(tri, j, i)
            else
                @test !DT.is_exterior_ghost_edge(tri, i, j)
                @test !DT.is_exterior_ghost_edge(tri, j, i)
                @inferred DT.is_exterior_ghost_edge(tri, i, j)
                @inferred DT.is_exterior_ghost_edge(tri, j, i)
            end
        end
    end
end

@testset "is_exterior_boundary_node" begin
    tri_outer_bnd = reduce(vcat, tri.boundary_nodes[1])
    tri2_outer_bnd = tri2.boundary_nodes[1][1]
    for i in each_vertex(tri)
        if i ∈ tri_outer_bnd
            @test DT.is_exterior_boundary_node(tri, i)
        else
            @test !DT.is_exterior_boundary_node(tri, i)
        end
    end
    for i in each_vertex(tri2)
        if i ∈ tri2_outer_bnd
            @test DT.is_exterior_boundary_node(tri2, i)
        else
            @test !DT.is_exterior_boundary_node(tri2, i)
        end
    end
end

@testset "edge_exists" begin
    @test DT.edge_exists(tri, -5, 38)
    @test !DT.edge_exists(tri, -2, 100)
    @test !DT.edge_exists(tri, (-2, 100))
    @test DT.edge_exists(5)
    @test !DT.edge_exists(DT.∅)
end

@testset "has_ghost_triangles" begin
    @test DT.has_ghost_triangles(tri)
    DT.delete_ghost_triangles!(tri)
    @test !DT.has_ghost_triangles(tri)
    DT.add_ghost_triangles!(tri)
    @test DT.has_ghost_triangles(tri)
end

@testset "has_boundary_nodes and is_constrained" begin
    @test DT.has_boundary_nodes(tri)
    @test DT.has_boundary_nodes(tri2)
    _tri = triangulate_rectangle(-3.0, 2.0, 5.0, 17.3, 23, 57; single_boundary=true)
    @test DT.has_boundary_nodes(_tri)
    __tri = triangulate(_tri.points)
    @test !DT.has_boundary_nodes(__tri)
    @test !DT.is_constrained(__tri)
    push!(__tri.all_segments, (1, 2), (2, 3), (4, 5))
    @test !DT.has_boundary_nodes(__tri)
    @test DT.is_constrained(__tri)
end

@testset "Orientation of a bare ghost triangle" begin
    p1 = @SVector[-3.32, 3.53]
    p2 = @SVector[-5.98, 2.17]
    p3 = @SVector[-6.36, -1.55]
    pts = [p1, p2, p3]
    tri = triangulate(pts; delete_ghosts=false)
    @test all(DT.is_positively_oriented(DT.triangle_orientation(tri, T)) for T in each_triangle(tri))
    DT.add_triangle!(get_triangles(tri), (2, 3, -1))
    @test !all(DT.is_positively_oriented(DT.triangle_orientation(tri, T)) for T in each_triangle(tri))
end

@testset "is_boundary_node" begin
    x, y = complicated_geometry()
    rng = StableRNG(99988)
    boundary_nodes, points = convert_boundary_points_to_indices(x, y)
    tri = triangulate(points; rng, boundary_nodes, delete_ghosts=false)
    for (ghost_vertex, segment_index) in get_ghost_vertex_map(tri)
        nodes = get_boundary_nodes(tri, segment_index)
        for node in nodes
            flag1, res1 = DT.is_boundary_node(tri, node)
            @test flag1
            @test res1 ∈ get_ghost_vertex_ranges(tri)[ghost_vertex]
            @test res1 ∈ DT.get_ghost_vertex_range(tri, ghost_vertex)
        end
    end
    reduced_bn = reduce(vcat, reduce(vcat, get_boundary_nodes(tri)))
    for node in each_vertex(tri)
        if node ∉ reduced_bn
            flag, res = DT.is_boundary_node(tri, node)
            @test !flag && res == DT.∅
        end
    end
    tri2, label_map, index_map = simple_geometry()
    for (ghost_vertex, segment_index) in get_ghost_vertex_map(tri2)
        nodes = get_boundary_nodes(tri2, segment_index)
        for node in nodes
            flag1, res1 = DT.is_boundary_node(tri2, node)
            @test flag1
            @test res1 ∈ get_ghost_vertex_ranges(tri2)[ghost_vertex]
            @test res1 ∈ DT.get_ghost_vertex_range(tri2, ghost_vertex)
        end
    end
    reduced_bn = reduce(vcat, reduce(vcat, get_boundary_nodes(tri2)))
    for node in each_vertex(tri)
        if node ∉ reduced_bn
            flag, res = DT.is_boundary_node(tri2, node)
            @test !flag && res == DT.∅
        end
    end
    tri3 = example_with_special_corners()
    ch = get_convex_hull_vertices(tri3)
    for node in ch
        flag2, res2 = DT.is_boundary_node(tri3, node)
        @test flag2
        @test res2 ∈ DT.get_ghost_vertex_ranges(tri3)[DT.𝒢]
        @test res2 ∈ DT.get_ghost_vertex_range(tri3, DT.𝒢)
    end
    for node in each_vertex(tri3)
        if node ∉ ch
            flag2, res2 = DT.is_boundary_node(tri3, node)
            @test !flag2 && res2 == DT.∅
        end
    end
end