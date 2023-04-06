using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using DataStructures

const BI = DT.BoundaryIndex

include("../helper_functions.jl")

x, y = complicated_geometry()
tri = generate_mesh(x, y, 2.0; convert_result=true, add_ghost_triangles=true)

tri2, label_map, index_map = simple_geometry()
add_ghost_triangles!(tri2)
DT.compute_representative_points!(tri2)
pts = get_points(tri2)
adj = get_adjacent(tri2)
adj2v = get_adjacent2vertex(tri2)
boundary_map = get_boundary_map(tri2)
tri2.representative_point_list[1].x = 10.0
tri2.representative_point_list[1].y = 10.0
graph = get_graph(tri2)
boundary_nodes = get_boundary_nodes(tri2)

@testset "is_boundary_index" begin
    @test DT.is_boundary_index(BI)
    @test DT.is_boundary_index(BI - 1)
    @test DT.is_boundary_index(BI - 5)
    @inferred DT.is_boundary_index(BI)
end

@testset "is_boundary_edge" begin
    adj = DT.Adjacent(DefaultDict(DT.DefaultAdjacentValue,
        Dict((1, 2) => 5,
            (2, 7) => 1,
            (5, 3) => BI,
            (3, 5) => 7,
            (5, 7) => BI - 1,
            (7, 5) => 13,
            (13, 3) => BI - 2)))
    results = [false, false, true, false, true, false, true, false, false]
    edges = [(1, 2), (2, 7), (5, 3), (3, 5), (5, 7), (7, 5), (13, 3), (23, 5), (36, 3)]
    for (ij, result) in zip(edges, results)
        i, j = ij
        @test DT.is_boundary_edge(ij, adj) ==
              DT.is_boundary_edge(i, j, adj) ==
              result
        @inferred DT.is_boundary_edge(ij, adj)
    end
end

T1 = (1, 5, 3)
T2 = (5, 1, BI)
T3 = (1, 9, 3)
T4 = (10, 11, 13)
T5 = (10, BI - 1, 11)
T6 = (23, 9, 5)
T7 = (BI - 2, 5, 9)
global T = (T1, T2, T3, T4, T5, T6, T7)

@testset "is_boundary_triangle" begin
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    [DT.add_triangle!(adj, T) for T in T]
    results = [true, false, false, true, false, true, false]
    for (T, result) in zip(T, results)
        i, j, k = T
        @test DT.is_boundary_triangle(i, j, k, adj) ==
              DT.is_boundary_triangle(T, adj) ==
              result
        @inferred DT.is_boundary_triangle(i, j, k, adj)
    end
end

@testset "is_ghost_edge" begin
    e = [(1, 2), (3, 5), (BI, 2), (5, BI)]
    results = [false, false, true, true]
    for ((i, j), result) in zip(e, results)
        @test DT.is_ghost_edge(i, j) == result
        @inferred DT.is_ghost_edge(i, j)
        @test DT.is_ghost_edge((i, j)) == result
    end
end

@testset "is_ghost_triangle" begin
    results = [false, true, false, false, true, false, true]
    for (T, result) in zip(T, results)
        i, j, k = T
        @test DT.is_ghost_triangle(i, j, k) ==
              DT.is_ghost_triangle(T) ==
              result
        @inferred DT.is_ghost_triangle(i, j, k)
    end
end

global idx = BI
global map1 = OrderedDict(idx => (1, 1), idx - 1 => (1, 2), idx - 2 => (1, 3), idx - 3 => (1, 4),
    idx - 4 => (2, 1), idx - 5 => (2, 2),
    idx - 6 => (3, 1), idx - 7 => (3, 2))

@testset "is_interior_curve" begin
    @test DT.is_interior_curve(2)
    @test !DT.is_interior_curve(1)
    @test !DT.is_interior_curve(idx, map1)
    @test !DT.is_interior_curve(idx - 1, map1)
    @test !DT.is_interior_curve(idx - 2, map1)
    @test !DT.is_interior_curve(idx - 3, map1)
    @test DT.is_interior_curve(idx - 4, map1)
    @test DT.is_interior_curve(idx - 5, map1)
    @test DT.is_interior_curve(idx - 6, map1)
    @test DT.is_interior_curve(idx - 7, map1)
end

@testset "is_outer_ghost_triangle" begin
    @test DT.is_outer_ghost_triangle(6, 7, idx, map1)
    @test DT.is_outer_ghost_triangle(6, idx, 5, map1)
    @test DT.is_outer_ghost_triangle(idx, 13, 8, map1)
    @test DT.is_outer_ghost_triangle(6, 7, idx - 2, map1)
    @test DT.is_outer_ghost_triangle(6, idx - 2, 5, map1)
    @test DT.is_outer_ghost_triangle(idx - 2, 13, 8, map1)
    @test !DT.is_outer_ghost_triangle(6, 7, idx - 4, map1)
    @test !DT.is_outer_ghost_triangle(6, idx - 4, 5, map1)
    @test !DT.is_outer_ghost_triangle(idx - 4, 13, 8, map1)
    @test !DT.is_outer_ghost_triangle(6, 7, idx - 6, map1)
    @test !DT.is_outer_ghost_triangle(6, idx - 6, 5, map1)
    @test !DT.is_outer_ghost_triangle(idx - 6, 13, 8, map1)
    @test !DT.is_outer_ghost_triangle(6, 18, 91, map1)
end

@testset "is_outer_ghost_edge" begin
    @test DT.is_outer_ghost_edge(6, idx, map1)
    @test DT.is_outer_ghost_edge(idx, 5, map1)
    @test DT.is_outer_ghost_edge(6, idx - 2, map1)
    @test DT.is_outer_ghost_edge(idx - 2, 5, map1)
    @test !DT.is_outer_ghost_edge(6, idx - 4, map1)
    @test !DT.is_outer_ghost_edge(idx - 4, 5, map1)
    @test !DT.is_outer_ghost_edge(6, idx - 6, map1)
    @test !DT.is_outer_ghost_edge(idx - 6, 5, map1)
    @test !DT.is_outer_ghost_edge(6, 4, map1)
    @test !DT.is_outer_ghost_edge(1, 5, map1)
end

@testset "is_outer_boundary_node" begin
    graph = get_graph(tri)
    boundary_ranges = get_boundary_index_ranges(tri)
    @test DT.is_outer_boundary_node(4, graph, boundary_ranges)
    @test DT.is_outer_boundary_node(6, graph, boundary_ranges)
    @test DT.is_outer_boundary_node(9, graph, boundary_ranges)
    @test DT.is_outer_boundary_node(10, graph, boundary_ranges)
    @test DT.is_outer_boundary_node(12, graph, boundary_ranges)
    @test DT.is_outer_boundary_node(1, graph, boundary_ranges)
    @test !DT.is_outer_boundary_node(13, graph, boundary_ranges)
    @test !DT.is_outer_boundary_node(15, graph, boundary_ranges)
    @test !DT.is_outer_boundary_node(69, graph, boundary_ranges)
    @test !DT.is_outer_boundary_node(115, graph, boundary_ranges)
    @test !DT.is_outer_boundary_node(78, graph, boundary_ranges)
    @test !DT.is_outer_boundary_node(326, graph, boundary_ranges)
    @test !DT.is_outer_boundary_node(284, graph, boundary_ranges)
    @test DT.is_outer_boundary_node(tri, 4)
    @test DT.is_outer_boundary_node(tri, 6)
    @test DT.is_outer_boundary_node(tri, 9)
    @test DT.is_outer_boundary_node(tri, 10)
    @test DT.is_outer_boundary_node(tri, 12)
    @test DT.is_outer_boundary_node(tri, 1)
    @test !DT.is_outer_boundary_node(tri, 13)
    @test !DT.is_outer_boundary_node(tri, 15)
    @test !DT.is_outer_boundary_node(tri, 69)
    @test !DT.is_outer_boundary_node(tri, 115)
    @test !DT.is_outer_boundary_node(tri, 78)
    @test !DT.is_outer_boundary_node(tri, 326)
    @test !DT.is_outer_boundary_node(tri, 284)
end

@testset "edge_exists" begin
    @test DT.edge_exists(-5, 38, tri.adjacent)
    @test !DT.edge_exists(-71, 100, tri.adjacent)
    @test !DT.edge_exists((-71, 100), tri.adjacent)
    @test DT.edge_exists((3, 254), tri.adjacent)
    @test DT.edge_exists(tri, -5, 38)
    @test !DT.edge_exists(tri, -2, 100)
    @test !DT.edge_exists(tri, (-2, 100))
    @test DT.edge_exists(tri, (3, 254))
    @test DT.edge_exists(5)
    @test !DT.edge_exists(DT.DefaultAdjacentValue)
end

@testset "has_ghost_triangles" begin
    @test DT.has_ghost_triangles(tri.adjacent, tri.adjacent2vertex)
    @test DT.has_ghost_triangles(tri)
    DT.delete_ghost_triangles!(tri)
    @test !DT.has_ghost_triangles(tri.adjacent, tri.adjacent2vertex)
    @test !DT.has_ghost_triangles(tri)
end

@testset "is_outer_boundary_index" begin
    @test all(i -> DT.is_outer_boundary_index(i, tri.boundary_map),
        [BI, BI - 1, BI - 2, BI - 3])
    @test !any(i -> DT.is_outer_boundary_index(i, tri.boundary_map),
        [1, 2, BI - 4, BI - 5, BI - 6])
    @test all(i -> DT.is_outer_boundary_index(tri, i), [BI, BI - 1, BI - 2, BI - 3])
    @test !any(i -> DT.is_outer_boundary_index(tri, i), [1, 2, BI - 4, BI - 5, BI - 6])
    @test DT.is_outer_boundary_index(BI, tri2.boundary_map)
    @test DT.is_outer_boundary_index(tri2, BI)
    @test !any(i -> DT.is_outer_boundary_index(i, tri2.boundary_map), [1, 2, BI - 2, BI - 3])
    @test !any(i -> DT.is_outer_boundary_index(tri2, i), [1, 2, BI - 1, BI - 2, BI - 3])
end

@testset "has_boundary_nodes and is_constrained" begin
    @test DT.has_boundary_nodes(tri)
    @test DT.has_boundary_nodes(tri2)
    _tri = triangulate_rectangle(-3.0, 2.0, 5.0, 17.3, 23, 57; single_boundary=true)
    @test DT.has_boundary_nodes(_tri)
    __tri = DT.triangulate_bowyer_watson(_tri.points)
    @test !DT.has_boundary_nodes(__tri)
    @test !DT.is_constrained(__tri)
    push!(__tri.constrained_edges, (1, 2), (2, 3), (4, 5))
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
    tri = generate_mesh(x, y, 2.0; convert_result=true, add_ghost_triangles=true)
    for (boundary_index, segment_index) in get_boundary_map(tri)
        nodes = get_boundary_nodes(tri, segment_index)
        for node in nodes
            flag1, res1 = DT.is_boundary_node(node, get_graph(tri), get_boundary_index_ranges(tri))
            flag2, res2 = DT.is_boundary_node(tri, node)
            @test flag1 && flag2 && res1 == res2
            @test res1 ∈ get_boundary_index_ranges(tri)[boundary_index]
        end
    end
    reduced_bn = reduce(vcat, reduce(vcat, get_boundary_nodes(tri)))
    for node in each_vertex(tri)
        if node ∉ reduced_bn
            flag1, res1 = DT.is_boundary_node(node, get_graph(tri), get_boundary_index_ranges(tri))
            flag2, res2 = DT.is_boundary_node(tri, node)
            @test !flag1 && !flag2 && res1 == res2 == DT.DefaultAdjacentValue
        end
    end
    tri2, label_map, index_map = simple_geometry()
    for (boundary_index, segment_index) in get_boundary_map(tri2)
        nodes = get_boundary_nodes(tri2, segment_index)
        for node in nodes
            flag1, res1 = DT.is_boundary_node(node, get_graph(tri2), get_boundary_index_ranges(tri2))
            flag2, res2 = DT.is_boundary_node(tri2, node)
            @test flag1 && flag2 && res1 == res2
            @test res1 ∈ get_boundary_index_ranges(tri2)[boundary_index]
        end
    end
    reduced_bn = reduce(vcat, reduce(vcat, get_boundary_nodes(tri2)))
    for node in each_vertex(tri2)
        if node ∉ reduced_bn
            flag1, res1 = DT.is_boundary_node(node, get_graph(tri2), get_boundary_index_ranges(tri2))
            flag2, res2 = DT.is_boundary_node(tri2, node)
            @test !flag1 && !flag2 && res1 == res2 == DT.DefaultAdjacentValue
        end
    end
    tri3 = example_with_special_corners()
    ch = get_convex_hull_indices(tri3)
    for node in ch
        flag1, res1 = DT.is_boundary_node(node, get_graph(tri3), get_boundary_index_ranges(tri3))
        flag2, res2 = DT.is_boundary_node(tri3, node)
        @test flag1 && flag2 && res1 == res2
        @test res1 ∈ get_boundary_index_ranges(tri3)[DT.BoundaryIndex]
    end
    for node in each_vertex(tri3)
        if node ∉ ch
            flag1, res1 = DT.is_boundary_node(node, get_graph(tri3), get_boundary_index_ranges(tri3))
            flag2, res2 = DT.is_boundary_node(tri3, node)
            @test !flag1 && !flag2 && res1 == res2 == DT.DefaultAdjacentValue
        end
    end
end