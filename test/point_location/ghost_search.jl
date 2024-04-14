using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using LinearAlgebra
using Random
using StatsBase

include("../helper_functions.jl")

tri, label_map, index_map = simple_geometry()
add_ghost_triangles!(tri)
DT.compute_representative_points!(tri)
rep = DT.get_representative_point_list(tri)
rep[1].x = 10.0
rep[1].y = 10.0
_pts = tri.points[[12, 11, 10, 9]]
rep[2].x = mean([8.0, 8.0, 4.0, 4.0])
rep[2].y = mean([16.0, 6.0, 6.0, 16.0])
_pts = tri.points[[18, 17, 16, 15, 14, 13]]
rep[3].x = mean([18.0, 18.0, 14.0, 12.0, 14.0, 14.0])
rep[3].y = mean([12.0, 6.0, 2.0, 4.0, 6.0, 10.0])
pts = get_points(tri)
adj = get_adjacent(tri)
adj2v = get_adjacent2vertex(tri)
ghost_vertex_ranges = get_ghost_vertex_ranges(tri)
ghost_vertex_map = get_ghost_vertex_map(tri)
graph = get_graph(tri)
boundary_nodes = get_boundary_nodes(tri)
k_list = unique(boundary_nodes[1][1])

@testset "Finding points outside" begin
    q = [-3.94232, 15.6067]
    for k in k_list
        i, j = DT.exterior_jump_and_march(tri, k, q)
        @test DT.is_ghost_vertex(get_adjacent(adj, i, j))
        @test i == index_map["h"] && j == index_map["g"]
    end
    q = [-5.0, 5.0]
    for k in k_list
        i, j = DT.exterior_jump_and_march(tri, k, q)
        @test DT.is_ghost_vertex(get_adjacent(adj, i, j))
        @test DT.is_ghost_vertex(get_adjacent(adj, i, j))
        @test i == index_map["a"] && j == index_map["h"]
    end
    q = [-4.8934362297993, -5.6330436693732]
    for k in k_list
        i, j = DT.exterior_jump_and_march(tri, k, q)
        @test DT.is_ghost_vertex(get_adjacent(adj, i, j))
        @test DT.is_ghost_vertex(get_adjacent(adj, i, j))
        @test i == index_map["b"] && j == index_map["a"]
    end
    q = [15.0, -5.0]
    for k in k_list
        i, j = DT.exterior_jump_and_march(tri, k, q)
        @test DT.is_ghost_vertex(get_adjacent(adj, i, j))
        @test DT.is_ghost_vertex(get_adjacent(adj, i, j))
        @test i == index_map["c"] && j == index_map["b"]
    end
    q = [25.0, 0.0]
    for k in k_list
        i, j = DT.exterior_jump_and_march(tri, k, q)
        @test DT.is_ghost_vertex(get_adjacent(adj, i, j))
        @test i == index_map["d"] && j == index_map["c"]
    end
    q = [25.0, 15.0]
    for k in k_list
        i, j = DT.exterior_jump_and_march(tri, k, q)
        @inferred DT.exterior_jump_and_march(tri, k, q)
        @test DT.is_ghost_vertex(get_adjacent(adj, i, j))
        @test i == index_map["e"] && j == index_map["d"]
    end
    q = [13.717, 22.3597]
    for k in k_list
        i, j = DT.exterior_jump_and_march(tri, k, q)
        @test DT.is_ghost_vertex(get_adjacent(adj, i, j))
        @test i == index_map["f"] && j == index_map["e"]
    end
    q = [2.3469, 22.61468]
    for k in k_list
        i, j = DT.exterior_jump_and_march(tri, k, q)
        @test DT.is_ghost_vertex(get_adjacent(adj, i, j))
        @test i == index_map["g"] && j == index_map["f"]
    end
    q = [-5.0, 25.0]
    for k in k_list
        i, j = DT.exterior_jump_and_march(tri, k, q)
        @test DT.is_ghost_vertex(get_adjacent(adj, i, j))
        f, g, h = index_map["f"], index_map["g"], index_map["h"]
        @test (j, i) == (f, g) || (j, i) == (g, h)
    end
    q = [-5.0, 1.0]
    for k in k_list
        i, j = DT.exterior_jump_and_march(tri, k, q)
        @test DT.is_ghost_vertex(get_adjacent(adj, i, j))
        g, h, a = index_map["g"], index_map["h"], index_map["a"]
        @test (j, i) == (g, h) || (j, i) == (h, a)
    end
    q = [-5.0, -5.0]
    for k in k_list
        i, j = DT.exterior_jump_and_march(tri, k, q)
        @test DT.is_ghost_vertex(get_adjacent(adj, i, j))
        h, a, b = index_map["h"], index_map["a"], index_map["b"]
        @test (j, i) == (h, a) || (j, i) == (a, b)
    end
    q = [10.0, -5.0]
    for k in k_list
        i, j = DT.exterior_jump_and_march(tri, k, q)
        @test DT.is_ghost_vertex(get_adjacent(adj, i, j))
        a, b, c = index_map["a"], index_map["b"], index_map["c"]
        @test (j, i) == (a, b) || (j, i) == (b, c)
    end
    q = [60.0, 60.0]
    for k in k_list
        i, j = DT.exterior_jump_and_march(tri, k, q)
        @test DT.is_ghost_vertex(get_adjacent(adj, i, j))
        d, e, f = index_map["d"], index_map["e"], index_map["f"]
        @test (j, i) == (d, e) || (j, i) == (e, f)
    end
end

@testset "Testing points that are already in the triangulation" begin
    for k in DT.each_point_index(pts)
        if DT.is_exterior_boundary_node(tri, k)
            i, j = DT.exterior_jump_and_march(tri, k, get_point(tri, k))
            @test k ∈ (i, j) &&
                  DT.is_exterior_ghost_triangle(tri, i, j, get_adjacent(tri, i, j))
            @test DT.is_on(DT.point_position_relative_to_triangle(tri, i, j, k, k))
        end
    end
end