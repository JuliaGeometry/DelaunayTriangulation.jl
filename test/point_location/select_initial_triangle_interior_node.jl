using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using LinearAlgebra
using Random
using StableRNGs
using StatsBase



tri, label_map, index_map = simple_geometry()
add_ghost_triangles!(tri)
pts = get_points(tri)
adj = get_adjacent(tri)
adj2v = get_adjacent2vertex(tri)
graph = DT.get_graph(tri)
boundary_nodes = get_boundary_nodes(tri)
ghost_vertex_ranges = get_ghost_vertex_ranges(tri)
boundary_map = get_ghost_vertex_map(tri)
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

@testset "Selecting a random edge" begin
      k = 7
      edges = DT.get_adjacent2vertex(tri, k)
      Random.seed!(2992881)
      rng = StableRNG(2992881)
      i, j = rand(rng, edges)
      pᵢ, pⱼ = pts[i], pts[j]
      rng = StableRNG(2992881)
      u, v, qᵢ, qⱼ = DT.select_random_edge(tri, edges, rng)
      @test u == i && v == j && qᵢ == pᵢ && pⱼ == qⱼ
      @inferred DT.select_random_edge(tri, edges, rng)
end

@testset "Testing for a non-exterior-boundary node" begin
      for _ in 1:350
            local i, j, pᵢ, pⱼ, line_cert_i, line_cert_j, k
            k = index_map["m"] # This point is on the boundary of an interior hole
            # Close
            q = (16.4445425, 19.8427)
            _p, _i, _j, _pᵢ, _pⱼ = DT.select_initial_triangle_interior_vertex(tri, k, q)
            @test _p == get_point(pts, k) &&
                  _i == index_map["v"] &&
                  _j == index_map["z"] &&
                  _pᵢ == get_point(pts, index_map["v"]) &&
                  _pⱼ == get_point(pts, index_map["z"])

            # Further away
            q = (4.008184, 1.077)
            _p, _i, _j, _pᵢ, _pⱼ = DT.select_initial_triangle_interior_vertex(tri, k, q)
            @test _p == get_point(pts, k) &&
                  _j == index_map["w"] &&
                  _i == index_map["u"] &&
                  _pⱼ == get_point(pts, index_map["w"]) &&
                  _pᵢ == get_point(pts, index_map["u"])

            # How are collinearities handled?
            q = (14.0, 0.0)
            _p, _i, _j, _pᵢ, _pⱼ = DT.select_initial_triangle_interior_vertex(tri, k, q)
            @test (_p == get_point(pts, k) &&
                   (_i == DT.𝒢 - 2 || _i == DT.𝒢 - 3) &&
                   _j == index_map["n"] &&
                   _pᵢ == get_point(tri, DT.𝒢 - 2) &&
                   _pⱼ == get_point(pts, index_map["n"])) ||
                  (_p == get_point(pts, k) &&
                   _i == index_map["n"] &&
                   _j == index_map["u"] &&
                   _pᵢ == get_point(tri, index_map["n"]) &&
                   _pⱼ == get_point(pts, index_map["u"]))

            # Now do the same test, but put the point above p
            q = (14.0, 19.0)
            _p, _i, _j, _pᵢ, _pⱼ = DT.select_initial_triangle_interior_vertex(tri, k, q)
            @test (_p == get_point(pts, k) &&
                   _i == index_map["w"] &&
                   _j == index_map["v"] &&
                   _pᵢ == get_point(tri, index_map["w"]) &&
                   _pⱼ == get_point(pts, index_map["v"])) ||
                  (_p == get_point(pts, k) &&
                   _i == index_map["v"] &&
                   _j == index_map["z"] &&
                   _pᵢ == get_point(tri, index_map["v"]) &&
                   _pⱼ == get_point(pts, index_map["z"]))

            # What happens if the point is on an edge? 
            q = (14.0, 7.0)
            _p, _i, _j, _pᵢ, _pⱼ = DT.select_initial_triangle_interior_vertex(tri, k, q)
            @test (_p == get_point(pts, k) &&
                   (_i == DT.𝒢 - 2 || _i == DT.𝒢 - 3) &&
                   _j == index_map["n"] &&
                   _pᵢ == get_point(tri, DT.𝒢 - 2) &&
                   _pⱼ == get_point(pts, index_map["n"])) ||
                  (_p == get_point(pts, k) &&
                   _i == index_map["n"] &&
                   _j == index_map["u"] &&
                   _pᵢ == get_point(tri, index_map["n"]) &&
                   _pⱼ == get_point(pts, index_map["u"]))

            # Now do the same test, but put the point above p 
            q = (14.0, 14.0)
            _p, _i, _j, _pᵢ, _pⱼ = DT.select_initial_triangle_interior_vertex(tri, k, q)
            @test (_p == get_point(pts, k) &&
                   _i == index_map["w"] &&
                   _j == index_map["v"] &&
                   _pᵢ == get_point(tri, index_map["w"]) &&
                   _pⱼ == get_point(pts, index_map["v"])) ||
                  (_p == get_point(pts, k) &&
                   _i == index_map["v"] &&
                   _j == index_map["z"] &&
                   _pᵢ == get_point(tri, index_map["v"]) &&
                   _pⱼ == get_point(pts, index_map["z"]))

            # What if it is an existing triangle?
            q = (12.0, 6.0)
            _p, _i, _j, _pᵢ, _pⱼ = DT.select_initial_triangle_interior_vertex(tri, k, q)
            @test _p == get_point(pts, k) &&
                  _i == index_map["n"] &&
                  _j == index_map["u"] &&
                  _pᵢ == get_point(tri, index_map["n"]) &&
                  _pⱼ == get_point(pts, index_map["u"])

            # Now do the same test, but put the point above p 
            q = (13.288, 15.01)
            _p, _i, _j, _pᵢ, _pⱼ = DT.select_initial_triangle_interior_vertex(tri, k, q)
            @test _p == get_point(pts, k) &&
                  _i == index_map["w"] &&
                  _j == index_map["v"] &&
                  _pᵢ == get_point(tri, index_map["w"]) &&
                  _pⱼ == get_point(pts, index_map["v"])
            q = (16.437, 15.42)
            _p, _i, _j, _pᵢ, _pⱼ = DT.select_initial_triangle_interior_vertex(tri, k, q)
            @test _p == get_point(pts, k) &&
                  _i == index_map["v"] &&
                  _j == index_map["z"] &&
                  _pᵢ == get_point(tri, index_map["v"]) &&
                  _pⱼ == get_point(pts, index_map["z"])
            q = (17.46287, 13.111)
            _p, _i, _j, _pᵢ, _pⱼ = DT.select_initial_triangle_interior_vertex(tri, k, q)
            @test _p == get_point(pts, k) &&
                  _i == index_map["z"] &&
                  _j == index_map["r"] &&
                  _pᵢ == get_point(tri, index_map["z"]) &&
                  _pⱼ == get_point(pts, index_map["r"])

            # Check that everything is fine when the point is outside 
            q = (23.068, 6.92)
            _p, _i, _j, _pᵢ, _pⱼ = DT.select_initial_triangle_interior_vertex(tri, k, q)
            @inferred DT.select_initial_triangle_interior_vertex(tri, k, q)
            @test _p == get_point(pts, k) &&
                  _i == index_map["r"] &&
                  (_j == DT.𝒢 - 2 || _j == DT.𝒢 - 3) &&
                  _pᵢ == get_point(tri, index_map["r"]) &&
                  _pⱼ == get_point(tri, DT.𝒢 - 3)

            # Can interior ghost edges be handled correctly?
            q = (15.5, 9.0)
            _p, _i, _j, _pᵢ, _pⱼ = DT.select_initial_triangle_interior_vertex(tri, k, q)
            @test _p == get_point(pts, k) &&
                  _i == index_map["r"] &&
                  (_j == DT.𝒢 - 2 || _j == DT.𝒢 - 3) &&
                  _pᵢ == get_point(tri, index_map["r"]) &&
                  _pⱼ == get_point(tri, DT.𝒢 - 3)
            q = (14.87, 4.01)
            _p, _i, _j, _pᵢ, _pⱼ = DT.select_initial_triangle_interior_vertex(tri, k, q)
            @test _p == get_point(pts, k) &&
                  _j == index_map["n"] &&
                  (_i == DT.𝒢 - 2 || _i == DT.𝒢 - 3) &&
                  _pⱼ == get_point(tri, index_map["n"]) &&
                  _pᵢ == get_point(tri, DT.𝒢 - 3)
      end
end

@testset "Testing points that are already in the triangulation" begin
      for k in DT.each_point_index(pts)
            local i, j, pᵢ, pⱼ
            if !DT.is_exterior_boundary_node(tri, k)
                  for (i, j) in get_adjacent2vertex(tri, k)
                        p1, i1, j1, pᵢ1, pⱼ1 = DT.select_initial_triangle_interior_vertex(tri, k, get_point(tri, i))
                        p2, i2, j2, pᵢ2, pⱼ2 = DT.select_initial_triangle_interior_vertex(tri, k, get_point(tri, j))
                        if DT.is_ghost_vertex(i)
                              @test i ∈ (i1, j1) || i - 1 ∈ (i1, j1) || i + 1 ∈ (i1, j1)
                        else
                              @test i ∈ (i1, j1)
                        end
                        if DT.is_ghost_vertex(j)
                              @test j ∈ (i2, j2) || j - 1 ∈ (i2, j2) || j + 1 ∈ (i2, j2)
                        else
                              @test j ∈ (i2, j2)
                        end
                  end
            end
            p, i, j, pᵢ, pⱼ = DT.select_initial_triangle_interior_vertex(tri, k, get_point(tri, k))
            @test get_adjacent(tri.adjacent, j, i) == k
      end
end