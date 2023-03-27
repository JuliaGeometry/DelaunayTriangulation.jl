using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using StableRNGs

include("./helper_functions.jl")

@testset "is_true" begin
      @test DT.is_true(true)
      @test DT.is_true(Val(true))
      @test !DT.is_true(false)
      @test !DT.is_true(Val(false))
      @test DT.is_true(Val{true})
      @test !DT.is_true(Val{false})
end

@testset "number_type" begin
      @test number_type([1, 2, 3]) == Int64
      @test number_type([1.0, 2.0, 3.0]) == Float64
      @test number_type([1.0 2.0; 3.0 3.5; 10.0 17.3]) == Float64
      @test number_type((1.0, 5.0)) == Float64
      @test number_type([(1.0f0, 2.0f0), (1.7f0, 2.5f0)]) == Float32
      @test number_type(2.4) == Float64
end

@testset "get_boundary_index" begin
      @test DT.get_boundary_index(1, 2, -3) == -3
      @test DT.get_boundary_index(1, 2, -1) == -1
      @test DT.get_boundary_index(1, -5, 2) == -5
      @test DT.get_boundary_index(-1, 2, 3) == -1
      @test_throws ArgumentError DT.get_boundary_index(2, 5, 7)
      @test DT.get_boundary_index(1, -2) == -2
      @test DT.get_boundary_index(-5, 1) == -5
      @test_throws ArgumentError DT.get_boundary_index(2, 5)
end

@testset "rotate_ghost_triangle_to_standard_form" begin
      @test DT.rotate_ghost_triangle_to_standard_form(1, 2, DT.BoundaryIndex) ==
            (1, 2, DT.BoundaryIndex)
      @test DT.rotate_ghost_triangle_to_standard_form(DT.BoundaryIndex - 2, 2, 3) ==
            (2, 3, DT.BoundaryIndex - 2)
      @test DT.rotate_ghost_triangle_to_standard_form(5, DT.BoundaryIndex - 1, 3) ==
            (3, 5, DT.BoundaryIndex - 1)
      @test DT.rotate_ghost_triangle_to_standard_form((1, 5, DT.BoundaryIndex)) ==
            (1, 5, DT.BoundaryIndex)
      @test DT.rotate_ghost_triangle_to_standard_form([1, DT.BoundaryIndex - 2, 7]) ==
            [7, 1, DT.BoundaryIndex - 2]
      @test DT.rotate_ghost_triangle_to_standard_form((DT.BoundaryIndex - 10, 5, 3)) ==
            (5, 3, DT.BoundaryIndex - 10)
end

@testset "get_left/right_boundary_node" begin
      x, y = complicated_geometry()
      tri = generate_mesh(x, y, 2.0; convert_result=true, add_ghost_triangles=true)
      nodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
      right = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1]
      left = [12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
      @inferred DT.get_right_boundary_node(tri.adjacent, 1, DT.BoundaryIndex,
            tri.boundary_index_ranges, Val(true))
      @inferred DT.get_left_boundary_node(tri.adjacent, 1, DT.BoundaryIndex,
            tri.boundary_index_ranges, Val(true))
      @inferred DT.get_right_boundary_node(tri, 1, DT.BoundaryIndex)
      @inferred DT.get_left_boundary_node(tri, 1, DT.BoundaryIndex)

      for (i, r, ℓ) in zip(nodes, right, left)
            @test DT.get_right_boundary_node(tri.adjacent, i, DT.BoundaryIndex,
                        tri.boundary_index_ranges, true) ==
                  DT.get_right_boundary_node(tri, i, DT.BoundaryIndex) == r
            @test DT.get_left_boundary_node(tri.adjacent, i, DT.BoundaryIndex,
                        tri.boundary_index_ranges, true) ==
                  DT.get_left_boundary_node(tri, i, DT.BoundaryIndex) == ℓ
            @inferred DT.get_right_boundary_node(tri.adjacent, i, DT.BoundaryIndex,
                  tri.boundary_index_ranges, false)
            @inferred DT.get_left_boundary_node(tri, i, DT.BoundaryIndex)
      end
      @test !DT.edge_exists(DT.get_left_boundary_node(tri.adjacent, 1, DT.BoundaryIndex,
            tri.boundary_index_ranges, Val(false)))
      @test !DT.edge_exists(DT.get_right_boundary_node(tri.adjacent, 12, DT.BoundaryIndex,
            tri.boundary_index_ranges, Val(false)))
      nodes = [13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61]
      right = [nodes[2:end]..., nodes[1]]
      left = [nodes[end]..., nodes[1:(end-1)]...]
      for (i, r, ℓ) in zip(nodes, right, left)
            @test DT.get_right_boundary_node(tri.adjacent, i, DT.BoundaryIndex - 4,
                        tri.boundary_index_ranges, true) ==
                  DT.get_right_boundary_node(tri, i, DT.BoundaryIndex - 4) == r
            @test DT.get_left_boundary_node(tri.adjacent, i, DT.BoundaryIndex - 4,
                        tri.boundary_index_ranges, true) ==
                  DT.get_left_boundary_node(tri, i, DT.BoundaryIndex - 4) == ℓ
      end
      nodes = [62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110]
      right = [nodes[2:end]..., nodes[1]]
      left = [nodes[end]..., nodes[1:(end-1)]...]
      for (i, r, ℓ) in zip(nodes, right, left)
            @test DT.get_right_boundary_node(tri.adjacent, i, DT.BoundaryIndex - 5,
                        tri.boundary_index_ranges, true) ==
                  DT.get_right_boundary_node(tri, i, DT.BoundaryIndex - 5) == r
            @test DT.get_left_boundary_node(tri.adjacent, i, DT.BoundaryIndex - 5,
                        tri.boundary_index_ranges, true) ==
                  DT.get_left_boundary_node(tri, i, DT.BoundaryIndex - 5) == ℓ
      end
      nodes = [111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122]
      right = [112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 111]
      left = [122, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121]
      for (i, r, ℓ) in zip(nodes, right, left)
            @test DT.get_right_boundary_node(tri.adjacent, i, DT.BoundaryIndex - 6,
                        tri.boundary_index_ranges, true) ==
                  DT.get_right_boundary_node(tri, i, DT.BoundaryIndex - 7) == r
            @test DT.get_left_boundary_node(tri.adjacent, i, DT.BoundaryIndex - 7,
                        tri.boundary_index_ranges, true) ==
                  DT.get_left_boundary_node(tri, i, DT.BoundaryIndex - 8) == ℓ
      end
      nodes = [123, 128, 124, 125, 126, 127]
      right = [128, 124, 125, 126, 127, 123]
      left = [127, 123, 128, 124, 125, 126]
      for (i, r, ℓ) in zip(nodes, right, left)
            @test DT.get_right_boundary_node(tri.adjacent, i, DT.BoundaryIndex - 10,
                        tri.boundary_index_ranges, true) ==
                  DT.get_right_boundary_node(tri, i, DT.BoundaryIndex - 10) == r
            @test DT.get_left_boundary_node(tri.adjacent, i, DT.BoundaryIndex - 10,
                        tri.boundary_index_ranges, true) ==
                  DT.get_left_boundary_node(tri, i, DT.BoundaryIndex - 10) == ℓ
      end
end

@testset "find_edge" begin
      points = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (0.5, 0.0), (0.5, 0.5), (0.0, 0.5)]
      T = (1, 2, 3)
      ℓ = 4
      @test DT.find_edge(T, points, ℓ) == (1, 2)
      T = (2, 3, 1)
      @test DT.find_edge(T, points, ℓ) == (1, 2)
      T = (1, 2, 3)
      ℓ = 5
      @test DT.find_edge(T, points, ℓ) == (2, 3)
      T = (2, 3, 1)
      @test DT.find_edge(T, points, ℓ) == (2, 3)
      T = (1, 2, 3)
      ℓ = 6
      @test DT.find_edge(T, points, ℓ) == (3, 1)
      T = (2, 3, 1)
      @test DT.find_edge(T, points, ℓ) == (3, 1)
      p1 = [2.0, 3.5]
      p2 = [0.0, 0.0]
      p3 = [3.0, 0.0]
      p4 = [17.2, -2.5]
      p5 = [0.0, 3.0]
      T = DT.construct_triangle(NTuple{3,Int64}, 2, 3, 5)
      pts = [p1, p2, p3, p4, p5]
      push!(pts, [1.0, 0.0])
      @test DT.find_edge(T, pts, length(pts)) == (2, 3)
      push!(pts, [2.0, 0.0])
      @test DT.find_edge(T, pts, length(pts)) == (2, 3)
      push!(pts, [1.5, 0.0])
      @test DT.find_edge(T, pts, length(pts)) == (2, 3)
      push!(pts, [1.0, 0.0])
      @test DT.find_edge(T, pts, length(pts)) == (2, 3)
      push!(pts, [0.5, 0.0])
      @test DT.find_edge(T, pts, length(pts)) == (2, 3)
      push!(pts, [2.5, 0.5])
      @test DT.find_edge(T, pts, length(pts)) == (3, 5)
      push!(pts, [2.0, 1.0])
      @test DT.find_edge(T, pts, length(pts)) == (3, 5)
      push!(pts, [1.5, 1.5])
      @test DT.find_edge(T, pts, length(pts)) == (3, 5)
      push!(pts, [1.0, 2.0])
      @test DT.find_edge(T, pts, length(pts)) == (3, 5)
      push!(pts, [0.5, 2.5])
      @test DT.find_edge(T, pts, length(pts)) == (3, 5)
      push!(pts, [0.0, 2.5])
      @test DT.find_edge(T, pts, length(pts)) == (5, 2)
      push!(pts, [0.0, 2.2])
      @test DT.find_edge(T, pts, length(pts)) == (5, 2)
      push!(pts, [0.0, 2.0])
      @test DT.find_edge(T, pts, length(pts)) == (5, 2)
      push!(pts, [0.0, 1.5])
      @test DT.find_edge(T, pts, length(pts)) == (5, 2)
      push!(pts, [0.0, 0.8])
      @test DT.find_edge(T, pts, length(pts)) == (5, 2)
      push!(pts, [0.0, 0.2])
      @test DT.find_edge(T, pts, length(pts)) == (5, 2)
end

@testset "choose_uvw" begin
      i, j, k = rand(Int64, 3)
      @test DT.choose_uvw(true, false, false, i, j, k) == (i, j, k)
      @test DT.choose_uvw(false, true, false, i, j, k) == (j, k, i)
      @test DT.choose_uvw(false, false, true, i, j, k) == (k, i, j)
end

@testset "is_circular" begin
      x = rand(10)
      @test !DT.is_circular(x)
      push!(x, x[begin])
      @test DT.is_circular(x)
      @test DT.is_circular([])
end

@testset "circular_equality" begin
      @test !DT.circular_equality([1, 2, 3, 4, 1], [3, 2, 4, 1, 3])
      @test !DT.circular_equality([1, 2, 3, 1], [3, 4, 5, 3])
      @test_throws AssertionError DT.circular_equality([1, 2, 3, 1], [3, 4, 5])
      @test_throws AssertionError DT.circular_equality([1, 2, 3], [5, 4, 3])
      @test !DT.circular_equality([1, 2, 3, 4, 1], [1, 2, 1])
      x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1]
      for i in 1:10
            local y
            y = @views circshift(x[begin:(end-1)], i)
            push!(y, y[begin])
            @test DT.circular_equality(x, y) && DT.circular_equality(y, x)
      end
      @test DT.circular_equality([3, 2, 1, 13, 12, 11, 5, 4, 3], [1, 13, 12, 11, 5, 4, 3, 2, 1])
      @test DT.circular_equality([], [])
end

@testset "get_surrounding_polygon" begin
      tri = example_triangulation()
      rng = StableRNG(999987897899)
      tri = triangulate(get_points(tri); delete_ghosts=false, rng)
      polys = Dict(
            1 => [6, 3, 7, 4, 6],
            2 => [8, 3, 6, DT.BoundaryIndex, 8],
            3 => [8, 7, 1, 6, 2, 8],
            4 => [6, 1, 7, 5, DT.BoundaryIndex, 6],
            5 => [4, 7, 8, DT.BoundaryIndex, 4],
            6 => [2, 3, 1, 4, DT.BoundaryIndex, 2],
            7 => [8, 5, 4, 1, 3, 8],
            8 => [5, 7, 3, 2, DT.BoundaryIndex, 5],
            DT.BoundaryIndex => [5, 8, 2, 6, 4, 5]
      )
      fnc_polys = Dict{Int64,Vector{Int64}}()
      for i in keys(polys)
            fnc_polys[i] = DT.get_surrounding_polygon(tri, i)
            push!(fnc_polys[i], fnc_polys[i][begin])
      end
      for (poly_true, poly_f) in zip(values(polys), values(fnc_polys))
            @test DT.circular_equality(poly_true, poly_f)
      end
      fnc_polys = Dict{Int64,Vector{Int64}}()
      for i in keys(polys)
            fnc_polys[i] = DT.get_surrounding_polygon(tri, i; skip_boundary_indices=true)
            push!(fnc_polys[i], fnc_polys[i][begin])
      end
      for (poly_true, poly_f) in zip(values(polys), values(fnc_polys))
            @test DT.circular_equality(filter(!DT.is_boundary_index, poly_true), poly_f)
      end
      tri, label_map, index_map = simple_geometry()
      DT.compute_representative_points!(tri)
      add_ghost_triangles!(tri)
      polys = Dict(
            1 => [[2, 20, 19, 8, DT.BoundaryIndex, 2]],
            2 => [[3, 16, 20, 1, DT.BoundaryIndex, 3]],
            3 => [[4, 17, 16, 2, DT.BoundaryIndex, 4]],
            4 => [[5, 18, 17, 3, DT.BoundaryIndex, 5]],
            5 => [[6, 22, 24, 18, 4, DT.BoundaryIndex, 6]],
            6 => [[7, 25, 9, 12, 23, 22, 5, DT.BoundaryIndex, 7]],
            7 => [[8, 26, 9, 25, 6, DT.BoundaryIndex, 8]],
            8 => [[1, 19, 10, 26, 7, DT.BoundaryIndex, 1]],
            9 => [[6, 25, 7, 26, 10, DT.BoundaryIndex - 1, 12, 6]],
            10 => [[11, DT.BoundaryIndex - 1, 9, 26, 8, 19, 20, 21, 11]],
            11 => [[12, DT.BoundaryIndex - 1, 10, 21, 23, 12]],
            12 => [[23, 6, 9, DT.BoundaryIndex - 1, 11, 23]],
            13 => [[18, 24, 22, 23, 21, 14, DT.BoundaryIndex - 2, 18],
                  [18, 24, 22, 23, 21, 14, DT.BoundaryIndex - 3, 18]],
            14 => [[13, 21, 15, DT.BoundaryIndex - 2, 13],
                  [13, 21, 15, DT.BoundaryIndex - 3, 13]],
            15 => [[14, 21, 20, 16, DT.BoundaryIndex - 2, 14],
                  [14, 21, 20, 16, DT.BoundaryIndex - 3, 14]],
            16 => [[15, 20, 2, 3, 17, DT.BoundaryIndex - 2, 15],
                  [15, 20, 2, 3, 17, DT.BoundaryIndex - 3, 15]],
            17 => [[16, 3, 4, 18, DT.BoundaryIndex - 2, 16],
                  [16, 3, 4, 18, DT.BoundaryIndex - 3, 16]],
            18 => [[17, 4, 5, 24, 13, DT.BoundaryIndex - 2, 17],
                  [17, 4, 5, 24, 13, DT.BoundaryIndex - 3, 17]],
            19 => [[1, 20, 10, 8, 1]],
            20 => [[16, 15, 21, 10, 19, 1, 2, 16]],
            21 => [[15, 14, 13, 23, 11, 10, 20, 15]],
            22 => [[24, 5, 6, 23, 13, 24]],
            23 => [[13, 22, 6, 12, 11, 21, 13]],
            24 => [[18, 5, 22, 13, 18]],
            25 => [[9, 6, 7, 9]],
            26 => [[10, 9, 7, 8, 10]],
            DT.BoundaryIndex => [[5, 4, 3, 2, 1, 8, 7, 6, 5]],
            DT.BoundaryIndex - 1 => [[9, 10, 11, 12, 9]],
            DT.BoundaryIndex - 2 => [[14, 15, 16, 17, 18, 13, 14]],
            DT.BoundaryIndex - 3 => [[14, 15, 16, 17, 18, 13, 14]]
      )
      fnc_polys = Dict{Int64,Vector{Int64}}()
      for i in keys(polys)
            fnc_polys[i] = DT.get_surrounding_polygon(tri, i)
            push!(fnc_polys[i], fnc_polys[i][begin])
      end
      for (poly_true, poly_f) in zip(values(polys), values(fnc_polys))
            @test any(DT.circular_equality(S, poly_f) for S in poly_true)
      end
      fnc_polys = Dict{Int64,Vector{Int64}}()
      for i in keys(polys)
            fnc_polys[i] = DT.get_surrounding_polygon(tri, i; skip_boundary_indices=true)
            push!(fnc_polys[i], fnc_polys[i][begin])
      end
      for (poly_true, poly_f) in zip(values(polys), values(fnc_polys))
            @test any(DT.circular_equality(filter(!DT.is_boundary_index, S), poly_f) for S in poly_true)
      end

      for _ in 1:1000
            local tri, polys, fnc_polys
            tri = example_with_special_corners()
            polys = Dict(
                  1 => [DT.BoundaryIndex, 5, 4, 16, 3, 2, DT.BoundaryIndex],
                  2 => [DT.BoundaryIndex, 1, 3, DT.BoundaryIndex],
                  3 => [2, 1, 16, 6, 7, DT.BoundaryIndex, 2],
                  4 => [5, 16, 1, 5],
                  5 => [13, 18, 16, 4, 1, DT.BoundaryIndex, 13],
                  6 => [3, 16, 17, 7, 3],
                  7 => [3, 6, 17, 15, 8, DT.BoundaryIndex, 3],
                  8 => [7, 15, 9, DT.BoundaryIndex, 7],
                  9 => [8, 15, 11, 10, DT.BoundaryIndex, 8],
                  10 => [9, 11, 12, 13, DT.BoundaryIndex, 9],
                  11 => [9, 15, 12, 10, 9],
                  12 => [18, 13, 10, 11, 15, 14, 18],
                  13 => [10, 12, 18, 5, DT.BoundaryIndex, 10],
                  14 => [15, 17, 18, 12, 15],
                  15 => [8, 7, 17, 14, 12, 11, 9, 8],
                  16 => [4, 5, 18, 17, 6, 3, 1, 4],
                  17 => [7, 6, 16, 18, 14, 15, 7],
                  18 => [5, 13, 12, 14, 17, 16, 5]
            )
            tri = triangulate(get_points(tri); delete_ghosts=false)
            fnc_polys = Dict{Int64,Vector{Int64}}()
            for i in keys(polys)
                  fnc_polys[i] = DT.get_surrounding_polygon(tri, i)
                  push!(fnc_polys[i], fnc_polys[i][begin])
            end
            for (k, S) in polys
                  @test DT.circular_equality(S, fnc_polys[k])
            end
            fnc_polys = Dict{Int64,Vector{Int64}}()
            for i in keys(polys)
                  fnc_polys[i] = DT.get_surrounding_polygon(tri, i; skip_boundary_indices=true)
                  push!(fnc_polys[i], fnc_polys[i][begin])
            end
            for (k, S) in polys
                  _S = filter(!DT.is_boundary_index, S)
                  DT.is_boundary_index(S[begin]) && push!(_S, _S[begin]) # might have removed a boundary index in the first entry, so we'd no longer have a circular array
                  @test DT.circular_equality(_S, fnc_polys[k])
            end
      end
end

@testset "sort_edge_by_degree" begin
      tri = triangulate(rand(2, 500); delete_ghosts=false)
      graph = get_graph(tri)
      for e in DT.get_edges(graph)
            new_e = DT.sort_edge_by_degree(e, graph)
            d1 = DT.num_neighbours(graph, e[1])
            d2 = DT.num_neighbours(graph, e[2])
            if d1 ≤ d2
                  @test new_e == e
            else
                  @test new_e == (e[2], e[1])
            end
      end
end