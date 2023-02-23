using ..DelaunayTriangulation
const DT = DelaunayTriangulation

include("./helper_functions.jl")

@test DT.is_true(true)
@test DT.is_true(Val(true))
@test !DT.is_true(false)
@test !DT.is_true(Val(false))
@test DT.is_true(Val{true})
@test !DT.is_true(Val{false})

@test number_type([1, 2, 3]) == Int64
@test number_type([1.0, 2.0, 3.0]) == Float64
@test number_type([1.0 2.0; 3.0 3.5; 10.0 17.3]) == Float64
@test number_type((1.0, 5.0)) == Float64
@test number_type([(1.0f0, 2.0f0), (1.7f0, 2.5f0)]) == Float32
@test number_type(2.4) == Float64

@test DT.get_boundary_index(1, 2, -3) == -3
@test DT.get_boundary_index(1, 2, -1) == -1
@test DT.get_boundary_index(1, -5, 2) == -5
@test DT.get_boundary_index(-1, 2, 3) == -1
@test_throws ArgumentError DT.get_boundary_index(2, 5, 7)
@test DT.get_boundary_index(1, -2) == -2
@test DT.get_boundary_index(-5, 1) == -5
@test_throws ArgumentError DT.get_boundary_index(2, 5)

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

i, j, k = rand(Int64, 3)
@test DT.choose_uvw(true, false, false, i, j, k) == (i, j, k)
@test DT.choose_uvw(false, true, false, i, j, k) == (j, k, i)
@test DT.choose_uvw(false, false, true, i, j, k) == (k, i, j)

x = rand(10)
@test !DT.is_circular(x)
push!(x, x[begin])
@test DT.is_circular(x)

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
