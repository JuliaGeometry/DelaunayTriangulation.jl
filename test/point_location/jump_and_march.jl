using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using LinearAlgebra
using Random
using StableRNGs
using StatsBase

include("../helper_functions.jl")

tri, label_map, index_map = simple_geometry()
add_ghost_triangles!(tri)
DT.compute_representative_points!(tri)
pts = get_points(tri)
adj = get_adjacent(tri)
adj2v = get_adjacent2vertex(tri)
boundary_index_ranges = get_boundary_index_ranges(tri)
boundary_map = get_boundary_map(tri)
graph = get_graph(tri)
boundary_nodes = get_boundary_nodes(tri)

@testset "Specific point and triangle examples" begin
    q1 = (2.0, 4.0)
    V1 = [(index_map["a"], index_map["s"], index_map["h"])]
    q2 = (1.1803475764382, 17.4209559058887)
    V2 = [(index_map["g"], index_map["b1"], index_map["i"])]
    q3 = (14.2726546767168, 18.7727571005127)
    V3 = [(index_map["f"], index_map["v"], index_map["e"])]
    q4 = (12.0, 2.0)
    V4 = [(index_map["t"], index_map["p"], index_map["o"]),
        (index_map["t"], index_map["b"], index_map["p"])]
    q5 = (13.0, 0.0)
    V5 = [(index_map["b"], index_map["c"], index_map["p"]),
        (index_map["c"], index_map["b"], DT.BoundaryIndex)]
    q6 = [6.0819947177817, 8.90410894457]
    V6 = [(index_map["j"], index_map["k"], DT.BoundaryIndex - 1)]
    q7 = (15.8, 9.2)
    V7 = [(index_map["m"], DT.BoundaryIndex - 2, index_map["r"]),
        (index_map["m"], DT.BoundaryIndex - 3, index_map["r"])]
    q8 = (14.5391680629781, 8.5135910023477)
    V8 = [(index_map["m"], index_map["n"], DT.BoundaryIndex - 2),
        (index_map["m"], index_map["n"], DT.BoundaryIndex - 3)]
    q9 = (22.0, 8.0)
    V9 = [(index_map["d"], index_map["c"], DT.BoundaryIndex)]
    q10 = (6.0, 11.0)
    V10 = [(index_map["j"], index_map["k"], DT.BoundaryIndex - 1),
        (index_map["k"], index_map["ℓ"], DT.BoundaryIndex - 1),
        (index_map["ℓ"], index_map["i"], DT.BoundaryIndex - 1),
        (index_map["i"], index_map["j"], DT.BoundaryIndex - 1)]
    q11 = (6.0819947177817, 8.90410894457)
    V11 = [(index_map["j"], index_map["k"], DT.BoundaryIndex - 1)]
    q12 = (-6.6465702688004, 19.4798146355492)
    V12 = [(index_map["h"], index_map["g"], DT.BoundaryIndex)]
    q13 = (-20.0, 10.0)
    V13 = [(index_map["h"], index_map["g"], DT.BoundaryIndex),
        (index_map["a"], index_map["h"], DT.BoundaryIndex)]
    q14 = (20.0, 6.0)
    V14 = [(index_map["c"], index_map["d"], index_map["q"]),
        (index_map["c"], index_map["d"], DT.BoundaryIndex)]
    allq = [q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14]
    allV = [V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14]
    for _ in 1:36
        for k in each_point_index(pts)
            for i in eachindex(allq)
                T1 = jump_and_march(pts, adj, adj2v, graph, boundary_index_ranges, boundary_map,
                    allq[i]; k)
                T2 = jump_and_march(pts, adj, adj2v, graph, boundary_index_ranges, boundary_map,
                    allq[i])
                @test DT.is_positively_oriented(DT.triangle_orientation(tri, T1))
                @test DT.is_positively_oriented(DT.triangle_orientation(tri, T2))
                @inferred jump_and_march(pts, adj, adj2v, graph, boundary_index_ranges,
                    boundary_map, allq[i]; k)
                @inferred jump_and_march(pts, adj, adj2v, graph, boundary_index_ranges,
                    boundary_map, allq[i])
                for V in [T1, T2]
                    if length(allV[i]) == 1
                        @test DT.compare_triangles(V, allV[i][1]) &&
                              DT.is_inside(DT.point_position_relative_to_triangle(tri, V,
                            allq[i]))
                    elseif length(allV[i]) == 2
                        @test (DT.compare_triangles(V, allV[i][1]) ||
                               DT.compare_triangles(V, allV[i][2])) &&
                              (DT.is_on(DT.point_position_relative_to_triangle(tri, V, allq[i])) ||
                               (DT.is_ghost_triangle(V) &&
                                DT.is_inside(DT.point_position_relative_to_triangle(tri, V,
                            allq[i]))))
                    else
                        bool1 = any(j -> DT.compare_triangles(V, allV[i][j]),
                            eachindex(allV[i]))
                        bool2 = (DT.is_on(DT.point_position_relative_to_triangle(tri, V,
                            allq[i])) ||
                                 (DT.is_ghost_triangle(V) &&
                                  DT.is_inside(DT.point_position_relative_to_triangle(tri, V,
                            allq[i]))))
                        @test bool1 && bool2
                    end
                end
            end
        end
    end
end

DT.RepresentativePointList[1].x = 10.0
DT.RepresentativePointList[1].y = 10.0
_pts = tri.points[[12, 11, 10, 9]]
DT.RepresentativePointList[2].x = mean([8.0, 8.0, 4.0, 4.0])
DT.RepresentativePointList[2].y = mean([16.0, 6.0, 6.0, 16.0])
_pts = tri.points[[18, 17, 16, 15, 14, 13]]
DT.RepresentativePointList[3].x = mean([18.0, 18.0, 14.0, 12.0, 14.0, 14.0])
DT.RepresentativePointList[3].y = mean([12.0, 6.0, 2.0, 4.0, 6.0, 10.0])

@testset "Tests with different types of triangulations" begin
    x, y = complicated_geometry()
    tri2 = generate_mesh(x, y, 2.0; convert_result=true, add_ghost_triangles=true)

    a, b, c, d = 2.0, 10.0, -5.0, 7.5
    nx = 20
    ny = 10
    tri3 = DT.triangulate_rectangle(a, b, c, d, nx, ny)

    for tri in (tri, tri2, tri3)
        DT.compute_representative_points!(tri)
        if !(tri === tri2 || tri === tri3)
            local _pts
            DT.RepresentativePointList[1].x = 10.0
            DT.RepresentativePointList[1].y = 10.0
            _pts = tri.points[[12, 11, 10, 9]]
            DT.RepresentativePointList[2].x = mean([8.0, 8.0, 4.0, 4.0])
            DT.RepresentativePointList[2].y = mean([16.0, 6.0, 6.0, 16.0])
            _pts = tri.points[[18, 17, 16, 15, 14, 13]]
            DT.RepresentativePointList[3].x = mean([18.0, 18.0, 14.0, 12.0, 14.0, 14.0])
            DT.RepresentativePointList[3].y = mean([12.0, 6.0, 2.0, 4.0, 6.0, 10.0])
        end

        @testset "Test that we can find a point in every triangle" begin
            for _ in 1:36
                for V in each_triangle(tri.triangles)
                    if !DT.is_outer_ghost_triangle(indices(V)..., tri.boundary_map)
                        i, j, k = indices(V)
                        p, q, r = get_point(tri.points, tri.boundary_map, i, j, k)
                        local c
                        c = (p .+ q .+ r) ./ 3
                        for k in each_point_index(tri.points)
                            _V1 = jump_and_march(tri.points, tri.adjacent, tri.adjacent2vertex,
                                tri.graph, tri.boundary_index_ranges,
                                tri.boundary_map, c; k)
                            _V2 = jump_and_march(tri, c; k)
                            for _V in [_V1, _V2]
                                @test DT.is_positively_oriented(DT.triangle_orientation(tri, _V))
                                if !DT.is_ghost_triangle(_V...)
                                    @test DT.compare_triangles(_V, V) &&
                                          DT.is_inside(DT.point_position_relative_to_triangle(tri,
                                        _V,
                                        c))
                                else
                                    local V1, V2
                                    V1 = DT.rotate_ghost_triangle_to_standard_form(V)
                                    V2 = DT.rotate_ghost_triangle_to_standard_form(_V)
                                    i1 = geti(V1)
                                    i2 = geti(V2)
                                    if i1 ≠ i2
                                        i1 = i1 - 1
                                    end
                                    if i1 ≠ i2
                                        i1 = i1 + 1
                                    end
                                    if i1 ≠ i2
                                        @test false
                                    end
                                    _V = DT.construct_triangle(typeof(V), i1, getj(V1), getk(V1))
                                    @test DT.compare_triangles(_V, V) &&
                                          DT.is_inside(DT.point_position_relative_to_triangle(tri,
                                        _V,
                                        c))
                                end
                            end
                        end
                    end
                end
            end

            @testset "Test that we don't break for points already in the triangulation" begin
                for _ in 1:36
                    for k in each_point_index(tri.points)
                        for j in each_point_index(tri.points)
                            _V1 = jump_and_march(tri, get_point(tri, k); k=j)
                            _V2 = jump_and_march(tri.points, tri.adjacent, tri.adjacent2vertex,
                                tri.graph, tri.boundary_index_ranges, tri.boundary_map,
                                get_point(tri.points, k))
                            for _V in [_V1, _V2]
                                @test k ∈ indices(_V)
                                @test DT.is_positively_oriented(DT.triangle_orientation(tri, _V))
                            end
                        end
                    end
                end
            end

            @testset "Finding points in ghost triangles" begin
                # Technically this will also find points in solid triangles, but by doing random testing with large random points, 
                # we ensure that we primarily find ghost triangles
                for _ in 1:36
                    q = (50randn(), 50rand())
                    for k in each_point_index(tri.points)
                        _V1 = jump_and_march(tri, q; k)
                        _V2 = jump_and_march(tri.points, tri.adjacent, tri.adjacent2vertex, tri.graph,
                            tri.boundary_index_ranges, tri.boundary_map, q)
                        @test DT.is_inside(DT.point_position_relative_to_triangle(tri, _V1, q))
                        @test DT.is_inside(DT.point_position_relative_to_triangle(tri, _V2, q))
                    end
                end
            end
        end
    end
end