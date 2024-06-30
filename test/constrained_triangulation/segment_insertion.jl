using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using StableRNGs
# Base.include(@__MODULE__, "../helper_functions.jl")

#=
We test constrained Delaunay triangulations by 
comparing results to MATLAB. For example, 
for fixed_shewchuk_example_constrained,
we use the MATLAB code for generating 
a triangulation with a constrained edge 
(2, 7):

    P = [
        0.0 0.0;
        0.0 1.0;
        0.0 2.5;
        2.0 0.0;
        6.0 0.0;
        8.0 0.0;
        8.0 0.5;
        7.5 1.0;
        4.0 1.0;
        4.0 2.5;
        8.0 2.5
    ];
    C = [2 7];
    tri = delaunayTriangulation(P, C);
    tri.ConnectivityList
    ans =

     1     4     2
     8    11    10
     8     7    11
     2     9     3
     3     9    10
     9     8    10
     9     7     8
     2     7     9
     7     5     6
     7     2     5
     5     2     4
=#

@testset "is_vertex_closer_than_neighbours" begin
    tri = fixed_shewchuk_example_constrained()
    @test !DT.is_vertex_closer_than_neighbours(tri, 2, 7, 10, 8, 9)
    @test !DT.is_vertex_closer_than_neighbours(tri, 2, 7, 10, 9, 8)
    @test DT.is_vertex_closer_than_neighbours(tri, 2, 7, 9, 8, 3)
    @test DT.is_vertex_closer_than_neighbours(tri, 2, 7, 9, 3, 8)
    @test !DT.is_vertex_closer_than_neighbours(tri, 2, 7, 8, 9, 11)
    @test DT.is_vertex_closer_than_neighbours(tri, 2, 7, 3, 10, 11)
    @test !DT.is_vertex_closer_than_neighbours(tri, 7, 2, 1, 4, 5)
    @test !DT.is_vertex_closer_than_neighbours(tri, 7, 2, 4, 1, 5)
    @test DT.is_vertex_closer_than_neighbours(tri, 7, 2, 5, 1, 4)
    @test DT.is_vertex_closer_than_neighbours(tri, 7, 2, 6, 5, 4)
end

@testset "Inserting segments into the Shewchuk example" begin
    tri = fixed_shewchuk_example_constrained()
    e = (2, 7)
    T, C, L, R = DT.locate_intersecting_triangles(tri, e)
    DT.delete_intersected_triangles!(tri, T)
    points = get_points(tri)
    _tri_1 = tri.cache
    DT.triangulate_cavity_cdt!(_tri_1.triangulation, L, _tri_1.triangulation_2, _tri_1.marked_vertices,_tri_1.fan_triangles)
    DT.add_new_triangles!(tri, _tri_1.triangulation)
    empty!(_tri_1)
    DT.triangulate_cavity_cdt!(_tri_1.triangulation, R, _tri_1.triangulation_2, _tri_1.marked_vertices,_tri_1.fan_triangles)
    DT.add_new_triangles!(tri, _tri_1.triangulation)
    empty!(_tri_1)
    true_tri = ([ # two triangles to account for the cocircular points (2, 9, 3, 10)
            (1, 4, 2)
            (8, 11, 10)
            (8, 7, 11)
            (2, 9, 3)
            (3, 9, 10)
            (9, 8, 10)
            (9, 7, 8)
            (2, 7, 9)
            (7, 5, 6)
            (7, 2, 5)
            (5, 2, 4)
            (2, 3, DT.𝒢)
            (10, 11, DT.𝒢)
            (3, 10, DT.𝒢)
            (11, 7, DT.𝒢)
            (7, 6, DT.𝒢)
            (6, 5, DT.𝒢)
            (5, 4, DT.𝒢)
            (4, 1, DT.𝒢)
            (1, 2, DT.𝒢)
        ],
        [
            (1, 4, 2)
            (8, 11, 10)
            (8, 7, 11)
            (2, 9, 10)
            (3, 2, 10)
            (9, 8, 10)
            (9, 7, 8)
            (2, 7, 9)
            (7, 5, 6)
            (7, 2, 5)
            (5, 2, 4)
            (2, 3, DT.𝒢)
            (3, 10, DT.𝒢)
            (10, 11, DT.𝒢)
            (11, 7, DT.𝒢)
            (7, 6, DT.𝒢)
            (6, 5, DT.𝒢)
            (5, 4, DT.𝒢)
            (4, 1, DT.𝒢)
            (1, 2, DT.𝒢)
        ]
    )
    @test any(T -> DT.compare_triangle_collections(get_triangles(tri), T), true_tri)
    push!(get_all_segments(tri), e)
    push!(get_interior_segments(tri), e)
    @test validate_triangulation(tri)

    e = (2, 11)
    T, C, L, R = DT.locate_intersecting_triangles(tri, e)
    DT.delete_intersected_triangles!(tri, T)
    points = get_points(tri)
    _tri_1 = tri.cache
    DT.triangulate_cavity_cdt!(_tri_1.triangulation, L, _tri_1.triangulation_2, _tri_1.marked_vertices,_tri_1.fan_triangles)
    DT.add_new_triangles!(tri, _tri_1.triangulation)
    empty!(_tri_1)
    DT.triangulate_cavity_cdt!(_tri_1.triangulation, R, _tri_1.triangulation_2, _tri_1.marked_vertices,_tri_1.fan_triangles)
    DT.add_new_triangles!(tri, _tri_1.triangulation)
    empty!(_tri_1)
    true_tri = [
        (1, 4, 2)
        (2, 10, 3)
        (2, 11, 10)
        (11, 9, 8)
        (11, 2, 9)
        (8, 7, 11)
        (9, 7, 8)
        (2, 7, 9)
        (7, 5, 6)
        (7, 2, 5)
        (5, 2, 4)
        (2, 3, DT.𝒢)
        (10, 11, DT.𝒢)
        (3, 10, DT.𝒢)
        (11, 7, DT.𝒢)
        (7, 6, DT.𝒢)
        (6, 5, DT.𝒢)
        (5, 4, DT.𝒢)
        (4, 1, DT.𝒢)
        (1, 2, DT.𝒢)
    ]
    @test DT.compare_triangle_collections(get_triangles(tri), true_tri)
    push!(get_all_segments(tri), e)
    push!(get_interior_segments(tri), e)
    @test validate_triangulation(tri)

    e = (1, 7)
    T, C, L, R = DT.locate_intersecting_triangles(tri, e)
    DT.delete_intersected_triangles!(tri, T)
    _tri_1 = tri.cache
    DT.triangulate_cavity_cdt!(_tri_1.triangulation, L, _tri_1.triangulation_2, _tri_1.marked_vertices,_tri_1.fan_triangles)
    DT.add_new_triangles!(tri, _tri_1.triangulation)
    empty!(_tri_1)
    DT.triangulate_cavity_cdt!(_tri_1.triangulation, R, _tri_1.triangulation_2, _tri_1.marked_vertices,_tri_1.fan_triangles)
    DT.add_new_triangles!(tri, _tri_1.triangulation)
    empty!(_tri_1)
    true_tri = [
        (2, 10, 3)
        (1, 7, 2)
        (7, 4, 5)
        (7, 1, 4)
        (11, 2, 9)
        (8, 7, 11)
        (9, 7, 8)
        (2, 7, 9)
        (7, 5, 6)
        (11, 9, 8)
        (2, 11, 10)
        (2, 3, DT.𝒢)
        (10, 11, DT.𝒢)
        (3, 10, DT.𝒢)
        (11, 7, DT.𝒢)
        (7, 6, DT.𝒢)
        (6, 5, DT.𝒢)
        (5, 4, DT.𝒢)
        (4, 1, DT.𝒢)
        (1, 2, DT.𝒢)
    ]
    @test DT.compare_triangle_collections(get_triangles(tri), true_tri)
    push!(get_all_segments(tri), e)
    push!(get_interior_segments(tri), e)
    @test validate_triangulation(tri)
end

@testset "Testing the special corner example with many edges" begin
    tri = example_with_special_corners()
    es = [(9, 13), (9, 12), (9, 18), (9, 15), (15, 3), (16, 2)]
    for e in es
        if !DT.edge_exists(tri, e)
            T, C, L, R = DT.locate_intersecting_triangles(tri, e)
            DT.delete_intersected_triangles!(tri, T)
            points = get_points(tri)
            _tri_1 = tri.cache
            DT.triangulate_cavity_cdt!(_tri_1.triangulation, L, _tri_1.triangulation_2, _tri_1.marked_vertices,_tri_1.fan_triangles)
            DT.add_new_triangles!(tri, _tri_1.triangulation)
            empty!(_tri_1)
            DT.triangulate_cavity_cdt!(_tri_1.triangulation, R, _tri_1.triangulation_2, _tri_1.marked_vertices,_tri_1.fan_triangles)
            DT.add_new_triangles!(tri, _tri_1.triangulation)
            empty!(_tri_1)
        end
        push!(get_all_segments(tri), e)
        push!(get_interior_segments(tri), e)
    end
    true_tri = [
        (14, 15, 17)
        (14, 9, 15)
        (15, 8, 7)
        (9, 13, 11)
        (17, 18, 14)
        (4, 1, 5)
        (15, 9, 8)
        (11, 10, 9)
        (9, 18, 12)
        (13, 9, 12)
        (15, 3, 17)
        (6, 15, 7)
        (12, 18, 13)
        (18, 5, 13)
        (18, 16, 5)
        (5, 16, 4)
        (18, 17, 16)
        (17, 3, 16)
        (6, 7, 3)
        (11, 13, 10)
        (16, 1, 4)
        (16, 2, 1)
        (2, 16, 3)
        (18, 9, 14)
        (3, 15, 6)
        (10, 13, DT.𝒢)
        (13, 5, DT.𝒢)
        (5, 1, DT.𝒢)
        (1, 2, DT.𝒢)
        (2, 3, DT.𝒢)
        (3, 7, DT.𝒢)
        (7, 8, DT.𝒢)
        (8, 9, DT.𝒢)
        (9, 10, DT.𝒢)
    ]
    @test DT.compare_triangle_collections(get_triangles(tri), true_tri)
    @test validate_triangulation(tri)
end