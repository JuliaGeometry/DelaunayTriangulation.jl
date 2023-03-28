using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
include("../helper_functions.jl")

tri = shewchuk_example_constrained()

fig, ax, sc = triplot(tri)
let vert = each_solid_vertex(tri)
    text!(ax, collect(get_point(tri, vert...)); text=string.(vert))
end
lines!(ax, [get_point(tri, 2, 7)...], color=:blue, linestyle=:dash)
fig


intersecting_triangles = DT.locate_intersecting_triangles(
    (2, 7),
    get_points(tri),
    get_adjacent(tri),
    get_adjacent2vertex(tri),
    get_graph(tri),
    get_boundary_index_ranges(tri),
    get_boundary_map(tri),
    NTuple{3,Int64},
    Vector{NTuple{3,Int64}}
)
@test all(T -> DT.is_positively_oriented(DT.triangle_orientation(tri, T)), intersecting_triangles)
@test all(DT.is_multiple, [DT.triangle_line_segment_intersection(tri, T, (2, 7)) for T in intersecting_triangles])
allT = [(2, 4, 3), (3, 4, 10), (10, 4, 9), (9, 5, 10), (5, 8, 10), (5, 6, 8), (8, 6, 7)]
@test DT.compare_triangle_collections(allT, intersecting_triangles)
@test allunique(intersecting_triangles)

for e in ((2, 10), (10, 2))
    intersecting_triangles1 = DT.locate_intersecting_triangles(tri, e)
    intersecting_triangles2 = DT.locate_intersecting_triangles(
        e,
        get_points(tri),
        get_adjacent(tri),
        get_adjacent2vertex(tri),
        get_graph(tri),
        get_boundary_index_ranges(tri),
        get_boundary_map(tri),
        NTuple{3,Int64},
        Vector{NTuple{3,Int64}}
    )
    for intersecting_triangles in (intersecting_triangles1, intersecting_triangles2)
        @test all(T -> DT.is_positively_oriented(DT.triangle_orientation(tri, T)), intersecting_triangles)
        @test all(DT.is_multiple, [DT.triangle_line_segment_intersection(tri, T, e) for T in intersecting_triangles])
        allT = [(2, 4, 3), (3, 4, 10)]
        @test DT.compare_triangle_collections(allT, intersecting_triangles)
        @test allunique(intersecting_triangles)
    end
end

for e in ((4, 11), (11, 4))
    intersecting_triangles1 = DT.locate_intersecting_triangles(tri, e)
    intersecting_triangles2 = DT.locate_intersecting_triangles(
        e,
        get_points(tri),
        get_adjacent(tri),
        get_adjacent2vertex(tri),
        get_graph(tri),
        get_boundary_index_ranges(tri),
        get_boundary_map(tri),
        NTuple{3,Int64},
        Vector{NTuple{3,Int64}}
    )
    for intersecting_triangles in (intersecting_triangles1, intersecting_triangles2)
        @test all(T -> DT.is_positively_oriented(DT.triangle_orientation(tri, T)), intersecting_triangles)
        @test all(DT.is_multiple, [DT.triangle_line_segment_intersection(tri, T, e) for T in intersecting_triangles])
        allT = [(4, 9, 10), (10, 9, 5), (10, 5, 8), (10, 8, 11)]
        @test DT.compare_triangle_collections(allT, intersecting_triangles)
        @test allunique(intersecting_triangles)
    end
end

e = (1, 6)
intersecting_triangles = DT.locate_intersecting_triangles(tri, e)
intersecting_triangles

V =DT.jump_and_march(tri, get_point(tri, 1); k=4) 
DT.point_position_relative_to_triangle(tri, V, get_point(tri,1))

get_point(tri,V...)