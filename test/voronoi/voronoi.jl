using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using ColorSchemes
using DataStructures
using LinearAlgebra

include("../helper_functions.jl")

@testset "Unconstrained test" begin
    A = (-1.0, 7.0)
    B = (4.0, 4.0)
    C = (-2.0, -1.0)
    D = (-1.0, 3.0)
    E = (3.0, -1.0)
    F = (1.0, 4.0)
    G = (-3.0, 5.0)
    pts = [A, B, C, D, E, F, G]
    tri = triangulate(pts; delete_ghosts=false, randomise=false)

    #fig, ax, sc = triplot(tri)
    #let vert = each_solid_vertex(tri)
    #    text!(ax, collect(get_point(tri, vert...)); text=string.(vert))
    #end
    #fig

    vorn = voronoi(tri)
    for (i, p) in DT.get_generators(vorn)
        @test get_point(tri, i) == get_generator(vorn, i) == p
    end
    @test validate_tessellation(vorn)
    @test DT.get_triangulation(vorn) == tri
    circumcenter_to_triangle = DT.get_circumcenter_to_triangle(vorn)
    triangle_to_circumcenter = DT.get_triangle_to_circumcenter(vorn)
    for V in DT.each_solid_triangle(DT.get_triangulation(vorn))
        V = DT.rotate_triangle_to_standard_form(V)
        c = DT.get_triangle_to_circumcenter(vorn, V)
        c = get_polygon_point(vorn, c)
        i, j, k = indices(V)
        p, q, r = get_point(DT.get_triangulation(vorn), i, j, k)
        cx, cy = DT.triangle_circumcenter(p, q, r)
        @test cx == c[1] && cy == c[2]
    end
    for (c, V) in circumcenter_to_triangle
        @test DT.get_circumcenter_to_triangle(vorn, c) == V
        @test DT.get_triangle_to_circumcenter(vorn, V) == c
    end
    for (V, c) in triangle_to_circumcenter
        @test DT.get_circumcenter_to_triangle(vorn, c) == V
        @test DT.get_triangle_to_circumcenter(vorn, V) == c
    end
    @test DT.circular_equality(
        get_polygon(vorn, 1),
        DT.get_triangle_to_circumcenter.(Ref(vorn), [
            (7, 4, 1), (4, 6, 1), (6, 2, 1), (1, 2, -1), (7, 1, -1), (7, 4, 1)
        ])
    )
    @test DT.circular_equality(
        get_polygon(vorn, 2),
        DT.get_triangle_to_circumcenter.(Ref(vorn), [
            (2, 5, -1), (1, 2, -1), (6, 2, 1), (6, 5, 2), (2, 5, -1)
        ])
    )
    @test DT.circular_equality(
        get_polygon(vorn, 3),
        DT.get_triangle_to_circumcenter.(Ref(vorn), [
            (5, 3, -1), (5, 4, 3), (4, 7, 3), (3, 7, -1), (5, 3, -1)
        ])
    )
    @test DT.circular_equality(
        get_polygon(vorn, 4),
        DT.get_triangle_to_circumcenter.(Ref(vorn), [
            (5, 6, 4), (4, 6, 1), (7, 4, 1), (4, 7, 3), (5, 4, 3), (5, 6, 4)
        ])
    )
    @test DT.circular_equality(
        get_polygon(vorn, 5),
        DT.get_triangle_to_circumcenter.(Ref(vorn), [
            (6, 5, 2), (5, 6, 4), (5, 4, 3), (5, 3, -1), (2, 5, -1), (6, 5, 2)
        ])
    )
    @test DT.circular_equality(
        get_polygon(vorn, 6),
        DT.get_triangle_to_circumcenter.(Ref(vorn), [
            (6, 5, 2), (6, 2, 1), (4, 6, 1), (5, 6, 4), (6, 5, 2)
        ])
    )
    @test DT.circular_equality(
        get_polygon(vorn, 7),
        DT.get_triangle_to_circumcenter.(Ref(vorn), [
            (7, 4, 1), (7, 1, -1), (3, 7, -1), (4, 7, 3), (7, 4, 1)
        ])
    )

    xmin, xmax, ymin, ymax = DT.polygon_bounds(vorn, 0.5)
    cmap = Makie.cgrad(:jet)
    colors = get_polygon_colors(vorn, cmap)
    fig, ax, sc = voronoiplot(vorn, polygon_color=colors)
    triplot!(ax, tri)
    xlims!(ax, xmin, xmax)
    ylims!(ax, ymin, ymax)
    fig
end

@testset "Voronoi point location" begin
    A = (-1.0, 7.0)
    B = (4.0, 4.0)
    C = (-2.0, -1.0)
    D = (-1.0, 3.0)
    E = (3.0, -1.0)
    F = (1.0, 4.0)
    G = (-3.0, 5.0)
    pts = [A, B, C, D, E, F, G]
    tri = triangulate(pts; delete_ghosts=false, randomise=false)
    vor = voronoi(tri)
    @test validate_tessellation(vor)
    xmin, xmax, ymin, ymax = DT.polygon_bounds(get_points(tri), get_convex_hull_indices(tri))
    p = NTuple{2,Float64}[]
    n = 100000
    while length(p) ≤ n # only going to test points that are inside the polygon
        pt = (xmin + rand() * (xmax - xmin), ymin + rand() * (ymax - ymin))
        if DT.distance_to_polygon(pt, get_points(tri), get_convex_hull_indices(tri)) ≥ 0
            push!(p, pt)
        end
    end
    for p in p
        u = jump_and_march(vor, p)
        all_dists = [norm(p .- get_generator(vor, i)) for i in sort(collect(each_generator(vor)))]
        @test findmin(all_dists)[2] == u
    end
end

@testset "Clipping a simple VoronoiTessellation" begin
    for _ in 1:1000
        A = (-1.0, 7.0)
        B = (4.0, 4.0)
        C = (-2.0, -1.0)
        D = (-1.0, 3.0)
        E = (3.0, -1.0)
        F = (1.0, 4.0)
        G = (-3.0, 5.0)
        pts = [A, B, C, D, E, F, G]
        tri = triangulate(pts; delete_ghosts=false, randomise=true)
        triplot(tri)
        #lock_convex_hull!(tri)
        vorn = voronoi(tri, true)
        for (i, p) in DT.get_generators(vorn)
            @test get_point(tri, i) == get_generator(vorn, i) == p
        end
        @test DT.get_triangulation(vorn) == tri
        circumcenter_to_triangle = DT.get_circumcenter_to_triangle(vorn)
        triangle_to_circumcenter = DT.get_triangle_to_circumcenter(vorn)
        for V in DT.each_solid_triangle(DT.get_triangulation(vorn))
            V = DT.rotate_triangle_to_standard_form(V)
            c = DT.get_triangle_to_circumcenter(vorn, V)
            c = get_polygon_point(vorn, c)
            i, j, k = indices(V)
            p, q, r = get_point(DT.get_triangulation(vorn), i, j, k)
            cx, cy = DT.triangle_circumcenter(p, q, r)
            @test cx == c[1] && cy == c[2]
        end
        for (c, V) in circumcenter_to_triangle
            @test DT.get_circumcenter_to_triangle(vorn, c) == V
            @test DT.get_triangle_to_circumcenter(vorn, V) == c
        end
        for (V, c) in triangle_to_circumcenter
            @test DT.get_circumcenter_to_triangle(vorn, c) == V
            @test DT.get_triangle_to_circumcenter(vorn, V) == c
        end
        voronoiplot(vorn)
        orig_pt = collect.([(0.5, 0.5)
            (2.5, 7.166666666666667)
            (1.1666666666666665, 1.1666666666666667)
            (-4.3, 1.7000000000000002)
            (-1.0, 5.0)
            (-0.75, 5.0)
            (2.5, 1.7000000000000002)
            (0.5, -1.0)
            (3.5, 1.5)
            (3.0, -1.0)
            (-2.7142857142857144, 3.2857142857142856)
            (-2.369565217391304, 1.2173913043478262)
            (2.5000000000000004, 4.8999999999999995)
            (0.7105263157894739, 5.973684210526315)
            (-2.0, 6.0)
            (-3.0, 5.0)
            (4.0, 4.0)
            (-2.0, -1.0)
            (-1.0, 7.0)])
        @test isempty(DT.get_unbounded_polygons(vorn))
        @test all([0 ≤ get_area(vorn, i) < Inf for i in each_polygon_index(vorn)])
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(vorn), get_polygon(vorn, 5))), getindex.(Ref(orig_pt), [8, 10, 9, 7, 3, 1, 8]), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(vorn), get_polygon(vorn, 4))), getindex.(Ref(orig_pt), [12, 1, 3, 6, 5, 11, 12]), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(vorn), get_polygon(vorn, 6))), getindex.(Ref(orig_pt), [3, 7, 13, 14, 6, 3]), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(vorn), get_polygon(vorn, 7))), getindex.(Ref(orig_pt), [11, 5, 15, 16, 11]), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(vorn), get_polygon(vorn, 2))), getindex.(Ref(orig_pt), [7, 9, 17, 13, 7]), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(vorn), get_polygon(vorn, 3))), getindex.(Ref(orig_pt), [18, 8, 1, 12, 18]), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(vorn), get_polygon(vorn, 1))), getindex.(Ref(orig_pt), [5, 6, 14, 19, 15, 5]), ≈)
        @test allunique(DT.get_polygon_points(vorn))
        for i in each_polygon_index(vorn)
            C = get_polygon(vorn, i)
            for (j, v) in pairs(C)
                δ = DT.distance_to_polygon(get_polygon_point(vorn, v), get_points(tri), get_convex_hull_indices(tri))
                @test δ ≥ 0
            end
        end
    end
end

@testset "Another example" begin
    for _ in 1:100
        tri = fixed_shewchuk_example_constrained()
        vorn = voronoi(tri, false)
        @test validate_tessellation(vorn)
        for (i, p) in DT.get_generators(vorn)
            @test get_point(tri, i) == get_generator(vorn, i) == p
        end
        @test DT.get_triangulation(vorn) == tri
        @test DT.get_unbounded_polygons(vorn) == Set((3, 10, 11, 7, 6, 5, 4, 1, 2))
        circumcenter_to_triangle = DT.get_circumcenter_to_triangle(vorn)
        triangle_to_circumcenter = DT.get_triangle_to_circumcenter(vorn)
        for V in DT.each_solid_triangle(DT.get_triangulation(vorn))
            V = DT.rotate_triangle_to_standard_form(V)
            c = DT.get_triangle_to_circumcenter(vorn, V)
            c = get_polygon_point(vorn, c)
            i, j, k = indices(V)
            p, q, r = get_point(DT.get_triangulation(vorn), i, j, k)
            cx, cy = DT.triangle_circumcenter(p, q, r)
            @test cx == c[1] && cy == c[2]
        end
        for (c, V) in circumcenter_to_triangle
            @test DT.get_circumcenter_to_triangle(vorn, c) == V
            @test DT.get_triangle_to_circumcenter(vorn, V) == c
        end
        for (V, c) in triangle_to_circumcenter
            @test DT.get_circumcenter_to_triangle(vorn, c) == V
            @test DT.get_triangle_to_circumcenter(vorn, V) == c
        end
        @test DT.circular_equality(
            get_polygon(vorn, 1),
            DT.get_triangle_to_circumcenter.(Ref(vorn), [
                (4, 2, 1), (1, 2, -1), (4, 1, -1), (4, 2, 1)
            ])
        )
        @test DT.circular_equality(
            get_polygon(vorn, 2),
            DT.get_triangle_to_circumcenter.(Ref(vorn), [
                (4, 2, 1), (4, 3, 2), (2, 3, -1), (1, 2, -1), (4, 2, 1)
            ])
        )
        @test DT.circular_equality(
            get_polygon(vorn, 3),
            DT.get_triangle_to_circumcenter.(Ref(vorn), [
                (4, 3, 2), (4, 10, 3), (3, 10, -1), (2, 3, -1), (4, 3, 2)
            ])
        )
        @test DT.circular_equality(
            get_polygon(vorn, 4),
            DT.get_triangle_to_circumcenter.(Ref(vorn), [
                (5, 9, 4), (9, 10, 4), (4, 10, 3), (4, 3, 2), (4, 2, 1), (4, 1, -1), (5, 4, -1), (5, 9, 4)
            ])
        )
        @test DT.circular_equality(
            get_polygon(vorn, 5),
            DT.get_triangle_to_circumcenter.(Ref(vorn), [
                (6, 8, 5), (8, 10, 5), (10, 9, 5), (5, 9, 4), (5, 4, -1), (6, 5, -1), (6, 8, 5)
            ])
        )
        @test DT.circular_equality(
            get_polygon(vorn, 6),
            DT.get_triangle_to_circumcenter.(Ref(vorn), [
                (7, 8, 6), (6, 8, 5), (6, 5, -1), (7, 6, -1), (7, 8, 6)
            ])
        )
        @test DT.circular_equality(
            get_polygon(vorn, 7),
            DT.get_triangle_to_circumcenter.(Ref(vorn), [
                (11, 8, 7), (7, 8, 6), (7, 6, -1), (11, 7, -1), (11, 8, 7)
            ])
        )
        @test DT.circular_equality(
            get_polygon(vorn, 8),
            DT.get_triangle_to_circumcenter.(Ref(vorn), [
                (7, 8, 6), (11, 8, 7), (11, 10, 8), (8, 10, 5), (6, 8, 5), (7, 8, 6)
            ])
        )
        @test DT.circular_equality(
            get_polygon(vorn, 9),
            DT.get_triangle_to_circumcenter.(Ref(vorn), [
                (10, 9, 5), (9, 10, 4), (5, 9, 4), (10, 9, 5)
            ])
        )
        @test DT.circular_equality(
            get_polygon(vorn, 10),
            DT.get_triangle_to_circumcenter.(Ref(vorn), [
                (9, 10, 4), (10, 9, 5), (8, 10, 5), (11, 10, 8), (10, 11, -1), (3, 10, -1), (4, 10, 3), (9, 10, 4)
            ])
        )

        _vorn = voronoi(tri, true)
        @test validate_tessellation(_vorn)
        for (i, p) in DT.get_generators(_vorn)
            @test get_point(tri, i) == get_generator(_vorn, i) == p
        end
        @test DT.get_triangulation(_vorn) == tri
        circumcenter_to_triangle = DT.get_circumcenter_to_triangle(_vorn)
        triangle_to_circumcenter = DT.get_triangle_to_circumcenter(_vorn)
        for V in DT.each_solid_triangle(DT.get_triangulation(_vorn))
            V = DT.rotate_triangle_to_standard_form(V)
            c = DT.get_triangle_to_circumcenter(_vorn, V)
            c = get_polygon_point(_vorn, c)
            i, j, k = indices(V)
            p, q, r = get_point(DT.get_triangulation(_vorn), i, j, k)
            cx, cy = DT.triangle_circumcenter(p, q, r)
            @test cx == c[1] && cy == c[2]
        end
        for (c, V) in circumcenter_to_triangle
            @test DT.get_circumcenter_to_triangle(_vorn, c) == V
            @test DT.get_triangle_to_circumcenter(_vorn, V) == c
        end
        for (V, c) in triangle_to_circumcenter
            @test DT.get_circumcenter_to_triangle(_vorn, c) == V
            @test DT.get_triangle_to_circumcenter(_vorn, V) == c
        end
        @test DT.circular_equality(get_polygon(_vorn, 1), [19, 18, 11, 17, 19])
        @test DT.circular_equality(get_polygon(_vorn, 2), [30, 17, 11, 7, 26, 30])
        @test DT.circular_equality(get_polygon(_vorn, 3), [26, 7, 4, 25, 27, 26])
        @test DT.circular_equality(get_polygon(_vorn, 4), [11, 18, 29, 24, 2, 4, 7, 11])
        @test DT.circular_equality(get_polygon(_vorn, 5), [13, 14, 12, 9, 5, 3, 13])
        @test DT.circular_equality(get_polygon(_vorn, 6), [12, 21, 20, 1, 9, 12])
        @test DT.circular_equality(get_polygon(_vorn, 7), [1, 20, 28, 15, 1])
        @test DT.circular_equality(get_polygon(_vorn, 8), [9, 1, 15, 16, 6, 5, 9])
        @test DT.circular_equality(get_polygon(_vorn, 9), [24, 13, 3, 2, 24])
        @test DT.circular_equality(get_polygon(_vorn, 10), [4, 2, 3, 5, 6, 22, 31, 25, 4])
        @test isempty(DT.get_unbounded_polygons(_vorn))
        @test all([0 ≤ get_area(_vorn, i) < Inf for i in each_polygon_index(_vorn)])
        @test allunique(DT.get_polygon_points(_vorn))
        for i in each_polygon_index(_vorn)
            C = get_polygon(_vorn, i)
            for (j, v) in pairs(C)
                δ = DT.distance_to_polygon(get_polygon_point(_vorn, v), get_points(tri), get_convex_hull_indices(tri))
                @test δ ≥ 0
            end
        end
    end
end

tri = example_triangulation()
tri = triangulate(get_points(tri))
vorn = voronoi(tri)
@test validate_tessellation(vorn)
bbox = DT.polygon_bounds(vorn, 0.1)
xmin, xmax, ymin, ymax = bbox
bbox = [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)]
bbox_order = [1, 2, 3, 4, 1]
c1 = DT.get_polygon_coordinates(vorn, 1, bbox, bbox_order)
c2 = DT.get_polygon_coordinates(vorn, 2, bbox, bbox_order)
c3 = DT.get_polygon_coordinates(vorn, 3, bbox, bbox_order)
c4 = DT.get_polygon_coordinates(vorn, 4, bbox, bbox_order)
c5 = DT.get_polygon_coordinates(vorn, 5, bbox, bbox_order)
c6 = DT.get_polygon_coordinates(vorn, 6, bbox, bbox_order)
c7 = DT.get_polygon_coordinates(vorn, 7, bbox, bbox_order)
@test DT.circular_equality(collect.(c1), collect.([(-1.5, 0.5)
    (0.16666666666666666, -1.1666666666666665)
    (1.0, 0.5)
    (1.0, 3.0)
    (-1.5, 0.5)]), ≈)
@test DT.circular_equality(collect.(c2), collect.([(3.5, 0.5)
    (0.5, -2.5)
    (0.5, -3.2000000001862645)
    (5.769999999552965, -1.7699999995529652)
    (3.5, 0.5)]), ≈)
@test DT.circular_equality(collect.(c3), collect.([(3.5, 0.5)
    (1.0, 0.5)
    (0.16666666666666666, -1.1666666666666665)
    (0.5, -2.5)
    (3.5, 0.5)]), ≈)
@test DT.circular_equality(collect.(c4), collect.([(1.5, 5.270000000484288)
    (-2.6999999997206032, 0.8999999999068677)
    (-1.5, 0.5)
    (1.0, 3.0)
    (1.5, 4.5)
    (1.5, 5.270000000484288)]), ≈)
@test DT.circular_equality(collect.(c5), collect.([(1.5, 4.5)
    (3.5, 0.5)
    (5.769999999552965, 2.769999999552965)
    (1.5, 5.270000000484288)
    (1.5, 4.5)]), ≈)
@test DT.circular_equality(collect.(c6), collect.([(-2.6999999997206032, 0.8999999999068677)
    (0.5, -3.2000000001862645)
    (0.5, -2.5)
    (0.16666666666666666, -1.1666666666666665)
    (-1.5, 0.5)
    (-2.6999999997206032, 0.8999999999068677)]), ≈)
@test DT.circular_equality(collect.(c7), collect.([(1.5, 4.5)
    (1.0, 3.0)
    (1.0, 0.5)
    (3.5, 0.5)
    (1.5, 4.5)]), ≈)

_vorn = voronoi(tri, true)



fig, ax, sc = triplot(tri)

lines!(ax, bbox[bbox_order], color=:blue)
fig
lines!(ax, [p, q], color=:red, linewidth=4)

fig

triplot(tri)

fig, ax, sc = triplot(tri)
voronoiplot!(ax, vorn)

j = 4
C = get_polygon(vorn, j)
F = Float64
coords = Vector{NTuple{2,F}}(undef, length(C) - 1)
i = 1
ghost_tri = DT.get_circumcenter_to_triangle(vorn, C[i])
u, v, _ = indices(ghost_tri) # w is the ghost vertex
p, q = get_generator(vorn, u, v)
px, py = getxy(p)
qx, qy = getxy(q)
m = (px + qx) / 2, (py + qy) / 2
is_first = DT.is_first_boundary_index(C, i)
prev_index = DT.previndex_circular(C, i)
r = get_polygon_point(vorn, C[prev_index])
r = getxy(r)
coords[i] = intersection = DT.intersection_of_ray_with_boundary(bbox, bbox_order, m, r)

i = 2
ghost_tri = DT.get_circumcenter_to_triangle(vorn, C[i])
u, v, _ = indices(ghost_tri) # w is the ghost vertex
p, q = get_generator(vorn, u, v)
px, py = getxy(p)
qx, qy = getxy(q)
m = (px + qx) / 2, (py + qy) / 2
is_first = DT.is_first_boundary_index(C, i)
next_index = DT.nextindex_circular(C, i)
r = get_polygon_point(vorn, C[next_index])
r = getxy(r)
if r == m # It's possible for the circumcenter to lie on the edge and exactly at the midpoint (e.g. [(0.0,1.0),(-1.0,2.0),(-2.0,-1.0)]). In this case, just rotate 
    mx, my = getxy(m)
    dx, dy = qx - mx, qy - my
    rotated_dx, rotated_dy = -dy, dx
    r = mx + rotated_dx, my + rotated_dy
    if DT.is_right(DT.point_position_relative_to_line(p, q, r))
        rotated_dx, rotated_dy = dy, -dx
        r = mx + rotated_dx, my + rotated_dy
    end
end
r = getxy(r)
intersection = DT.intersection_of_ray_with_boundary(bbox, bbox_order, m, r)



px, py = getxy(p)
qx, qy = getxy(q)
t1 = zero(px)
t2 = one(px)
points = bbox
boundary_nodes = bbox_order
δ1 = DT.distance_to_polygon(p, points, boundary_nodes)
δ2 = DT.distance_to_polygon(q, points, boundary_nodes)
while sign(δ2) == 1
    t2 *= 2
    r = px + t2 * (qx - px), py + t2 * (qy - py)
    δ2 = DT.distance_to_polygon(r, points, boundary_nodes)
end
t = (t1 + t2) / 2
r = px + t * (qx - px), py + t * (qy - py)
δ = DT.distance_to_polygon(r, points, boundary_nodes)



tri = example_triangulation()
tri = triangulate(get_points(tri), delete_ghosts=false)
lock_convex_hull!(tri)
vorn = voronoi(tri, false)
boundary_edges = (keys ∘ get_boundary_edge_map)(tri)
boundary_sites = Set{Int64}()
E = DT.edge_type(tri)
boundary_sites = Dict{Int64,E}()
for e in boundary_edges
    i, j = DT.edge_indices(e)
    p, q = get_point(tri, i, j)
    px, py = getxy(p)
    qx, qy = getxy(q)
    m = (px + qx) / 2, (py + qy) / 2
    incident_polygon = jump_and_march(vorn, m; k=i)
    boundary_sites[i] = e
    boundary_sites[j] = e
    boundary_sites[incident_polygon] = e
end
E = NTuple{2,Int64}
F = DT.number_type(vorn)
segment_intersections = NTuple{2,F}[]
intersected_edge_cache = Vector{E}(undef, 2)
boundary_site_additions = Dict{Int64,Set{Int64}}()
boundary_site_deletions = Dict{Int64,Set{Int64}}()
for (incident_polygon, e) in boundary_sites
    left_edge, right_edge, left_bnd, right_bnd = DT.get_neighbouring_edges(vorn, e)
    polygon = DT.get_polygon(vorn, incident_polygon)
    nedges = DT.num_boundary_edges(polygon)
    num_intersections = 0
    for ℓ in 1:nedges
        num_intersections == 2 && break
        u = DT.get_boundary_nodes(polygon, ℓ)
        v = DT.get_boundary_nodes(polygon, ℓ + 1)
        if DT.is_boundary_index(u) && DT.is_boundary_index(v)
            continue
        elseif DT.is_boundary_index(u) && !DT.is_boundary_index(v)
            num_intersections = DT.process_ray_intersection!(vorn, u, v, incident_polygon, num_intersections, intersected_edge_cache, segment_intersections, boundary_site_additions)
        elseif !DT.is_boundary_index(u) && DT.is_boundary_index(v)
            num_intersections = DT.process_ray_intersection!(vorn, v, u, incident_polygon, num_intersections, intersected_edge_cache, segment_intersections, boundary_site_additions)
        else
            for e in (e, left_edge, right_edge)
                _num_intersections = DT.process_segment_intersection!(vorn, u, v, e, incident_polygon, num_intersections, intersected_edge_cache, segment_intersections, boundary_site_additions, boundary_site_deletions)
                if _num_intersections > num_intersections
                    num_intersections = _num_intersections
                    # break < -- Actually, don't break because a single line could go past multiple parts of the boundary (e.g. the fixed_shewchuk_example_constrained example).
                    num_intersections == 2 && break
                end
            end
        end
    end
    interior_intersection, index, adjacent_incident_polygon, shared_vertex = DT.segment_intersection_type(e, intersected_edge_cache)
    if !interior_intersection
        other_e = intersected_edge_cache[index]
        m = get_generator(vorn, shared_vertex)
        DT.add_segment_intersection!(segment_intersections, boundary_site_additions, m, incident_polygon)
    end
end
n = num_polygon_vertices(vorn)
for p in each_point(segment_intersections)
    DT.push_polygon_point!(vorn, p)
end
boundary_sites = keys(boundary_sites)

polygon = 4
polygon_vertices = get_polygon(vorn, polygon)
pop!(polygon_vertices)
filter!(!DT.is_boundary_index, polygon_vertices)




polygon = 5 
polygon_vertices = get_polygon(vorn, polygon)
pop!(polygon_vertices)
filter!(!DT.is_boundary_index, polygon_vertices)
polygon_points = [get_polygon_point(vorn, polygon_vertices...)...]








polygon = 5
polygon_vertices = get_polygon(vorn, polygon)
unique!(polygon_vertices)
filter!(!is_boundary_index, polygon_vertices)
for new_vert in boundary_site_additions[polygon]
    push!(polygon_vertices, n + new_vert)
end
cx, cy = DT.arithmetic_average_unsorted(vorn, polygon_vertices)
θ = zeros(F, length(polygon_vertices))
for (j, i) in pairs(polygon_vertices)
    p = get_polygon_point(vorn, i)
    px, py = getxy(p)
    θ[j] = atan(py - cy, px - cx)
end
idx = sortperm(θ)
permute!(polygon_vertices, idx)
push!(polygon_vertices, polygon_vertices[begin])
lines!(ax, [get_polygon_point(vorn, polygon_vertices...)...], color=:red)

polygon = 4
polygon_vertices = get_polygon(vorn, polygon)
pop!(polygon_vertices)
filter!(!is_boundary_index, polygon_vertices)
for new_vert in boundary_site_additions[polygon]
    push!(polygon_vertices, n + new_vert)
end
if haskey(boundary_site_deletions, polygon)
    deleteat!(polygon_vertices, findfirst(isequal(boundary_site_deletions[polygon]), polygon_vertices))
end
cx, cy = DT.arithmetic_average_unsorted(vorn, polygon_vertices)
θ = zeros(F, length(polygon_vertices))
for (j, i) in pairs(polygon_vertices)
    p = get_polygon_point(vorn, i)
    px, py = getxy(p)
    θ[j] = atan(py - cy, px - cx)
end
idx = sortperm(θ)
permute!(polygon_vertices, idx)
push!(polygon_vertices, polygon_vertices[begin])
lines!(ax, [get_polygon_point(vorn, polygon_vertices...)...], color=:red)

polygon = 6
polygon_vertices = get_polygon(vorn, polygon)
pop!(polygon_vertices)
filter!(!is_boundary_index, polygon_vertices)
for new_vert in boundary_site_additions[polygon]
    push!(polygon_vertices, n + new_vert)
end
if haskey(boundary_site_deletions, polygon)
    deleteat!(polygon_vertices, findfirst(isequal(boundary_site_deletions[polygon]), polygon_vertices))
end
cx, cy = DT.arithmetic_average_unsorted(vorn, polygon_vertices)
θ = zeros(F, length(polygon_vertices))
for (j, i) in pairs(polygon_vertices)
    p = get_polygon_point(vorn, i)
    px, py = getxy(p)
    θ[j] = atan(py - cy, px - cx)
end
idx = sortperm(θ)
permute!(polygon_vertices, idx)
push!(polygon_vertices, polygon_vertices[begin])
lines!(ax, [get_polygon_point(vorn, polygon_vertices...)...], color=:red)

polygon = 7
polygon_vertices = get_polygon(vorn, polygon)
pop!(polygon_vertices)
filter!(!is_boundary_index, polygon_vertices)
for new_vert in boundary_site_additions[polygon]
    push!(polygon_vertices, n + new_vert)
end
if haskey(boundary_site_deletions, polygon)
    deleteat!(polygon_vertices, findfirst(isequal(boundary_site_deletions[polygon]), polygon_vertices))
end
cx, cy = DT.arithmetic_average_unsorted(vorn, polygon_vertices)
θ = zeros(F, length(polygon_vertices))
for (j, i) in pairs(polygon_vertices)
    p = get_polygon_point(vorn, i)
    px, py = getxy(p)
    θ[j] = atan(py - cy, px - cx)
end
idx = sortperm(θ)
permute!(polygon_vertices, idx)
push!(polygon_vertices, polygon_vertices[begin])
lines!(ax, [get_polygon_point(vorn, polygon_vertices...)...], color=:red)

polygon = 2
polygon_vertices = get_polygon(vorn, polygon)
pop!(polygon_vertices)
filter!(!is_boundary_index, polygon_vertices)
for new_vert in boundary_site_additions[polygon]
    push!(polygon_vertices, n + new_vert)
end
if haskey(boundary_site_deletions, polygon)
    deleteat!(polygon_vertices, findfirst(isequal(boundary_site_deletions[polygon]), polygon_vertices))
end
cx, cy = DT.arithmetic_average_unsorted(vorn, polygon_vertices)
θ = zeros(F, length(polygon_vertices))
for (j, i) in pairs(polygon_vertices)
    p = get_polygon_point(vorn, i)
    px, py = getxy(p)
    θ[j] = atan(py - cy, px - cx)
end
idx = sortperm(θ)
permute!(polygon_vertices, idx)
push!(polygon_vertices, polygon_vertices[begin])
lines!(ax, [get_polygon_point(vorn, polygon_vertices...)...], color=:cyan)

polygon = 3
polygon_vertices = get_polygon(vorn, polygon)
pop!(polygon_vertices)
filter!(!is_boundary_index, polygon_vertices)
for new_vert in boundary_site_additions[polygon]
    push!(polygon_vertices, n + new_vert)
end
if haskey(boundary_site_deletions, polygon)
    deleteat!(polygon_vertices, findfirst(isequal(boundary_site_deletions[polygon]), polygon_vertices))
end
cx, cy = DT.arithmetic_average_unsorted(vorn, polygon_vertices)
θ = zeros(F, length(polygon_vertices))
for (j, i) in pairs(polygon_vertices)
    p = get_polygon_point(vorn, i)
    px, py = getxy(p)
    θ[j] = atan(py - cy, px - cx)
end
idx = sortperm(θ)
permute!(polygon_vertices, idx)
push!(polygon_vertices, polygon_vertices[begin])
lines!(ax, [get_polygon_point(vorn, polygon_vertices...)...], color=:red)

polygon = 1
polygon_vertices = get_polygon(vorn, polygon)
pop!(polygon_vertices)
filter!(!is_boundary_index, polygon_vertices)
for new_vert in boundary_site_additions[polygon]
    push!(polygon_vertices, n + new_vert)
end
if haskey(boundary_site_deletions, polygon)
    deleteat!(polygon_vertices, findfirst(isequal(boundary_site_deletions[polygon]), polygon_vertices))
end
cx, cy = DT.arithmetic_average_unsorted(vorn, polygon_vertices)
θ = zeros(F, length(polygon_vertices))
for (j, i) in pairs(polygon_vertices)
    p = get_polygon_point(vorn, i)
    px, py = getxy(p)
    θ[j] = atan(py - cy, px - cx)
end
idx = sortperm(θ)
permute!(polygon_vertices, idx)
push!(polygon_vertices, polygon_vertices[begin])
lines!(ax, [get_polygon_point(vorn, polygon_vertices...)...], color=:blue, linewidth=6)
