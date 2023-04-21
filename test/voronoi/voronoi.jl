using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using ColorSchemes
using DataStructures
using LinearAlgebra

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


A = (-1.0, 7.0)
B = (4.0, 4.0)
C = (-2.0, -1.0)
D = (-1.0, 3.0)
E = (3.0, -1.0)
F = (1.0, 4.0)
G = (-3.0, 5.0)
pts = [A, B, C, D, E, F, G]
tri = triangulate(pts; delete_ghosts=false, randomise=false)
lock_convex_hull!(tri)
vorn = voronoi(tri, true)
voronoiplot(vorn)

tri = DT.get_triangulation(vorn)
boundary_edges = (keys ∘ get_boundary_edge_map)(tri)
I = DT.integer_type(tri)
boundary_sites = Set{I}()
E = DT.edge_type(tri)
boundary_sites = Dict{I,E}()
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
F = DT.number_type(vorn)
segment_intersections = NTuple{2,F}[]
intersected_edge_cache = Vector{E}(undef, 2)
boundary_site_additions = Dict{I,Set{I}}()
boundary_site_deletions = Dict{I,I}()
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
        elseif !is_boundary_index(u) && DT.is_boundary_index(v)
            num_intersections = DT.process_ray_intersection!(vorn, v, u, incident_polygon, num_intersections, intersected_edge_cache, segment_intersections, boundary_site_additions)
        else
            for e in (e, left_edge, right_edge)
                _num_intersections = DT.process_segment_intersection!(vorn, u, v, e, incident_polygon, num_intersections, intersected_edge_cache, segment_intersections, boundary_site_additions, boundary_site_deletions)
                if _num_intersections > num_intersections
                    num_intersections = _num_intersections
                    break
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
