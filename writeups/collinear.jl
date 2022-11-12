include("functions.jl")

a = 0.0
b = 10.0
c = 0.0
d = 10.0
nx = 11
ny = 11
function DT._get_point(pts::AbstractMatrix, i)
    return @view pts[:, i]
end
function DT._eachindex(pts::AbstractMatrix)
    return axes(pts, 2)
end
T, adj, adj2v, DG, pts = triangulate_structured(a, b, c, d, nx, ny)
DT.compute_centroid!(pts)

fig = Figure()
ax = Axis(fig[1, 1])
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0), plot_ghost_edges=true, DG=DG, recompute_centroid=false)
xlims!(ax, a - 1 / 2, b + 1 / 2)
ylims!(ax, c - 1 / 2, d + 1 / 2)
fig

Random.seed!(299756756791)
T, adj, adj2v, DG, pts = triangulate_structured(a, b, c, d, nx, ny)

IntegerType = I = Int64
EdgeType = E = NTuple{2,IntegerType}
TriangleType = V = NTuple{3,IntegerType}
EdgesType = Es = Set{EdgeType}
TrianglesType = Ts = Set{TriangleType}
randomise = false
trim = true
try_last_inserted_point = true
skip_pts = Set{Int64}()

pt_order = randomise ? shuffle(DT._eachindex(pts)) : collect(DT._eachindex(pts))
setdiff!(pt_order, skip_pts)
T = Ts()
adj = DT.Adjacent{I,E}()
adj2v = DT.Adjacent2Vertex{I,Es,E}()
DG = DT.DelaunayGraph{I}()
initial_triangle = DT.construct_positively_oriented_triangle(TriangleType, I(pt_order[begin]), I(pt_order[begin+1]), I(pt_order[begin+2]), pts)
while DT.isoriented(initial_triangle, pts) == 0 # We cannot start with a degenerate triangle. We need to shift the insertion order
    circshift!(pt_order, 1)
    initial_triangle = DT.construct_positively_oriented_triangle(TriangleType, I(pt_order[begin]), I(pt_order[begin+1]), I(pt_order[begin+2]), pts)
end
u, v, w = indices(initial_triangle)
DT.add_triangle!(u, v, w, T, adj, adj2v, DG; update_ghost_edges=true)
DT.compute_centroid!((
    get_point(pts, pt_order[begin]),
    get_point(pts, pt_order[begin+1]),
    get_point(pts, pt_order[begin+2])
))
for (num_points, new_point) in enumerate(@view pt_order[(begin+3):end])
    if num_points < 119
        last_inserted_point_number = num_points + 3 - 1 # + 3 for the first three points already inserted
        last_inserted_point = pt_order[last_inserted_point_number]
        pt_idx = @view pt_order[begin:last_inserted_point_number] # + 3 for the first three points already inserted
        m = ceil(Int64, length(pt_idx)^(1 / 3))
        initial_search_point = I(DT.select_initial_point(pts, new_point; pt_idx, m, try_points=try_last_inserted_point ? last_inserted_point : ()))
        DT.add_point_bowyer!(T, adj, adj2v, DG, pts, new_point;
            pt_idx, m, initial_search_point)
        DT.update_centroid_after_new_point!(pts, new_point)
    else
        break
    end
end
check_ghost = true
num_points = 116
new_point = pt_order[(begin+3):end][num_points]
r = new_point
r = I(r)

boundary_neighbours = I[]
for u in DT.get_neighbour(DG, r)
    DT.is_boundary_point(u, adj, DG) && push!(boundary_neighbours, u)
end
!DT.is_boundary_edge(boundary_neighbours[1], boundary_neighbours[2], adj) && reverse!(boundary_neighbours)
i, j = boundary_neighbours
DT.delete_triangle!(j, i, r, T, adj, adj2v, DG; protect_boundary=true, update_ghost_edges=false)
DT.add_triangle!(r, j, I(DT.BoundaryIndex), T, adj, adj2v, DG; update_ghost_edges=false)
DT.add_triangle!(i, r, I(DT.BoundaryIndex), T, adj, adj2v, DG; update_ghost_edges=false)

point_neighbours = DT.get_neighbour(DG, r)
DT.is_boundary_point.(point_neighbours, Ref(adj), Ref(DG))


last_inserted_point_number = num_points + 3 - 1 # + 3 for the first three points already inserted
last_inserted_point = pt_order[last_inserted_point_number]
pt_idx = @view pt_order[begin:last_inserted_point_number]
m = ceil(Int64, length(pt_idx)^(1 / 3))
initial_search_point = I(DT.select_initial_point(pts, new_point; pt_idx, m, try_points=try_last_inserted_point ? last_inserted_point : ()))
k = initial_search_point
r = new_point
r = I(r)
tri_type = DT.triangle_type(Ts)
V = jump_and_march(r, adj, adj2v, DG, pts; k=initial_search_point, TriangleType=tri_type)
flag = DT.isintriangle(V, pts, r)
i, j, k = indices(V)
ℓ = DT.get_edge(adj, j, i) # (j, i, ℓ) is the triangle on the other side of the edge (i, j) from r 
DT.delete_triangle!(i, j, k, T, adj, adj2v, DG; protect_boundary=false, update_ghost_edges=true)
incirc = DT.isincircle(pts, r, i, j, ℓ)
pᵣ, pᵢ, pⱼ = get_point(pts, r, i, j)
points_are_collinear = DT.orient(pᵣ, pᵢ, pⱼ) == 0
DT.add_triangle!(i, r, k, T, adj, adj2v, DG, update_ghost_edges=true)
DT.add_triangle!(r, j, k, T, adj, adj2v, DG, update_ghost_edges=true)



fig = Figure()
ax = Axis(fig[1, 1])
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0), plot_ghost_edges=true, DG=DG, recompute_centroid=false)
xlims!(ax, a - 1 / 2, b + 1 / 2)
ylims!(ax, c - 1 / 2, d + 1 / 2)
fig


color = rand(size(pts, 2))
fig, ax, pl = poly(pts, [collect(T)[i][j] for i in 1:length(T), j in 1:3], color=color, strokewidth=2)


T, adj, adj2v, DG, pts = triangulate_structured(a, b, c, d, nx, ny)


Random.seed!(668591)
IntegerType = I = Int64
EdgeType = E = NTuple{2,IntegerType}
TriangleType = V = NTuple{3,IntegerType}
EdgesType = Es = Set{EdgeType}
TrianglesType = Ts = Set{TriangleType}
randomise = true
trim = true
try_last_inserted_point = true
skip_pts = Set{Int64}()
pt_order = randomise ? shuffle(DT._eachindex(pts)) : collect(DT._eachindex(pts))
setdiff!(pt_order, skip_pts)
T = Ts()
adj = DT.Adjacent{I,E}()
adj2v = DT.Adjacent2Vertex{I,Es,E}()
DG = DT.DelaunayGraph{I}()
initial_triangle = DT.construct_positively_oriented_triangle(TriangleType, I(pt_order[begin]), I(pt_order[begin+1]), I(pt_order[begin+2]), pts)
while DT.isoriented(initial_triangle, pts) == 0 # We cannot start with a degenerate triangle. We need to shift the insertion order
    circshift!(pt_order, 1)
    initial_triangle = DT.construct_positively_oriented_triangle(TriangleType, I(pt_order[begin]), I(pt_order[begin+1]), I(pt_order[begin+2]), pts)
end
u, v, w = indices(initial_triangle)
DT.add_triangle!(u, v, w, T, adj, adj2v, DG; update_ghost_edges=true)
DT.compute_centroid!((
    get_point(pts, pt_order[begin]),
    get_point(pts, pt_order[begin+1]),
    get_point(pts, pt_order[begin+2])
))
for (num_points, new_point) in enumerate(@view pt_order[(begin+3):end])
    @show num_points
    if num_points < 3
        last_inserted_point_number = num_points + 3 - 1 # + 3 for the first three points already inserted
        last_inserted_point = pt_order[last_inserted_point_number]
        pt_idx = @view pt_order[begin:last_inserted_point_number] # + 3 for the first three points already inserted
        m = ceil(Int64, length(pt_idx)^(1 / 3))
        initial_search_point = I(DT.select_initial_point(pts, new_point; pt_idx, m, try_points=try_last_inserted_point ? last_inserted_point : ()))
        DT.add_point_bowyer!(T, adj, adj2v, DG, pts, new_point;
            pt_idx, m, initial_search_point)
        DT.update_centroid_after_new_point!(pts, new_point)
        DT.validate_triangulation(T, adj, adj2v, DG, pts)
    else
        break
    end
end
num_points = 2
new_point = pt_order[(begin+3):end][num_points]
last_inserted_point_number = num_points + 3 - 1 # + 3 for the first three points already inserted
last_inserted_point = pt_order[last_inserted_point_number]
pt_idx = @view pt_order[begin:last_inserted_point_number]
m = ceil(Int64, length(pt_idx)^(1 / 3))
initial_search_point = I(DT.select_initial_point(pts, new_point; pt_idx, m, try_points=try_last_inserted_point ? last_inserted_point : ()))
r = new_point
r = I(r)
k = initial_search_point
tri_type = DT.triangle_type(Ts)
i, j = DT.straight_line_search_ghost_triangles(q, adj, k, pts)


V = jump_and_march(r, adj, adj2v, DG, pts; k=initial_search_point, TriangleType=tri_type)
flag = DT.isintriangle(V, pts, r)
i, j, k = indices(V)
ℓ₁ = DT.get_edge(adj, j, i)
ℓ₂ = DT.get_edge(adj, k, j)
ℓ₃ = DT.get_edge(adj, i, k)
DT.delete_triangle!(i, j, k, T, adj, adj2v, DG; protect_boundary=true, update_ghost_edges=false)
DT.dig_cavity!(r, i, j, ℓ₁, T, adj, adj2v, DG, pts)
DT.dig_cavity!(r, j, k, ℓ₂, T, adj, adj2v, DG, pts)
DT.dig_cavity!(r, k, i, ℓ₃, T, adj, adj2v, DG, pts)

flag == 0 && (DT.is_boundary_triangle(V, adj) || DT.is_ghost_triangle(V)) && !DT.is_boundary_point(r, adj, DG)


p = get_point(pts, k)
i = DT.get_edge(adj, k, I(DT.BoundaryIndex)) # Note that this is at the left
pᵢ = get_point(pts, i)
o1 = DT.orient(p, q, pᵢ) # Is pᵢ to the left of pq?
j = DT.get_edge(adj, k, i)
pⱼ = get_point(pts, j)
o2 = DT.orient(p, q, pⱼ) # Is pⱼ to the left of pq?
o1, i, pᵢ = o2, j, pⱼ # Step onto the next triangle

j = DT.get_edge(adj, k, i)
pⱼ = get_point(pts, j)
o2 = DT.orient(p, q, pⱼ) # Is pⱼ to the left of pq?
o1, i, pᵢ = o2, j, pⱼ # Step onto the next triangle

j = DT.get_edge(adj, k, i)
pⱼ = get_point(pts, j)
o2 = DT.orient(p, q, pⱼ) # Is pⱼ to the left of pq?
o2 == 0 && !DT.is_boundary_point(j, adj, DG)
DT.closest(p, q, pⱼ) == -1


o1, i, pᵢ = o2, j, pⱼ # Step onto the next triangle


p, i, j, pᵢ, pⱼ = DT.select_initial_triangle(q, adj, adj2v, DG, k, pts)

k = DT.get_edge(adj, i, j)
pₖ = get_point(pts, k)
i = k
pᵢ = pₖ

k = DT.get_edge(adj, i, j)
pₖ = get_point(pts, k)



V = jump_and_march(r, adj, adj2v, DG, pts; k=initial_search_point, TriangleType=tri_type)


k = initial_search_point
q = get_point(pts, r)
p = get_point(pts, k)
i, j = rand(DT.get_edge(adj2v, k))
pᵢ = get_point(pts, i)
pⱼ = get_point(pts, j)

o1 = DT.orient(p, q, pⱼ)
i = j
pᵢ = pⱼ
j = DT.get_edge(adj, k, j)
pⱼ = get_point(pts, j)

o1 = DT.orient(p, q, pⱼ)
i = j
pᵢ = pⱼ
j = DT.get_edge(adj, k, j)
pⱼ = get_point(pts, j)

o1 = DT.orient(p, q, pⱼ)
i = j
pᵢ = pⱼ
j = DT.get_edge(adj, k, j)
pⱼ = get_point(pts, j)

o1 = DT.orient(p, q, pⱼ)
i = j
pᵢ = pⱼ
j = DT.get_edge(adj, k, j)
pⱼ = get_point(pts, j)

q = get_point(pts, r)
p = get_point(pts, k)
i = DT.get_edge(adj, k, I(DT.BoundaryIndex)) # Note that this is at the left
pᵢ = get_point(pts, i)
o1 = DT.orient(p, q, pᵢ) # Is pᵢ to the left of pq?
DT.sameside(pᵢ, p, q) == -1
j = DT.get_edge(adj, i, k)
pⱼ = get_point(pts, j)
DT.orient(pⱼ, pᵢ, p) ≥ 0 && DT.orient(pᵢ, pₖ, p) ≥ 0 && DT.orient(pₖ, pᵢ, p) ≥ 0

j = DT.get_edge(adj, k, i)
pⱼ = get_point(pts, j)
o2 = DT.orient(p, q, pⱼ) # Is pⱼ to the left of pq?
o1, i, pᵢ = o2, j, pⱼ # Step onto the next triangle
j = DT.get_edge(adj, k, i)
pⱼ = get_point(pts, j)
o2 = DT.orient(p, q, pⱼ) # Is pⱼ to the left of pq?

DT.sameside(pⱼ, p, q) == -1

i, j, intersects_edge, inside_triangle, k = DT.check_interior_edge_intersections(q, adj, DG, k, pts)
p, pᵢ, pⱼ = get_point(pts, k, i, j)

k = DT.get_edge(adj, i, j)

V = jump_and_march(r, adj, adj2v, DG, pts; k=initial_search_point, TriangleType=tri_type)


a = 0.0
b = 10.0
c = 0.0
d = 10.0
nx = 11
ny = 11
function DT._get_point(pts::AbstractMatrix, i)
    return @view pts[:, i]
end
function DT._eachindex(pts::AbstractMatrix)
    return axes(pts, 2)
end
_T, _adj, _adj2v, _DG, pts = triangulate_structured(a, b, c, d, nx, ny)
Random.seed!(20288888)
T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
@test DT.validate_triangulation(T, adj, adj2v, DG, pts)