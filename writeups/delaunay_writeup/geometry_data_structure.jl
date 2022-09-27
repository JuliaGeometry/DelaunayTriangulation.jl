using Meshes, MeshViz
import CairoMakie: Figure, Axis, poly!, xlims!, ylims!, scatter!, text!, @L_str, save

## Define the Delaunay triangulation using basic objects
a = Point(5, 5);    b = Point(4.5, 2.5);    c = Point(2.5, 1.5); 
d = Point(3, 3.5);  e = Point(0, 2);        f = Point(1, 5)
g = Point(1, 3);    h = Point(4, -1);       i = Point(-1, 4)
pts = [a, b, c, d, e, f, g, h, i]
tris = [7 3 4; 7 4 6; 7 6 9; 7 9 5
        7 5 3; 2 8 1; 2 3 8; 2 4 3
        2 1 4; 3 5 8; 4 1 6] # Make sure it's counter-clockwise
## The triangles need to be connected using connect 
tris = [connect(Tuple(tri), Triangle) for tri in eachrow(tris)]
## Using Meshes.jl, we can construct a SimpleMesh
mesh = SimpleMesh(pts, tris)

## Using MeshViz.jl, we can easily visualise the mesh
fig, ax, p = viz(mesh, showfacets=true)
## Harder to customise, though, so here's the original Makie version 
fig = Figure(fontsize=32)
ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y",
    xticks=(-2:2:6, [L"%$s" for s in -2:2:6]),
    yticks=(-2:2:6, [L"%$s" for s in -2:2:6]))
xlims!(ax, -2, 6); ylims!(ax, -2, 6)
poly!(ax, coordinates.(pts), 
    reshape(reinterpret(Int64, indices.(tris)), (3, :))', # https://stackoverflow.com/a/68550861
    color=(:white, 0), strokewidth=2)
scatter!(ax, coordinates.(pts))
text!(ax, [L"p_%$i" for i in 0:8],
    position=[pts[1] - Point(-0.05, -0.05), pts[2] - Point(0.25, -0.2), pts[3] - Point(0.2, 0.5),
        pts[4] - Point(0.0, -0.2), pts[5] - Point(0.4, 0.4), pts[6] - Point(0.0, -0.1),
        pts[7] - Point(-0.1, -0.1), pts[8] - Point(-0.1, 0.0), pts[9] - Point(0.0, -0.1)],
    space=:data, textsize=32)
save("$(@__DIR__)/figures/example_delaunay.pdf", fig)

topo = convert(HalfEdgeTopology, topology(mesh))
edge_to_face = Coboundary{1,2}(topo)    # edges  ⟹ triangles
edge_to_points = Boundary{1,0}(topo)    # points ⟹ vertices 
for k in 1:nfacets(topo)
    simplices = edge_to_face(k) # Triangles incident to the kth edge
    points = edge_to_points(k)
end

for edge in facets(topo)
    println(edge)
end