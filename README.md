# DelaunayTriangulation

[![DOI](https://zenodo.org/badge/540660309.svg)](https://zenodo.org/badge/latestdoi/540660309)

- [DelaunayTriangulation](#delaunaytriangulation)
- [Getting started](#getting-started)
  - [de Berg's method](#de-bergs-method)
  - [Bowyer-Watson algorithm](#bowyer-watson-algorithm)
  - [Structured triangulation](#structured-triangulation)
  - [Gmsh](#gmsh)
- [Customisation](#customisation)
  - [Triangle interface](#triangle-interface)
  - [Edge interface](#edge-interface)
  - [Point interface](#point-interface)
  - [Example](#example)
- [A note on Voronoi tessellations](#a-note-on-voronoi-tessellations)

This is a package for creating (unconstrained) two-dimensional Delaunay triangulations. The great package [ExactPredicates.jl](https://github.com/lairez/ExactPredicates.jl) is used for all geometric predicates. There also exist routines for building the convex hull (from the triangulation), and currently commented-out Voronoi tessellation code (see some discussion at the end). Much of the work in this package is derived from the book *Delaunay Mesh Generation* by Cheng, Dey, and Shewchuk (2013). Point location is implemented using the jump-and-march algorithm of Mücke, Saias, and Zhu (1999); see `jump_and_march`. I hope to eventually build up to constrained Delaunay triangulations, weighted Delaunay triangulations, etc., when I eventually find the time. Mesh refinement is also a priority.

The package has two algorithms for computing Delauanay triangulations, namely de Berg's method from de Berg et al. (1999), and the Bowyer-Watson algorithm as presented by Cheng, Dey, and Shewchuk (2013). de Berg's method is the slowest of the two, and has less features (and is somehow more problematic as evident by the three existing issues for de Berg's method), so we recommend the Bowyer-Watson algorithm in the function `triangulate_bowyer`. Keep reading for more examples, and see the tests for more features. You can also see the `writeups` folder for some (very) rough code that I used while testing some of this code -- I make no promises that all the code there still works.

I provide a PDF in the main directory, called `DelaunayTriangulation.pdf`, that outlines some of my working for the algorithms used in this package. I have tried to keep up with it, but feel free to ask about it or raise any issues that are in the document. I plan on rewriting it once I eventually add constrained/weighted Delaunay triangulations in the future.

Feel free to use the issues tab for any questions / feedback / etc (or email me at vandenh2@qut.edu.au).

# Getting started

## de Berg's method

It is easy to construct a triangulation. Here is a simple example, where we use the method outlined in the book by de Berg et al. (1999) to construct the triangulation; this triangulation adds points one at a time, constructing the Delaunay triangulation $\mathcal D\mathcal T(\mathcal P_n)$ from $\mathcal D\mathcal T(\mathcal P_{n-1})$, $\mathcal P_n = \mathcal P_{n-1} \cup {p_n}$ using edge flips. Using `triangulate_berg`, this triangulation can be computed:

```julia
using DelaunayTriangulation
p1 = [0.0, 1.0]
p2 = [3.0, -1.0]
p3 = [2.0, 0.0]
p4 = [-1.0, 2.0]
p5 = [4.0, 2.0]
p6 = [-2.0, -1.0]
pts = [p1, p2, p3, p4, p5, p6]
tri, HG = triangulate_berg(pts)
T, adj, adj2v, DG = tri # tri is a Triangulation struct, see ?Triangulation
```
This function `triangulate_berg` creates the triangulation for the given `pts`, returning:
- `T`: This is the set of triangles, of default type `Set{NTuple{3, Int64}}`. By default, triangles are represented as tuples of indices `(i, j, k)` so that `(pts[i], pts[j], pts[k])` are the points representing the triangle. All triangles `T` are positively oriented.
- `adj`: This is the `Adjacent` map, mapping edges `(i, j)` to a vertex `k` such that `(i, j, k)` is a positively oriented (counter-clockwise) triangle that exists in `T`. If no such triangle exists, `-4` is returned (this constant is defined in `DelaunayTriangulation.DefaultAdjacentValue`). You can extract this vertex using `get_edge(adj, i, j)`. By default, edges are represented as `NTuple{2, Int64}`. These empty keys (those that map to `-4`) can be removed using `DelaunayTriangulation.clear_empty_keys!(adj)`.
- `adj2v`: This is the `Adjacent2Vertex` map, mapping vertices `i` to a set of edges `(j, k)` such that `(i, j, k)` is a positively oriented triangle in `T`. In particular, the set of all triangles with `i` as a vertex is `[(i, j, k) for (j, k) in get_edge(adj2v, i)]`.
- `DG`: This is the `DelaunayGraph` object, mapping vertices `i` to the set of all neighbouring vertices. You can access all these neighbours using `get_neighbour(DG, i)`.
- `HG`: Unique to `triangulate_berg`, this is a `HistoryGraph` object, primarily used for point location. For example, if a triangulation is constructed for points `pts[1:(n-1)]` and a new point `pts[n]` is to be added, the triangle that `pts[n]` is inside of can be found using `locate_triangle(HG, pts, n)`.

Some key points to note are that if an edge `(i, j)` is a boundary edge, defined to be an edge `(i, j)` with no vertex `k` such that `(i, j, k)` is a positively oriented triangle in `T`, then `get_edge(adj, i, j) = 0`, where the constant `0` is determined by `DelaunayTriangulation.BoundaryIndex`. The entire set of boundary edges can be determined using `get_edge(adj2v, 0)`, and the entire set of boundary nodes can be determined using `get_neighbour(DG, 0)`. These facts make it easy to extract the convex hull:

```julia
ch = convex_hull(DG, pts)
```

The returned value `ch` is a `Vector{Int64}` and is such that the boundary nodes are listed counter-clockwise.

A key limitation of this algorithm is that points cannot be removed (although this is not yet implemented for any method -- this is planned).

## Bowyer-Watson algorithm 

The Bowyer-Watson algorithm is different to de Berg's method, although points are still added one at a time. Instead of flipping edges to reach a Delaunay triangulation at each step, triangles are instead evacuated and then the resulting polygonal cavity is re-filled to reach a Delaunay triangulation. Our implementation follows that of Cheng, Dey, and Shewchuk (2013). The function for this is `triangulate_bowyer`:
```julia
using StaticArraysCore
p1 = SVector{2, Float64}([0.0, 1.0])
p2 = SVector{2, Float64}([3.0, -1.0])
p3 = SVector{2, Float64}([2.0, 0.0])
p4 = SVector{2, Float64}([-1.0, 2.0])
p5 = SVector{2, Float64}([4.0, 2.0])
p6 = SVector{2, Float64}([-2.0, -1.0])
pts = [p1, p2, p3, p4, p5, p6]
T, adj, adj2v, DG = triangulate_bowyer(pts)
```
(Points don't have to be represented as vectors, they could be `Tuple`s or `SVector`s, or any type as long as you implement the interface correctly; see the Customisation section.) Points are added in random order, but you can use the order defined by the order in `pts` by setting the keyword `randomise` to `false` in `triangulate_bowyer`.

The definitions of `T`, `adj`, `adj2v`, and `DG` are the same as above. One nice feature of the Bowyer-Watson algorithm as implemented is the use of ghost triangles, i.e. a ghost triangle adjoined to a boundary edge with a vertex out at infinity. This makes it easy to locate points in the boundary, defining a unique ghost triangle that they reside in, and makes it easy to traverse the boundary. You can keep these ghost triangles by setting `trim=false`:

```julia
T, adj, adj2v, DG = triangulate_bowyer(pts; trim=false)
```

Alternatively, you can add it after the fact as follows:

```julia
DelaunayTriangulation.add_ghost_triangles!(T, adj, adj2v, DG)
```

(Use `DelaunayTriangulation.remove_ghost_triangles!(T, adj, adj2v, DG)` to remove them.) The vertex at infinity is represented by `0`.

Visualising these triangulations is easy:

```julia
pts = 15rand(2, 500)
T, adj, adj2v, DG = triangulate_bowyer(pts)
ch = convex_hull(DG, pts)
using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1])
poly!(ax, pts, [collect(T)[i][j] for i in 1:length(T), j in 1:3], color = (:white, 0.0), strokewidth = 2)
lines!(ax, hcat(pts[:, ch], pts[:, ch[1]]), color = :red, linewidth = 4)
```

![A triangulation](https://github.com/DanielVandH/DelaunayTriangulation.jl/blob/main/test_fig.png?raw=true)

## Structured triangulation

It may sometimes be convenient to have a structured triangulation for a rectangular geometry. A structured triangulation is simply one where all the triangles have the same shape. Here is an example:
```julia 
a = 0.0
b = 10.0
c = 0.0
d = 10.0
nx = 20
ny = 10
T, adj, adj2v, DG, pts = triangulate_structured(a, b, c, d, nx, ny)
fig = Figure()
ax = Axis(fig[1, 1])
poly!(ax, pts, [collect(T)[i][j] for i in 1:length(T), j in 1:3], color = (:white, 0.0), strokewidth = 2)
```

![A structured triangulation](https://github.com/DanielVandH/DelaunayTriangulation.jl/blob/main/test_fig2.png?raw=true)

## Gmsh 

Support is also added for a simple mesh generator with Gmsh (see https://gmsh.info/), tested up to v4.9.4 on Windows 64. The function for this is `generate_mesh`, and accepts inputs of boundary points that are in counter-clockwise order. This is especially useful for e.g. finite volume codes. Currently I only have code working for simply connected domains - it would be nice to have an alternative, but this is the best I can do with the time I have (the alternative would require me to think a lot more about ghost nodes, boundary edges, etc. when the domain has holes, and the impact this would have on the existing code and existing data structures).

Let me give an example. In my directory, I have downloaded `gmsh` and saved it as `gmsh-4.9.4-Windows64`, so I define
```julia
GMSH_PATH = "./gmsh-4.9.4-Windows64/gmsh.exe"
```
Then, for the meshing algorithm, I want to use a frontal Delaunay algorithm, so I set 
```julia
mesh_algo = 6
```
(These numbers are defined in the [Gmsh documentation](https://gmsh.info/doc/texinfo/gmsh.html).) Next, let's define our geometry. Let us start with a simple geometry, the circle:
```julia
θ = LinRange(0, 2π, 250)
x = cos.(θ)
y = sin.(θ)
T, adj, adj2v, DG, pts, BN = generate_mesh(x, y, 0.1; mesh_algorithm=mesh_algo, gmsh_path=GMSH_PATH)
```
The argument `0.1` is the refinement parameter for Gmsh; please see the Gmsh documentation for more information about this. This function call starts by processing the Gmsh script, and then builds up a list of elements, nodes, and boundary nodes. Then, using `triangulate`, the result is converting into our data structures. The outputs `T`, `adj`, `adj2v`, `DG`, and `pts` are as before. The only new result is `BN`. This is a `Vector{Vector{Int64}}`, with `BN[1]` representing the boundary nodes. The reason that this is a vector of vectors is so that a boundary can be represented as having multiple boundary segments, i.e. any boundary $\Gamma = \cup_i \Gamma_i$ such that $\Gamma$ is a simple closed curve and $\Gamma_i \cap \Gamma_j = \emptyset$, $i \neq j$.

Let us demonstrate this boundary feature with a more pathological example. We consider a boundary with five segments. It is important that the boundarys are all given in a counter-clockwise orientation, and that the endpoints of the boundaries all line up, with the last segment's endpoint existing in the next segment as its initial point.
```julia
## The first segment 
t = LinRange(0, 1 / 4, 250)
x1 = cos.(2π * t)
y1 = sin.(2π * t)
## The second segment 
t = LinRange(0, -3, 250)
x2 = collect(t)
y2 = repeat([1.0], length(t))
## The third segment 
t = LinRange(1, 0, 250)
x3 = -3.0 .+ (1 .- t) .* sin.(t)
y3 = collect(t)
## The fourth segment 
t = LinRange(0, 1, 250)
x4 = collect(-3.0(1 .- t))
y4 = collect(0.98t)
## The fifth segment 
x5 = [0.073914, 0.0797, 0.1522, 0.1522, 0.2, 0.28128, 0.3659, 0.4127, 0.3922, 0.4068, 0.497, 0.631, 0.728, 0.804, 0.888, 1.0]
y5 = [0.8815, 0.8056, 0.80268, 0.73258, 0.6, 0.598, 0.5777, 0.525, 0.4346, 0.3645, 0.3032, 0.2886, 0.2623, 0.1367, 0.08127, 0.0]
## Now combine the vectors 
x = [x1, x2, x3, x4, x5]
y = [y1, y2, y3, y4, y5]
## Mesh 
tri, BN = generate_mesh(x, y, 0.05; gmsh_path=GMSH_PATH)
T = tri.triangles 
pts = tri.points
## Visualise 
fig = Figure()
ax = Axis(fig[1, 1])
poly!(ax, pts, [collect(T)[i][j] for i in 1:length(T), j in 1:3], color=(:white, 0.0), strokewidth=1)
lines!(ax, pts[:, BN[1]], color=:red, linewidth=4)
lines!(ax, pts[:, BN[2]], color=:blue, linewidth=4)
lines!(ax, pts[:, BN[3]], color=:orange, linewidth=4)
lines!(ax, pts[:, BN[4]], color=:purple, linewidth=4)
lines!(ax, pts[:, BN[5]], color=:darkgreen, linewidth=4)
```

![A Gmsh triangulation](https://github.com/DanielVandH/DelaunayTriangulation.jl/blob/main/test_fig3.png?raw=true)

As you can see, the triangulation successfully works on this non-convex geometry, and we have been able to clearly identify where in `pts` all the boundary nodes are. This makes it especially useful for example for applying boundary conditions in a finite volume code.

# Customisation 

The package is built to allow for customisation. For either `triangulate_bowyer` or `triangulate_berg`, we have the following keywords (these are not all of them):

- `IntegerType`: This is the type used for representing integers (default is `Int64`).
- `EdgeType`: The type used to represent an edge (default is `NTuple{2, IntegerType})`.
- `TriangleType`: The type used to represent a triangle (default is `NTuple{3, IntegerType})`.
- `EdgesType`: The type used to represent a set of edges (default is `Set{EdgeType}`).
- `TrianglesType`: The type used to represent a set of triangles (default is `Set{TriangleType}`).

With these keywords we can allow for a generic interface for these geometric primitives. Below we list these interfaces, followed by an example.

## Triangle interface

Any definition of a triangle must implement the following functions, where we let `T` denote a triangle and `Triangle` is the triangle type, and we assume that `T` has vertices `i`, `j`, and `k`:

- `geti(T::Triangle)`: Returns `i`.
- `getj(T::Triangle)`: Returns `j`.
- `getk(T::Triangle)`: Returns `k`.
- `construct_triangle(::Type{Triangle}, i, j, k)`: Given vertices `i`, `j`, and `k`, and the type `Triangle`, constructs a triangle of type `Triangle` with vertices `i`, `j`, and `k`.
- `integer_type(::Type{Triangle})`: Returns the integer type for the vertices.

For a collection of triangles, it is needed is that the type is iterable and `length(T)` returns the number of triangles (so that `num_triangles` returns the correct number). You also need to be able to `push!` into `T` and use `delete!` to remove a specific triangle. Finally, you `eltype` should remove `Triangle` when applied to `Triangles`.

## Edge interface 

Any definition of an edge must implement the following function, where we let `E` denote an edge and `Edge` is the edge type, and we assume that `E` has vertices `i` and `j`:

- `construct_edge(::Type{Edge}, i, j)`: Given the type `Edge` and vertices `i` and `j`, returns the `Edge` with vertices `i` and `j`.
- `Base.iterate(e::Edge, state...)`: `Edge` must be iterable.

For a collection of edges, it is needed is that the type is iterable, and you must define `Edges()` which should return an empty structure that edges can be added into. Moreover, it must be possible to `push!` `Edge`s into the `Edges` object. Additionally, it must also possible to use `delete!` to remove a specific edge; when `delete!` is used, it should make sure that all instances of the edges are removed, in case you use e.g. a `Vector` data type where duplicates are allowed (unlike a `Set`).. Finally, it must be possible to sample an edge using `rand`.
## Point interface

Any definition of a point must implement the following functions, where we let `p` denote a point and `Point` is the point type, and we assume that `p` has coordinates `(x, y)`:

- `getx(p::Point)`: Returns `x`. The default implementation for this is `getx(p) = p[1]`.
- `gety(p::Point)`: Returns `y`. The default implementation for this is `gety(p) = p[2]`.
- `number_type(p::Point)`: Returns the number type used for representing the coordinates. The default implementation for this is `number_type(p) = number_type(p[1])`, `number_type(::T) where {T<:Number} = T`.

For a collection of points, we require:

- `_eachindex(pts::Points)`: Returns an iterator over the indices for the points. The default implementation is `_eachindex(pts) = eachindex(pts)`.
- `_get_point(pts::Points, i)`: Returns the `i`th point.
- `num_points(pts::Points)`: Returns the number of points. The default implementation is `num_points(pts) = length(pts)`.
- `number_type(pts::Points)`: Returns the number type used for representing the coordinates. The default implementation for this is `number_type(pts) = number_type(pts[1])`.
- `Base.firstindex(pts::Points)`: This should return the index of the first point.

You can also use a matrix of points (as illustrated when we plot the triangulation in the last section).

## Example 

Let us now give an example of how we can implement these interfaces.

First, we must define the types.
```julia 
struct CustomEdge
    i::Int32
    j::Int32
end
struct CustomTriangle
    i::Int32
    j::Int32
    k::Int32
end
struct CustomEdges
    e::Vector{CustomEdge}
end
struct CustomTriangles
    V::Vector{CustomTriangle}
end    
struct CustomPoint
    px::Float64 
    py::Float64 
end
struct CustomPoints{N} 
    pts::NTuple{N, CustomPoint} 
end
```
Note that for `CustomPoint`, `Float64` is the only type that is compatible with ExactPredicates.jl. Now, let us define all the functions needed.
```julia
DT = DelaunayTriangulation 

# Triangle
DT.geti(CT::CustomTriangle) = CT.i
DT.getj(CT::CustomTriangle) = CT.j
DT.getk(CT::CustomTriangle) = CT.k
DT.construct_triangle(::Type{CustomTriangle}, i, j, k) = CustomTriangle(i, j, k)
DT.integer_type(::Type{CustomTriangle}) = Int32

# Triangles 
DT.construct_triangles(::Type{CustomTriangles}) = CustomTriangles(CustomTriangle[])
DT.triangle_type(::Type{CustomTriangles}) = CustomTriangle 
Base.push!(T::CustomTriangles, V::Vararg{CustomTriangle, N}) where {N} = push!(T.V, V...)
function Base.delete!(T::CustomTriangles, V::CustomTriangle) 
    triangle_idx = findfirst(==(V), T.V)
    if triangle_idx !== nothing # V might not be in T - don't do anything in that case. Don't worry about checking the triangle up to rotation, that is already done in the main code.
        deleteat!(T.V, triangle_idx)
    end 
    return nothing
end
Base.length(T::CustomTriangles) = length(T.V)
Base.iterate(T::CustomTriangles, state...) = Base.iterate(T.V, state...)

# Edge
DT.construct_edge(::Type{CustomEdge}, i, j) = CustomEdge(i, j)
Base.length(::CustomEdge) = 2 # Not needed, but helps the iterator
Base.iterate(E::CustomEdge, state=1) = state > 2 ? nothing : (state == 1 ? (E.i, state+1) : (E.j, state+1))

# Edges 
CustomEdges() = CustomEdges(CustomEdge[])
DT.push!(Es::CustomEdges, E::Vararg{CustomEdge, N}) where {N} = push!(Es.e, E...)
Base.length(E::CustomEdges) = length(E.e) # Not needed, but helps the iterator 
Base.iterate(E::CustomEdges, state...) = Base.iterate(E.e, state...)
Base.delete!(Es::CustomEdges, E::CustomEdge) = deleteat!(Es.e, findall(==(E), Es.e)) # No need to check if the edge is not there - an error should given in that case since it should not be possible. Use findall to avoid the problem of duplicate edges.
using Random 
Random.rand(Es::CustomEdges) = rand(Es.e)

# Point
DT.getx(p::CustomPoint) = p.px 
DT.gety(p::CustomPoint) = p.py 
DT.number_type(::CustomPoint) = Float64 

# Points 
DT._eachindex(::CustomPoints{N}) where {N} = Base.OneTo(N)
DT._get_point(pts::CustomPoints, i) = pts.pts[i]
DT.number_type(::CustomPoints) = Float64
Base.firstindex(::CustomPoints) = 1
```

Now having defined all these new methods, let's triangulate.
```julia
x_data = 15rand(250)
y_data = 15rand(250)
pts = CustomPoints(Tuple(CustomPoint(x, y) for (x, y) in zip(x_data, y_data)))
_seed = 922881
Random.seed!(_seed) # Keep the same insertion order for later
T, adj, adj2v, DG = triangulate_bowyer(pts; 
    IntegerType = Int32,
    EdgeType = CustomEdge, 
    TriangleType = CustomTriangle, 
    EdgesType = CustomEdges,
    TrianglesType = CustomTriangles)
```

To make sure this is correct, we can first convert all the structures to standard structures and compare with the triangulation computed with the default structures.
```julia
# Convert T, adj, adj2v, DG to standard structures for comparison
converted_T = Set{NTuple{3, Int64}}(indices(T) for T in T)
converted_adj = DT.Adjacent{Int64, NTuple{2, Int64}}()
for ((u, v), k) in adjacent(adj)
    DT.add_edge!(converted_adj,u,v,k)
end
converted_adj2v = DT.Adjacent2Vertex{Int64, Set{NTuple{2,Int64}}, NTuple{2, Int64}}()
for (k, es) in adjacent2vertex(adj2v)
    for (i, j) in es
        DT.add_edge!(converted_adj2v, k, i, j)
    end
end
converted_DG = DT.DelaunayGraph{Int64}()
for i in DT._eachindex(pts)
    DT.add_neighbour!(converted_DG, i, get_neighbour(DG, i)...)
end
DT.add_neighbour!(converted_DG, DT.BoundaryIndex, get_neighbour(DG, DT.BoundaryIndex)...)

# Now compare
using Test
_pts = [x_data'; y_data']
Random.seed!(_seed)
_T, _adj, _adj2v, _DG = triangulate_bowyer(_pts)
@test DT.compare_unconstrained_triangulations(converted_T, converted_adj, converted_adj2v, converted_DG, _T, _adj, _adj2v, _DG)
```
This all applies equally as well to `triangulate_berg`.

```julia
Random.seed!(_seed) # Keep the same insertion order for later
(T, adj, adj2v, DG), _ = triangulate_berg(pts; 
    IntegerType = Int32,
    EdgeType = CustomEdge, 
    TriangleType = CustomTriangle, 
    EdgesType = CustomEdges,
    TrianglesType = CustomTriangles)
converted_T = Set{NTuple{3, Int64}}(indices(T) for T in T)
converted_adj = DT.Adjacent{Int64, NTuple{2, Int64}}()
for ((u, v), k) in adjacent(adj)
    DT.add_edge!(converted_adj,u,v,k)
end
converted_adj2v = DT.Adjacent2Vertex{Int64, Set{NTuple{2,Int64}}, NTuple{2, Int64}}()
for (k, es) in adjacent2vertex(adj2v)
    for (i, j) in es
        DT.add_edge!(converted_adj2v, k, i, j)
    end
end
converted_DG = DT.DelaunayGraph{Int64}()
for i in DT._eachindex(pts)
    DT.add_neighbour!(converted_DG, i, get_neighbour(DG, i)...)
end
DT.add_neighbour!(converted_DG, DT.BoundaryIndex, get_neighbour(DG, DT.BoundaryIndex)...)
Random.seed!(_seed)
_T, _adj, _adj2v, _DG = triangulate_berg(_pts)
@test DT.compare_unconstrained_triangulations(converted_T, converted_adj, converted_adj2v, converted_DG, _T, _adj, _adj2v, _DG)
```

# A note on Voronoi tessellations 

There is some existing code for creating Voronoi tessellations in `src/voronoi.jl`, but it is currently not used since I am still trying to enable polygons to be chopped off so that they live only inside the convex hull. You can still use `voronoi` though, it just won't trim any polygons. See also the tests (that may not all work yet) in `test/voronioi.jl`. Feel free to ask any questions about the implementation here, it may not all be fully commented yet. See also 
`writeups/voronoi.jl` for some work that I was doing with it, it's not so organised unfortunately (nothing in the `writeups` folder is).
