# DelaunayTriangulation


This is a package for creating (unconstrained) two-dimensional Delaunay triangulations. The great package [ExactPredicates.jl](https://github.com/lairez/ExactPredicates.jl) is used for all geometric predicates. There also exist routines for building the convex hull (from the triangulation), and currently commented-out Voronoi tessellation code (see some discussion at the end). Much of the work in this package is derived from the book *Delaunay Mesh Generation* by Cheng, Dey, and Shewchuk (2013). Point location is implemented using the jump-and-march algorithm of Mücke, Saias, and Zhu (1999); see `jump_and_march`.

The package has two algorithms for computing Delauanay triangulations, namely de Berg's method from de Berg et al. (1999), and the Bowyer-Watson algorithm as presented by Cheng, Dey, and Shewchuk (2013). de Berg's method is the slowest of the two, and has less features, so we recommend the Bowyer-Watson algorithm in the function `triangulate_bowyer`. Keep reading for more examples, and see the tests for more features.

# Getting started

## de Berg's method

It is easy to construct a triangulation. Here is a simple example, where we use the method outlined in the book by de Berg et al. (1999) to construct the triangulation; this triangulation adds points one at a time, constructing the Delaunay triangulation $\mathcal DT(\mathcal P_n)$ from $\mathcal DT(\mathcal P_{n-1})$, $\mathcal P_n = \mathcal P_{n-1} \cup \{p_n\}$ using edge flips. Using `triangulate_berg`, this triangulation can be computed:

```julia
using DelaunayTriangulation
p1 = [0.0, 1.0]
p2 = [3.0, -1.0]
p3 = [2.0, 0.0]
p4 = [-1.0, 2.0]
p5 = [4.0, 2.0]
p6 = [-2.0, -1.0]
pts = [p1, p2, p3, p4, p5, p6]
T, adj, adj2v, DG, HG = triangulate_berg(pts)
```
This function `triangulate_berg` creates the triangulation for the given `pts`, returning:
- `T`: This is the set of triangles, of default type `Set{NTuple{3, Int64}}`. By default, triangles are represented as tuples of indices `(i, j, k)` so that `(pts[i], pts[j], pts[k])` are the points representing the triangle. All triangles `T` are positively oriented.
- `adj`: This is the `Adjacent` map, mapping edges `(i, j)` to a vertex `k` such that `(i, j, k)` is a positively oriented (counter-clockwise) triangle that exists in `T`. If no such triangle exists, `-4` is returned (this constant is defined in `DelaunayTriangulation.DefaultAdjacentValue`). You can extract this vertex using `get_edge(adj, i, j)`. By default, edges are represented as `NTuple{2, Int64}`.
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
p1 = [0.0, 1.0]
p2 = [3.0, -1.0]
p3 = [2.0, 0.0]
p4 = [-1.0, 2.0]
p5 = [4.0, 2.0]
p6 = [-2.0, -1.0]
pts = [p1, p2, p3, p4, p5, p6]
T, adj, adj2v, DG = triangulate_bowyer(pts)
```

The definitions of `T`, `adj`, `adj2v`, and `DG` are the same as above. One nice feature of the Bowyer-Watson algorithm as implemented is the use of ghost triangles, i.e. a ghost triangle adjoined to a boundary edge with a vertex out at infinity. This makes it easy to locate points in the boundary, defining a unique ghost triangle that they reside in, and makes it easy to traverse the boundary. You can keep these ghost triangles by setting `trim=false`:

```julia
T, adj, adj2v, DG = triangulate_bowyer(pts; trim=false)
```

Alternatively, you can add it after the fact as follows:

```julia
DelaunayTriangulation.add_ghost_triangles!(T, adj, adj2v, DG)
```

(Use `DelaunayTriangulation.remove_ghost_triangles!(T, adj, adj2v, DG)` to remove them.) 

Visualising these triangulations is easy:

```julia
pts = [15rand(2) for _ in 1:500]
T, adj, adj2v, DG = triangulate_bowyer(pts)
using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1])
poly!(ax, Tuple.(pts), [collect(T)[i][j] for i in 1:length(T), j in 1:3], color = (:white, 0.0), strokewidth = 2)
```

![A triangulation](https://github.com/DanielVandH/DelaunayTriangulation.jl/blob/main/test_fig.png?raw=true)