# DelaunayTriangulation

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://DanielVandH.github.io/DelaunayTriangulation.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://DanielVandH.github.io/DelaunayTriangulation.jl/dev/)
[![Coverage](https://codecov.io/gh/DanielVandH/DelaunayTriangulation.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/DanielVandH/DelaunayTriangulation.jl)
[![DOI](https://zenodo.org/badge/540660309.svg)](https://zenodo.org/badge/latestdoi/540660309)

This is a package for constructing Delaunay triangulations of planar point sets, with partial support for constrained Delaunay triangulations. The support is partial in that we support them when creating via Gmsh (see the `generate_mesh` function), but we do not yet have a method for constructing constrained Delaunay triangulations with a given piecewise linear complex; such a method is in development. All geometric predicates are computed via ExactPredicates.jl.

Much of the work in this package is derived from the book *Delaunay Mesh Generation* by Cheng, Dey, and Shewchuk (2013).

## Example

Here we give a simple example, but you should refer to the docs or tests for more details.

```julia
using DelaunayTriangulation, CairoMakie
pts = rand(2, 500) # This is a matrix, but we could just as well do a vector of vectors or a vector of tuples
tri = triangulate(pts)
fig, ax, sc = triplot(tri; show_convex_hull = true, convex_hull_linewidth=4)
```

![A triangulation](https://github.com/DanielVandH/DelaunayTriangulation.jl/blob/main/test/figures/custom_interface_testing.png?raw=true)

You can also use Gmsh for constrained triangulations (support for this without having to rely on Gmsh is coming). Here's an example of triangulationg a multiply-connected domain - please see the relevant section of the docs for more information.

```julia
using DelaunayTriangulation, CairoMakie
x1 = [collect(LinRange(0, 2, 4)),
    collect(LinRange(2, 2, 4)),
    collect(LinRange(2, 0, 4)),
    collect(LinRange(0, 0, 4))]
y1 = [collect(LinRange(0, 0, 4)),
    collect(LinRange(0, 6, 4)),
    collect(LinRange(6, 6, 4)),
    collect(LinRange(6, 0, 4))]
r = 0.5
h = k = 0.6
θ = LinRange(2π, 0, 50)
x2 = [h .+ r .* cos.(θ)]
y2 = [k .+ r .* sin.(θ)]
r = 0.2
h = 1.5
k = 0.5
x3 = [h .+ r .* cos.(θ)]
y3 = [k .+ r .* sin.(θ)]
x4 = reverse(reverse.([collect(LinRange(1, 1.5, 4)),
    collect(LinRange(1.5, 1.5, 4)),
    collect(LinRange(1.5, 1, 4)),
    collect(LinRange(1, 1, 4))]))
y4 = reverse(reverse.([collect(LinRange(2, 2, 4)),
    collect(LinRange(2, 5, 4)),
    collect(LinRange(5, 5, 4)),
    collect(LinRange(5, 2, 4))]))
x5 = [reverse([0.2, 0.5, 0.75, 0.75, 0.2, 0.2])]
y5 = [reverse([2.0, 2.0, 3.0, 4.0, 5.0, 2.0])]
x = [x1, x2, x3, x4, x5]
y = [y1, y2, y3, y4, y5]
tri = generate_mesh(x, y, 0.2)
fig, ax, sc = triplot(tri; show_ghost_edges=true, convex_hull_linestyle=:solid, convex_hull_linewidth=4)
xlims!(ax, -0.5, 2.5)
ylims!(ax, -0.5, 6.5)
```

![A triangulation](https://github.com/DanielVandH/DelaunayTriangulation.jl/blob/main/test/figures/gmsh_example_3.png?raw=true)