# DelaunayTriangulation

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGeometry.github.io/DelaunayTriangulation.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGeometry.github.io/DelaunayTriangulation.jl/dev/)
[![Coverage](https://codecov.io/gh/JuliaGeometry/DelaunayTriangulation.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaGeometry/DelaunayTriangulation.jl)
[![DOI](https://zenodo.org/badge/540660309.svg)](https://zenodo.org/badge/latestdoi/540660309)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

This is a package for constructing Delaunay triangulations and Voronoi tessellations of planar point sets. Supports unconstrained and constrained triangulations, mesh refinement, triangulation of curve bounded domains, Voronoi tessellations, and clipped and centroidal Voronoi tessellations. To install the package, do

```julia
julia>] add DelaunayTriangulation
```

Many features are available, some of these being:

- Unconstrained and constrained Delaunay triangulations, supporting many types of domains.
- Computation of Voronoi tessellations, clipped Voronoi tessellations where the Voronoi tiles get clipped to the convex hull, and centroidal Voronoi tessellations where each Voronoi tile's generator is the tile's centroid.
- Mesh refinement, with support for custom angle and area constraints, as well as refinement of curve-bounded domains.
- Dynamic point insertion, point deletion, and segment insertion, amongst many other operations.
- Computation of convex hulls.
- Triangulation of convex polygons.
- Point location.
- Computation of the pole of inaccessibility.
- The interface for defining geometric primitives is fully customisable.

To ensure that the algorithms are robust, by default we make use of [AdaptivePredicates.jl](https://github.com/JuliaGeometry/AdaptivePredicates.jl) to use 
adaptive arithmetic for all geometric predicates in this package. This choice can be configured, allowing for the additional choices of [ExactPredicates.jl](https://github.com/lairez/ExactPredicates.jl) or no adaptive or exact arithmetic at all; see the documentation. Much of the work in this package is derived from the book [*Delaunay Mesh Generation* by Cheng, Dey, and Shewchuk (2013)](https://people.eecs.berkeley.edu/~jrs/meshbook.html). Please see the documentation for much more information.

Some examples are below (and in the documentation), but if you would also like to see how DelaunayTriangulation.jl is used in other packages, see [FiniteVolumeMethod.jl](https://github.com/SciML/FiniteVolumeMethod.jl) (solving 2D PDEs) and [NaturalNeighbours.jl](https://github.com/DanielVandH/NaturalNeighbours.jl) (scattered data interpolation).

If you would like to make an issue or contribute, please see the [CONTRIBUTING.md](CONTRIBUTING.md) file. For feature requests, please see the discussions tab.

## Example Usage 
Here is a quick example of some ways the package can be used. As mentioned, see the docs for many examples.

```julia
# using Pkg; Pkg.add(["DelaunayTriangulation", "CairoMakie"])
using DelaunayTriangulation, CairoMakie

# Unconstrained 
points = rand(2, 50)
tri1 = triangulate(points) # default predicate kernel is AdaptiveKernel()

# Voronoi example 
vorn2 = voronoi(tri1)

# Clipped Voronoi 
vorn3 = voronoi(tri1, clip=true, predicates = ExactKernel())

# Smoothed Voronoi 
vorn4 = centroidal_smooth(vorn3, predicates = FastKernel())

# Constrained example with refinement 
boundary_points = [(0.0, 0.0), (1.0, 0.0), (1.0, 0.3), (0.5, 0.3),
    (0.3, 0.7), (0.1, 1.0), (0.0, 1.0), (0.0, 0.0)]
boundary_nodes, points = convert_boundary_points_to_indices(boundary_points)
tri5 = triangulate(points; boundary_nodes)
refine!(tri5; max_area=1e-2get_area(tri5))

# Disjoint constrained example with refinement 
boundary_points = [
    [[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]],
    [[(0.3, 0.3), (0.3, 0.7), (0.7, 0.7), (0.7, 0.3), (0.3, 0.3)]],
    [[(1.2, 0.0), (1.4, 0.0), (1.4, 1.2), (0.0, 1.2), (0.0, 1.1),
        (1.1, 1.1), (1.1, 0.0), (1.2, 0.0)]]
]
boundary_nodes, points = convert_boundary_points_to_indices(boundary_points)
tri6 = triangulate(points; boundary_nodes)
refine!(tri6; max_area=1e-2get_area(tri6))

# Curve-bounded example
using DelaunayTriangulation: EllipticalArc
ellipse = EllipticalArc((1.0, 0.0), (1.0, 0.0), (0.0, 0.0), 1.0, 2.0, 0.0)
tri7 = triangulate(NTuple{2,Float64}[]; boundary_nodes=[ellipse])
refine!(tri7; max_area=1e-2get_area(tri7))

# Disjoint curve-bounded example 
ellipse = EllipticalArc((7.0, 3.5), (7.0, 3.5), (0.0, 3.5), 7.0, 10.0, 0.0)
catmull_cp = [(0.0, 0.0), (-2.0, -1.0), (-4.0, 0.0), (-5.0, 2.0), (-1.0, 4.0), (0.0, 3.0),
    (1.0, 4.0), (5.0, 2.0), (4.0, 0.0), (2.0, -1.0), (0.0, 0.0)]
catmull_spl = CatmullRomSpline(catmull_cp)
circle = CircularArc((0.5, 1.5), (0.5, 1.5), (0.0, 1.0))
circle2 = CircularArc((0.1, 1.5), (0.1, 1.5), (0.0, 1.0), positive=false)
points = [(-4.0, -10.0), (4.0, -10.0), (4.0, -7.0), (-4.0, -7.0)]
square = [1, 2, 3, 4, 1]
boundary_nodes = [[square], [[ellipse]], [[catmull_spl]], [[circle]], [[circle2]]]
tri8 = triangulate(points; boundary_nodes)
refine!(tri8; max_area=1e-2get_area(tri8)) # could also use find_polygon to help define a custom refinement function for each shape

# Plotting 
fig = Figure(fontsize = 42, size = (2800, 1480))
ax = Axis(fig[1, 1], title="Unconstrained", width=600,height=600);            triplot!(ax, tri1)
ax = Axis(fig[1, 2], title="Voronoi", width=600,height=600);                  voronoiplot!(ax, vorn2)
ax = Axis(fig[1, 3], title="Clipped Voronoi", width=600,height=600);          voronoiplot!(ax, vorn3)
ax = Axis(fig[1, 4], title="Centroidal Voronoi", width=600,height=600);       voronoiplot!(ax, vorn4)
ax = Axis(fig[2, 1], title="Constrained", width=600,height=600);              triplot!(ax, tri5)
ax = Axis(fig[2, 2], title="Disjoint Constrained", width=600,height=600);     triplot!(ax, tri6)
ax = Axis(fig[2, 3], title="Curve-Bounded", width=600,height=600);            triplot!(ax, tri7)
ax = Axis(fig[2, 4], title="Disjoint Curve-Bounded", width=600,height=600);   triplot!(ax, tri8)
```

![](readme.png)

## Similar Packages

This is not the only Delaunay triangulation package available. Some others are:

- [VoronoiDelaunay.jl](https://github.com/JuliaGeometry/VoronoiDelaunay.jl): A pure Julia library that constructs planar triangulations and tessellations like in this package, although no support for constrained triangulations / mesh refinement or clipped / centroid tessellations. Restricts points to $[1, 2] \times [1, 2]$.
- [VoronoiCells.jl](https://github.com/JuliaGeometry/VoronoiCells.jl): A pure Julia library that extends VoronoiDelaunay.jl. This package provides useful tools for constructing and working with Voronoi tessellations. Supports clipping Voronoi cells to a specified rectangle. Like VoronoiDelaunay.jl, restricts points to $[1, 2] \times [1, 2]$.
- [Delaunay.jl](https://github.com/eschnett/Delaunay.jl): Wraps Python's main Delaunay triangulation library, [`scipy.spatial.Delaunay`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Delaunay.html), for computing Delaunay triangulations in $\mathbb R^N$. I don't believe constrained triangulations or mesh refinement is available here.
- [MiniQhull.jl](https://github.com/gridap/MiniQhull.jl): Wraps [Qhull](http://www.qhull.org/) for computing unconstrained Delaunay triangulations in $\mathbb R^N$. No support is provided for mesh refinement.
- [DirectQhull.jl](https://github.com/JuhaHeiskala/DirectQhull.jl/): Similar to MiniQhull.jl, although also provides support for convex hulls and Voronoi tessellations from Qhull.
- [Delaunator.jl](https://github.com/JuliaGeometry/Delaunator.jl): A pure Julia library modelled after the [JavaScript Delaunator library](https://github.com/mapbox/delaunator). This package can construct unconstrained triangulations of planar point sets. No support is available for constrained triangulations or mesh refinement, although support exists for computing the dual Voronoi tessellation. Centroidal tessellations are not implemented, although the Voronoi cells can be clipped to a bounding box. 
- [TriangleMesh.jl](https://github.com/konsim83/TriangleMesh.jl), [Triangulate.jl](https://github.com/JuliaGeometry/Triangulate.jl), [Triangle.jl](https://github.com/cvdlab/Triangle.jl): Interfaces to [Shewchuk's Triangle library](https://www.cs.cmu.edu/~quake/triangle.html).
- [TetGen.jl](https://github.com/JuliaGeometry/TetGen.jl): This is for Delaunay tetrahedralisation, wrapping [TetGen](https://wias-berlin.de/software/index.jsp?id=TetGen).
- [GMT.jl](https://github.com/GenericMappingTools/GMT.jl): A wrapper of [GMT](https://github.com/GenericMappingTools/gmt), allowing for [unconstrained Delaunay triangulations in two dimensions](https://www.generic-mapping-tools.org/GMTjl_doc/documentation/modules/triangulate/index.html#triangulate), and for [spherical triangulation, i.e. triangulation of points lying on a sphere](https://www.generic-mapping-tools.org/GMTjl_doc/documentation/modules/sphtriangulate/index.html#sphtriangulate).
- [Quickhull.jl](https://github.com/augustt198/Quickhull.jl): A pure Julia library for unconstrained triangulations, Voronoi tessellations, and convex hulls in $N$ dimensions.

