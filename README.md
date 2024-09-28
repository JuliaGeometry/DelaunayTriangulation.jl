# DelaunayTriangulation

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGeometry.github.io/DelaunayTriangulation.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGeometry.github.io/DelaunayTriangulation.jl/dev/)
[![Coverage](https://codecov.io/gh/JuliaGeometry/DelaunayTriangulation.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaGeometry/DelaunayTriangulation.jl)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.07174/status.svg)](https://doi.org/10.21105/joss.07174)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

This is a package for constructing Delaunay triangulations and Voronoi tessellations of planar point sets. Supports unconstrained and constrained triangulations, weighted triangulations, mesh refinement, triangulation of curve bounded domains, Voronoi tessellations, power diagrams, and clipped and centroidal Voronoi tessellations. To install the package, do

```julia
julia>] add DelaunayTriangulation
```

Many features are available, some of these being:

- Unconstrained and constrained Delaunay triangulations, supporting many types of domains.
- Weighted Delaunay triangulations.
- Computation of Voronoi tessellations, clipped Voronoi tessellations where the Voronoi tiles get clipped to a convex polygon, and centroidal Voronoi tessellations where each Voronoi tile's generator is the tile's centroid.
- Power diagrams.
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
# using Pkg; Pkg.add(["DelaunayTriangulation", "CairoMakie", "StableRNGs"])
using DelaunayTriangulation, CairoMakie, StableRNGs

# Unconstrained 
points = rand(2, 50)
tri1 = triangulate(points) # default predicate kernel is AdaptiveKernel()

# Voronoi example 
vorn2 = voronoi(tri1)

# Clipped Voronoi 
vorn3 = voronoi(tri1, clip=true, predicates = ExactKernel())

# Smoothed Voronoi 
vorn4 = centroidal_smooth(vorn3; predicates = FastKernel()) # or do voronoi(tri1, clip = true, smooth = true)

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

# Weighted triangulation example
rng = StableRNG(1) # You can pass rngs for reproducibility. Alternatively, just use Random.seed!(n)
points = tuple.(rand(rng, 20), rand(rng, 20))
weights = 3randn(rng, 20)
tri9 = triangulate(points; weights, rng)

# Power diagram example 
points = tuple.(rand(rng, 100), rand(rng, 100))
weights = rand(rng, 100)
tri10 = triangulate(points; weights, rng)
vorn10 = voronoi(tri10) # can also use clip/smooth here 

# Clipped Voronoi example with a generic convex polygon
points = 10randn(rng, 2, 100)
weights = rand(rng, 100)
circ = CircularArc((0.0, 5.0), (0.0, 5.0), (0.0, 0.0)) # clip to a circle 
clip_points = [circ(t) for t in LinRange(0, 1, 100)]
clip_vertices = [1:(length(clip_points)-1); 1]
tri11 = triangulate(points; weights, rng)
vorn11 = voronoi(tri11, clip=true, clip_polygon=(clip_points, clip_vertices), rng=rng)

# Plotting 
fig = Figure(fontsize = 42, size = (2800, 2200))
wh = (width = 600, height = 600)
ax = Axis(fig[1, 1]; title = "Unconstrained", wh...);            triplot!(ax, tri1)
ax = Axis(fig[1, 2]; title = "Voronoi", wh...);                  voronoiplot!(ax, vorn2)
ax = Axis(fig[1, 3]; title = "Clipped Voronoi", wh...);          voronoiplot!(ax, vorn3)
ax = Axis(fig[1, 4]; title = "Centroidal Voronoi", wh...);       voronoiplot!(ax, vorn4)
ax = Axis(fig[2, 1]; title = "Constrained", wh...);              triplot!(ax, tri5)
ax = Axis(fig[2, 2]; title = "Disjoint Constrained", wh...);     triplot!(ax, tri6)
ax = Axis(fig[2, 3]; title = "Curve-Bounded", wh...);            triplot!(ax, tri7)
ax = Axis(fig[2, 4]; title = "Disjoint Curve-Bounded", wh...);   triplot!(ax, tri8)
ax = Axis(fig[3, 1]; title = "Weighted", wh...);                 triplot!(ax, tri9)
ax = Axis(fig[3, 2]; title = "Power Diagram", wh...);            voronoiplot!(ax, vorn10)
ax = Axis(fig[3, 3]; title = "Generic Clipped Voronoi", wh...);  voronoiplot!(ax, vorn11)
```

![](readme.png)

## Citing DelaunayTriangulation.jl

If you use DelaunayTriangulation.jl, please cite it. There is a JOSS paper published at [https://doi.org/10.21105/joss.07174](https://doi.org/10.21105/joss.07174). The BibTeX entry for this paper is:

```bibtex
@article{VandenHeuvel_DelaunayTriangulation_jl_A_Julia_2024,
    author = {VandenHeuvel, Daniel J.},
    doi = {10.21105/joss.07174},
    journal = {Journal of Open Source Software},
    month = sep,
    number = {101},
    pages = {7174},
    title = {{DelaunayTriangulation.jl: A Julia package for Delaunay triangulations and Voronoi tessellations in the plane}},
    url = {https://joss.theoj.org/papers/10.21105/joss.07174},
    volume = {9},
    year = {2024}
}
```

## License 

DelaunayTriangulation.jl is provided under a [MIT license](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/blob/main/LICENSE).

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

In addition to these Julia packages, software packages in other programming languages are available, such as:

- [Triangle](https://www.cs.cmu.edu/~quake/triangle.html): A C library for generating (constrained) Delaunay triangulations, Voronoi tessellations, and refined meshes. Also supports weighted triangulations. As mentioned above, there is a Julia interface to this library through TriangleMesh.jl, Triangulate.jl, or Triangle.jl. 
- [CGAL](https://www.cgal.org/): A C++ library which, among many other features relevant in computational geometry, supports Delaunay triangulations, spherical triangulations, periodic triangulations, constrained triangulations, 3D Delaunay triangulations, Voronoi tessellations, mesh refinement, surface triangulations, and more.
- [Gmsh](https://gmsh.info/): Gmsh is a finite element mesh generator with many features, providing constrained Delaunay triangulations (or, in addition to triangles, elements can also be mixed with e.g. squares) and mesh refinement with many options. Provides a set of meshing algorithms to choose from. Supports both 2D and 3D meshing. In addition to simple segments, you can also provide curves for boundaries and surfaces in the 3D case. Has interfaces in several languages, including [Gmsh.jl](https://github.com/JuliaFEM/Gmsh.jl).
- [SciPy](https://docs.scipy.org/doc/scipy/tutorial/spatial.html): A Python library providing Delaunay triangulations, Voronoi tessellations, and convex hulls in $\mathbb R^N$. Does not support constrained triangulations or mesh refinement.
- [MATLAB](https://uk.mathworks.com/help/matlab/computational-geometry.html?s_tid=CRUX_lftnav): MATLAB has built-in functions for Delaunay triangulations and Voronoi tessellations in two and three dimensions. Additionally, it supports constrained triangulations and mesh refinement.
- [Qhull](http://www.qhull.org/): A C library for computing Delaunay triangulations, Voronoi tessellations, and convex hulls in $\mathbb R^N$. No support is provided for constrained triangulations or mesh refinement, but it does provide furthest-site Delaunay triangulations and furthest-site Voronoi tessellations. The packages MiniQhull.jl and DirectQhull.jl wrap Qhull for use in Julia.
- [Delaunator](https://github.com/mapbox/delaunator): A JavaScript library that provides support for fast construction of two-dimensional unconstrained Delaunay triangulations.
- [Spade](https://github.com/Stoeoef/spade): A Rust library providing support for Delaunay triangulations and constrained Delaunay triangulations, mesh refinement, and Voronoi tessellations in the plane. 
- [Acute](https://www.cise.ufl.edu/~ungor/aCute/index.html): A C library that builds on top of Shewchuk's Triangle library, being the first of its kind to not only allow for minimum angle constraints in mesh refinement, but also for maximum angle constraints.
- [Tinfour](https://github.com/gwlucastrig/Tinfour): A Java library for Delaunay triangulations and constrained Delaunay triangulations in the plane, with a lot of features for interpolation.
- [JTS](https://github.com/locationtech/jts): A Java library which, among many other features for working with vectory geometries, supports Delaunay triangulations, constrained Delaunay triangulations, mesh refinement (at least to make a conformal triangulation), and Voronoi tessellations.
- [Voro++](https://math.lbl.gov/voro++/): A C++ library for constructing Voronoi tessellations, power diagrams, and clipped tessellations. 
- [Stellar](https://people.eecs.berkeley.edu/~jrs/stellar/): A C library for constructing three-dimensional Delaunay triangulations, providing the ability to efficiently refine the mesh.

Compared to all these other libraries, and only in the context of planar triangulations, DelaunayTriangulation.jl is one of the most developed in terms of the features provided, except possibly with the exception of CGAL and Gmsh who provide many features although neither are in the public domain (CGAL being GPL v3+ and Gmsh being GPL v2+), unlike DelaunayTriangulation.jl which is MIT licensed. A tabular comparison of all these packages is given below (focusing only on two dimensional meshing). If there are any errors in this comparison, please let me know. Also please note that the features listed are not intended to be exhaustive, but rather to give a general idea of the capabilities of each package, and certainly not all software packages are listed here. 

<!--
<div style="font-size: 14px; background-color: white;">
<table><thead><tr><th></th><th></th><th></th><th colspan="10">Delaunay Triangulations</th><th colspan="5">Voronoi Tessellations</th><th></th><th></th><th></th></tr></thead><tbody><tr><td>Package</td><td>Language</td><td>Open Source License</td><td>Unconstrained</td><td>Constrained</td><td>Weighted</td><td>Refinement<br></td><td>Curve-Bounded</td><td>Dynamic Updates</td><td>Convex Polygons</td><td>Periodic</td><td>Spherical</td><td>Furthest-site</td><td>Standard</td><td>Clipped (Rect. or Hull)</td><td>Generic Clipped</td><td>Centroidal</td><td>Power</td><td>Point Location</td><td>Robust Predicates</td><td>Customisable Interface</td></tr><tr><td>DelaunayTriangulation.jl</td><td>Julia</td><td>✅ (MIT)</td><td>✅</td><td>✅</td><td>✅</td><td>✅</td><td>✅</td><td>✅</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>✅</td><td>✅</td><td>✅</td><td>✅</td><td>✅</td><td>✅</td><td>✅</td></tr>
<tr><td>VoronoiDelaunay.jl</td><td>Julia</td><td>✅ (MIT)</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>✅</td><td>❌</td></tr><tr><td>VoronoiCells.jl</td><td>Julia</td><td>✅ (MIT)</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td></tr><tr><td>Delaunay.jl</td><td>Julia</td><td>✅ (MIT)</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td></tr><tr><td>MiniQhull.jl</td><td>Julia</td><td>✅ (MIT)</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td></tr>
<tr><td>DirectQhull.jl</td><td>Julia</td><td>✅ (MIT)</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td></tr><tr><td>Delaunator.jl</td><td>Julia</td><td>✅ (MIT)</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td></tr><tr><td>TriangleMesh.jl</td><td>Julia</td><td>✅ (MIT)</td><td>✅</td><td>✅</td><td>❌</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td></tr><tr><td>Triangulate.jl</td><td>Julia</td><td>✅ (MIT)</td><td>✅</td><td>✅</td><td>❌</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td></tr>
<tr><td>Triangle.jl</td><td>Julia</td><td>✅ (MIT)</td><td>✅</td><td>✅</td><td>❌</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td></tr><tr><td>GMT.jl</td><td>Julia</td><td>✅ (MIT)</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td></tr><tr><td>Quickhull.jl</td><td>Julia</td><td>✅ (MIT)</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td></tr><tr><td>Triangle</td><td>C</td><td>❌ (Non-standard)</td><td>✅</td><td>✅</td><td>✅</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>✅</td><td>❌</td></tr>
<tr><td>CGAL</td><td>C++</td><td>❌ (GPL v3+/GNU GPL/Commerical)</td><td>✅</td><td>✅</td><td>✅</td><td>✅</td><td>❌</td><td>✅</td><td>❌</td><td>✅</td><td>✅</td><td>❌</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>✅</td><td>✅</td><td>❌</td></tr><tr><td>Gmsh</td><td>Own/C++</td><td>❌ (GPL v2+)</td><td>✅</td><td>✅</td><td>❌</td><td>✅</td><td>✅</td><td>❌</td><td>❌</td><td>✅</td><td>✅</td><td>❌</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td></tr><tr><td>SciPy</td><td>Python</td><td>✅ (BSD 3-Clause)</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td></tr><tr><td>MATLAB</td><td>MATLAB</td><td>❌ (Commercial)</td><td>✅</td><td>✅</td><td>❌</td><td>✅</td><td>❌</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>✅</td><td>❌</td></tr>
<tr><td>Qhull</td><td>C</td><td>✅ (Permissive)</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td></tr><tr><td>Delaunator</td><td>JavaScript</td><td>✅ (ISC)</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td></tr><tr><td>Spade</td><td>Rust</td><td>✅ (MIT/Apache 2.0)</td><td>✅</td><td>✅</td><td>❌</td><td>✅</td><td>❌</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td></tr><tr><td>Acute</td><td>C</td><td>❌ (Non-standard)</td><td>✅</td><td>✅</td><td>✅</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>✅</td><td>❌</td></tr>
<tr><td>Tinfour</td><td>Java</td><td>✅ (Apache 2.0)</td><td>✅</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>✅</td><td>❌</td></tr><tr><td>JTS</td><td>Java</td><td>✅ (EPL 2.0/EDL 1.0)</td><td>✅</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>❌</td><td>✅</td><td>❌</td></tr></tbody></table>
</div>
-->

![Comparison](softwarecomparison.png)
