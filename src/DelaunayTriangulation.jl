"""
    DelaunayTriangulation

Module for computing Delaunay triangulations and Voronoi tessellations in two dimensions. There 
are many features available, including:

- Unconstrained and constrained Delaunay triangulations - see `triangulate`.
- Computation of Voronoi tessellations, clipped Voronoi tessellations (clipped to the convex hull), and centroidal Voronoi tessellations - see `voronoi`.
- Mesh refinement, with support custom angles and area constraints, as well as refinement of curve-bounded domains - see `refine!`.
- Dynamic updates such as point insertion and segment insertion - see e.g. `add_point!`, `delete_vertex!`, and `add_segment!`.
- Computation of convex hulls - see `convex_hull` (or use `get_convex_hull` on a computed triangulation).
- Triangulation of convex polygons - see `triangulate_convex`. Lattices can also be triangulated, see `triangulate_rectangle`.
- Point location - see `find_triangle` and `find_polygon`.
- Computation of the pole of inaccessibility - see `DelaunayTriangulation.pole_of_inaccessibility`.

See the documentation for more, and the package's associated GitHub repository's README.
"""
module DelaunayTriangulation

include("setup.jl")

import ExactPredicates as EP
import AdaptivePredicates as AP
import EnumX
import Random

abstract type AbstractPredicateKernel end # needs to be defined early for use in data_structures.jl

include("data_structures.jl")
include("predicates.jl")
include("geometric_primitives.jl")
include("utils.jl")
include("algorithms.jl")
include("validation.jl")
include("exports.jl")
include("public.jl")

end
