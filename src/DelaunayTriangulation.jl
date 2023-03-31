module DelaunayTriangulation

####################################################
##
## CONSTANTS 
##
####################################################
const DelaunayTriangulation = DelaunayTriangulation
const DT = DelaunayTriangulation

export DelaunayTriangulation

const DefaultAdjacentValue = 0
const BoundaryIndex = -1
const FirstPointIndex = DefaultAdjacentValue + 1

####################################################
##
## DEPENDENCIES
##
####################################################
using DataStructures
using SimpleGraphs
using ExactPredicates
using EnumX
using MakieCore
using Random

####################################################
##
## FILES AND EXPORTS
##
####################################################
include("interfaces/points.jl")
include("interfaces/triangles.jl")
include("interfaces/edges.jl")
include("interfaces/boundary_nodes.jl")

export indices
export num_triangles
export each_triangle
export geti
export getj
export getk
export initial
export terminal
export num_edges
export each_edge
export getx
export gety
export getxy
export each_point
export each_point_index
export get_point
export num_points
export get_boundary_nodes
export num_boundary_edges

include("data_structures/adjacent.jl")
include("data_structures/adjacent2vertex.jl")
include("data_structures/graph.jl")
include("data_structures/convex_hull.jl")
include("data_structures/triangulation.jl")
include("data_structures/representative.jl")

export get_adjacent
export get_adjacent2vertex
export get_graph
export get_edges
export get_neighbours
export get_points
export get_triangles
export get_boundary_map
export get_constrained_edges
export get_boundary_nodes
export get_convex_hull
export get_boundary_index_ranges
export Triangulation
export ConvexHull
export convex_hull
export convex_hull!
export each_solid_triangle
export each_ghost_triangle
export get_vertices
export clear_empty_features!
export get_indices
export get_convex_hull_indices
export each_vertex
export num_vertices
export each_solid_edge 
export each_ghost_edge 
export each_solid_vertex 
export each_ghost_vertex

include("predicates/certificate.jl")
include("predicates/boundaries_and_ghosts.jl")
include("predicates/general.jl")
include("predicates/index_and_ghost_handling.jl")

export Certificate

include("operations/add_triangle.jl")
include("operations/add_boundary_information.jl")
include("operations/add_ghost_triangles.jl")
include("operations/delete_triangle.jl")
include("operations/delete_ghost_triangles.jl")
include("operations/add_point.jl")
include("operations/flip_edge.jl")
include("operations/split_edge.jl")
include("operations/split_triangle.jl")
include("operations/legalise_edge.jl")
include("operations/delete_point.jl")

export add_ghost_triangles!
export delete_ghost_triangles!
export add_point!
export add_triangle! 
export delete_triangle! 
export flip_edge!
export add_boundary_information!
export split_edge!
export split_triangle!
export legalise_edge!
export delete_point!

include("triangulation/gmsh.jl")
include("triangulation/rectangle.jl")
include("triangulation/bowyer_watson.jl")
include("triangulation/triangulate.jl")
include("triangulation/convex_triangulation.jl")

export generate_mesh
export triangulate_rectangle
export triangulate
export triangulate_convex

include("point_location/brute_force.jl")
include("point_location/initialisers.jl")
include("point_location/jump_and_march.jl")

export brute_force_search
export jump_and_march

include("plotting.jl")

export triplot
export triplot!

include("utils.jl")
include("polygon_utils.jl")

export number_type
export polygon_features
export pole_of_inaccessibility
const polylabel = pole_of_inaccessibility
export polylabel

####################################################
##
## EXTRA DOCSTRINGS 
##
####################################################
"""
    Interfaces 

The main objects used throughout this package, namely 
points, edges, triangles, and boundaries, have interfaces 
that are fully customisable. We describe these interfaces 
below, along with all the methods that need to be 
defined. 

For these interfaces, we define default methods that 
support `Tuple`s, `AbstractVector`s, etc., so you
do not need to extend any of these for normal use. Moreover, 
while we do allow for different index types to be using, 
all still have to subtypes of `Integer` (for now?).

# Triangle Interface 

Here we describe the interface used for describing 
triangles and collections of triangles. 

## Individual Triangles 

A triangle is assumed to be of the form `(i, j, k)`,
but the object used for storing these three indices 
can be customised. The following methods are used for  
working with triangles; see their docstrings for more 
information.

- [`construct_triangle`](@ref): Must be defined.
- [`geti`](@ref): Must be defined.
- [`getj`](@ref): Must be defined.
- [`getk`](@ref): Must be defined.
- [`indices`](@ref): Uses the [`geti`](@ref), [`getj`](@ref), and [`getk`](@ref) definitions.
- [`integer_type`](@ref): Must be defined.
- [`triangle_edges`](@ref): Uses the [`geti`](@ref), [`getj`](@ref), and [`getk`](@ref) definitions.
- [`rotate_triangle`](@ref): Uses the [`geti`](@ref), [`getj`](@ref), and [`getk`](@ref) definitions.
- [`construct_positively_oriented_triangle`](@ref): Uses existing definitions.
- [`compare_triangles`](@ref): Uses the [`geti`](@ref), [`getj`](@ref), and [`getk`](@ref) definitions.
- [`sort_triangle`](@ref): Uses the [`indices`](@ref) and [`construct_triangle`](@ref) definitions.

## Collection of Triangles 

A collection of triangles simply stores many triangles. The 
following methods are used for working with these collections; see
their docstrings for more information. Note that these collections 
must be mutable.

- [`initialise_triangles`](@ref): Must be defined.
- [`triangle_type`](@ref): Must be defined. 
- [`num_triangles`](@ref): Must be defined.
- [`contains_triangle`](@ref): This calls `Base.in` and [`contains_triangle`](@ref).
- [`add_to_triangles!`](@ref): Must be defined. 
- [`add_triangle!`](@ref): This simply calls [`add_to_triangles!`](@ref).
- [`delete_from_triangles!`](@ref): Must be defined. 
- [`delete_triangle!`](@ref): This simply calls [`delete_from_triangles!`](@ref).
- [`each_triangle`](@ref): Must be defined. 
- `Base.in`: Must be defined.
- `Base.sizehint!`: Must be defined.
- `Base.unique!`: Must be defined, unless your collection is a `Set`.
- [`compare_triangle_collections`](@ref): Calls [`num_triangles`](@ref), [`each_triangle`](@ref), and [`contains_triangle`](@ref).
- [`sort_triangles`](@ref): Uses the [`initialise_triangles`](@ref), [`each_triangle`](@ref), [`sort_triangle`](@ref), and [`add_triangle!`](@ref) definitions.
- [`remove_duplicate_triangles`](@ref): Uses `Base.unique!`.

Note that [`Triangulation`](@ref)s also define [`each_solid_triangle`](@ref) and [`each_ghost_triangle`](@ref).
# Edge Interface 

Here we describe the interface used for describing 
edges and collections of edges. 

## Individual Edges 

An edge is assumed to be of the form `(i, j)`,
but the object used for storing these two indices 
can be customised. The following methods are used for  
working with edges; see their docstrings for more 
information.

- [`construct_edge`](@ref): Must be defined. 
- [`initial`](@ref): Must be defined. 
- [`terminal`](@ref): Must be defined. 
- [`edge_indices`](@ref): Uses the [`initial`](@ref) and [`terminal`](@ref) definitions.

## Collection of Edges 

A collection of edges simply stores many edges. The 
following methods are used for working with these collections; see
their docstrings for more information. Note that these collections 
must be mutable.

- [`initialise_edges`](@ref): Must be defined. 
- [`edge_type`](@ref): Must be defined. 
- [`num_edges`](@ref): Must be defined.
- [`contains_edge`](@ref): Must be defined.
- [`add_to_edges!`](@ref): Must be defined. 
- [`add_edge!`](@ref): This simply calls [`add_to_edges!`](@ref).
- [`delete_from_edges!`](@ref): Must be defined. 
- [`delete_edge!`](@ref): This simply calls [`delete_from_edges!`](@ref).
- [`each_edge`](@ref): Must be defined. 
- [`random_edge`](@ref): Uses the [`each_edge`](@ref) definition.
- [`is_empty`](@ref): Simply uses the `isempty` definition from Base.

# Point Interface 

Here we describe the interface used for describing 
points and collections of points.

## Individual Points 

A point is assumed to be of the form `(x, y)`,
but the object used for storing these two coordinates 
can be customised. The following methods are used for  
working with points; see their docstrings for more 
information.

- [`getx`](@ref): Must be defined.
- [`gety`](@ref): Must be defined.
- [`getxy`](@ref): Uses the [`getx`](@ref) and [`gety`](@ref) definitions.

## Collection of Points

A collection of points simply stores many points. The 
following methods are used for working with these collections; see
their docstrings for more information.

- [`getpoint`](@ref): Must be defined.
- [`get_point`](@ref): This simply calls [`getpoint`](@ref).
- [`each_point_index`](@ref): Must be defined. 
- [`each_point`](@ref): Must be defined.
- [`num_points`](@ref): Must be defined.
- [`points_are_unique`](@ref): Makes use of the existing methods.
- [`lexicographic_order`](@ref): Makes use of the existing methods.
- [`number_type`](@ref): Must be defined.

# Boundary Nodes Interface 

Here we describe the interface used for representing the boundaries. As 
described in [`Triangulation`](@ref), this interface can be used for 
representing either a collection of curves each made up of multiple 
segments, a collection of segments, or a single continuous curve. 
The following functions facilitate these possibilities; see the
corresponding docstrings for more information.

- [`has_multiple_curves`](@ref): Must be defined. 
- [`has_multiple_segments`](@ref): Must be defined.
- [`num_curves`](@ref): Must be defined. 
- [`num_segments`](@ref): Must be defined.
- [`num_boundary_edges`](@ref): Must be defined. 
- [`getboundarynodes`](@ref): Must be defined. 
- [`get_boundary_nodes`](@ref): This just calls [`getboundarynodes`](@ref).
- [`each_boundary_node`](@ref): Must be defined.
- [`construct_boundary_map`](@ref): Makes use of the existing methods.
- [`construct_boundary_index_ranges`](@ref): Makes use of the existing methods.
- [`map_boundary_index`](@ref): Makes use of the result from [`construct_boundary_map`](@ref).
- [`get_curve_index`](@ref): Makes use of the result from [`construct_boundary_map`](@ref).
- [`get_segment_index`](@ref): Makes use of the result from [`construct_boundary_map`](@ref).
- [`num_outer_boundary_segments`](@ref): Makes use of the existing methods.
"""
function Interfaces end

export Interfaces

end
