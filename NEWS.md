# Changelog

## 1.6.0
- Feature: Define `reverse` for `AbstractParametricCurve`s, making it easier to reverse the orientation of a curve. See [#195](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/195).
- Bugfix: Fixed an issue with `LineSegment` not returning the exact endpoints at `t=1`, which can be problematic when joining boundary nodes. This has been fixed. See [#195](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/195).
- Bugfix: Introduced `is_linear` to fix issues with boundary enrichment of domains with `LineSegment`s. In particular, `LineSegment`s are no longer enriched. See [#195](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/195).
- Bugfix: `orientation_markers` now uses `uniquetol` instead of `unique` for the final set of markers (it already did it for the intermediate markers). See [#195](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/195).
- Bugfix: For large `Tuple`s, functions like `eval_fnc_at_het_tuple_two_elements` are problematic and allocate more than their non-type-stable counterparts. To get around this, for `Tuple`s of length `N > 32`, the non-type-stable version is used. See [#195](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/195).
- Bigfix: Fixed issue with `use_barriers` when a ghost edge is selected at random during point location. See [#196](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/196).
- Feature: Introduced the (currently internal) function `get_positive_curve_indices` for finding curves with positive orientation in a `Triangulation`. [#196](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/196).
- `is_exterior_curve`, `is_interior_curve`, `num_exterior_curves`, and `is_disjoint` are now defined based on `get_positive_curve_indices` rather than `get_exterior_curve_indices`. See [#196](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/196).
- Bugfix: `PointLocationHistory` was not marked as public. This has been fixed. 

## 1.5.0

- Introduced the ability to reconstruct unconstrained triangulations from an existing set of points and triangles using `Triangulation(points, triangles)`. See [#192](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/192)

## 1.4.0

- Updated to AdaptivePredicates.jl v1.2, now allowing caches to be passed to the predicates involving `incircle` and `orient3`. These are only useful when using the `AdaptiveKernel()` kernel. Outside of triangulating, these caches are not passed by default, but can be provided. The functions `get_incircle_cache` and `get_orient3_cache` can be used for this purpose on a triangulation (without a triangulation, refer to AdaptivePredicate.jl's `incircleadapt_cache` and `orient3adapt_cache`). See [#185](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/185).

## 1.3.1

- Fix an issue with a weighted triangulation where the lifted points' convex hull was entirely coplanar. See [#184](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/184)

## 1.3.0

This release finally introduces weighted triangulations and power diagrams, and also allows for users to provide a generic convex polygon to for clipping a Voronoi tessellation instead of only the convex hull.

- Weighted triangulations have now been implemented, as have power diagrams. The weights are also no longer restricted to `Float64` type. See [#180](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/180).
- `intersection_of_edge_and_bisector_ray` now accepts a `project` keyword argument. See [#180](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/180).
- `get_weight(w, i)` now returns, when `i` is not an integer, either `i[3]` if it represents a point in space or `0`. See [#180](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/180).
- Define `project_onto_line(p, q, r)` for projecting a point `r` onto the line defined by `p` and `q`. See [#180](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/180).
- Fixed a bug with clipping Voronoi tessellations in cases where there are no intersections of any Voronoi polygon with the convex hull. See [#180](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/180).
- `voronoi` now accepts an optional `clip_polygon` keyword argument, defaulting to `nothing` (corresponding to the convex hull), allowing for a convex clip polygon to be used instead of the convex hull. The `clip_polygon` should be a `Tuple` of the form `(points, boundary_nodes)` where the `boundary_nodes` give vertices of `points` adhering to the usual convention. Note that this could be used as an alternative to looping over `get_polygon_coordinates` for clipping to a rectangle. See [#180](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/180).
- `centroidal_smooth` now accepts `clip_points` and `clip_vertices` as keyword arguments, defaulting to `nothing` (corresponding to the convex hull), to accommodate the new `clip_polygon` keyword argument in `voronoi`. See [#180](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/180).
- `has_multiple_curves`, `has_multiple_sections`, and `num_boundary_edges` now have methods for `Tuple`s of integers. A bug was also fixed with `number_type` of a `Tuple` of `Tuple`s of coordinates returning the `Tuple` type instead of the coordinate type. See [#180](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/180).

## v1.2.0

- Warnings are now thrown when you try and triangulate point sets not in the plane. The `is_planar` function has been introduced for this. See [#178](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/178).

## v1.1.4

- Fixed a bug with curve-bounded refinement with custom edge(s) structs. See [#175](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/175).

## v1.1.3

- `sort_triangle` is now public. See [#174](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/174).

## v1.1.2

- Clarified type stability of `triangulate` in docstring, and notes about field access and public API. See [#171](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/171).

## v1.1.1

- Fixed issue on nightly with symbols being marked as both public and exported. See [#168](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/168).

## v.1.1.0

There are a lot of changes in this release, most of them irrelevant for the user. The most important change is the following:

- We now support a choice between fast, exact, and adaptive predicates via `FastKernel()`, `ExactKernel()`, and `AdaptiveKernel()`, respectively. The default is now `AdaptiveKernel()`. Moreover, triangle areas are now computed using the adaptive `orient` predicate to be more robust. See [#165](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/165).

Previously, ExactPredicates.jl was used everywhere, which can be slow and not necessary for certain point sets. The `FastKernel()` option 
has no exact or adaptive arithmetic and so should be used with caution. The documentation discusses these choices in more detail. 

To actually configure the choice of predicate, you can e.g. in `triangulate` use the `predicates` keyword argument and pass one of 
`DelaunayTriangulation.FastKernel()`, `DelaunayTriangulation.ExactKernel()`, or `DelaunayTriangulation.AdaptiveKernel()`. If you are computing a predicate manually, then the predicate is instead passed as the first argument.

Some other changes: 

- Added `DelauanyTriangulation.validate_triangulation` for validating triangulations. See [#131](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/131).
- Fixed a bug with the currently unused `orient(p, q, r, s)` predicate. See [#131](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/131).
- Added private functions `getz`, `_getz`, `getxyz`, and `_getxyz`. See [#131](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/131).
- `jump_and_march` has now been renamed to `find_triangle`. For compatibility, `jump_and_march` still works and is simply an alias of `find_triangle`. See [#133](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/133).
- Mutable structs now use `const` on fields that aren't changed. For compatibility with older versions, this is implemented using a macro that is a no-op where this is not supported. See [#140](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/140).
- We now use the `public` word to define public functions. This is only included on Julia versions v1.11 and above. See [#140](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/140).
- We now test on the pre-release. See [#140](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/140).
- The module `DelaunayTriangulation` now has a docstring. See [#140](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/140).
- The `.md` files for tutorials and applications in the docs have been properly updated to match their literate counterparts. See [#140](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/140).
- We now use a workflow to enforce changes to `NEWS.md` for any PRs. See [#140](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/140).
- Improved the error message for an incorrect orientation. See [#144](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/144).
- Added a CONTRIBUTING.md file and issue templates. See [#160](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/160).
- Added `is_point2` and `is_point3` to detect if a given input is a point. This allows vector coordinates to be passed to `convert_boundary_points_to_indices`. See [#161](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/161).
- Removed an allocation from `add_vertex!`. See [#163](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/163).
- Fixed an issue with the user-supplied `rng` not being passed to `lock_convex_hull!`.

## v1.0.5 

- Disabled `deepcopy` on `PolygonTree`s and made it a no-op. See [#129](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/issues/129).

## v1.0.4

Nothing breaking. Main changes:

- Fixes some issue with type instabilities
- Adds `construct_polygon_hierarchy` to the public API list (it was already intended to be there)
- All computations are now in the provided precision except for the triangle area and circumcenter which are done in `Float64`. A warning is given when a non-Float64 precision is given.
- Closes #118 
- Redesigns the `Polygon` struct to have an `is_circular` field, avoiding the need for `views`
- Completely refactors `validate_triangulation`. Still only lives in the test files though. Maybe one day it can live inside the package itself incase users somehow have a use for it..
- Clean up the runtests.jl file, and make sure that the README/docs are fully tested
- Use Aqua and fix ambiguities
- Remove accidental piracy of `minimum(::Nothing)`
- Remove CI spam from method redefinitions
- Fix some issues with doc images and some typos
- Fix #109

## v1.0.3

- Removed some old function definitions that were no longer needed anymore following the new Makie release.

## v1.0.0

In addition to the changes below, note that many bugs have been fixed. Feel free to make any issues or PRs if you encounter any problems.

### Added
- Triangulation and refinement of curve-bounded domains has now been added.
- The `Triangulation` struct has some more fields. One of these is `weights`, although this is not to be used just yet; it is possible to computed weighted Delaunay triangulations of convex polygons, but not for general domains currently.
- Where reasonable, default methods for the geometric primitive interfaces have now been added. For example, `getx` now has the default definition `getx(p) = p[1]`. 
- You can now use `each_boundary_edge` to iterate over the boundary edges (in no specific order), rather than having to use `keys(get_boundary_edge_map(tri))` as you had to before 1.0.
- The documentation has had a complete overhaul.
- The `retriangulate` function has now been added, allowing for the retriangulation of a given `Triangulation`.
- The function `dist` can now be used for computing the distance between a point and the triangulation, rather than having to use `distance_to_polygon`.
- You can now use `find_polygon` to find what polygon in a provided boundary of a `Triangulation` contains a given point.
- Mesh refinement now, by default, defaults to the diametral lens definition of encroachment rather than the diametral circle definition. You can toggle this behaviour using the `use_lens` keyword in `refine!`.
- Now whenever encroached edges are split during mesh refinement, all free vertices inside that edge's diametral circle will be deleted prior to splitting.

### Breaking
Please note that a lot of these changes that follow are only _technically_ breaking, because I failed to properly specify what the public API was for this package (this has now been corrected in v1.0 and onwards).

- All references to constrained edges are now referred to as _segments_. For example, the field `constrained_edges` is now `interior_segments`, and `all_constrained_edges` is now `all_segments`. The same goes for the associated getters.
- The keyword `edges` in `triangulate` is now `segments`.
- `lock_convex_hull!` will now delete from `interior_segments` field any edges it encounters along the convex hull during locking.
- `add_edge!` is now `add_segment!`.
- `triangulate_convex` no longer accepts the `recompute_centers` keyword.
- `RefinementTargets` is now `RefinementConstraints`.
- `maxiters` is no longer a valid keyword argument of `refine!`. Instead, you should just pass `max_points` appropriately. With this change, the default for `max_points` is no longer `typemax(Int)` but is instead `max(1000, num_solid_vertices(tri))^2`.
- `refine!` no longer returns `statistics(tri)`. It now only returns `tri`. If you want the statistics, just use `statistics(tri)` afterwards.
- `refine!` no longer accepts `max_radius_edge_ratio`. Instead, you must provide `min_angle`. (Internally, it is still `max_radius_edge_ratio` that gets primarily used, but `min_angle` is more interpretable for a user.)
- `lock_convex_hull` is no longer a keyword of argument of `refine!` - if it has to be locked, it will be automatically.
- `edge_indices` has been removed and is now `edge_vertices`.
- `indices` has been removed and is now `triangle_vertices`.
- `initialise_edges` and `initialise_triangles` have been removed. Just use the types instead, e.g. `initialise_edges(::Type{S})` is now just `S()`.
- `is_empty` has been removed. Just use `isempty`.
- `compare_unoriented_edge` is now `compare_unoriented_edges`.
- `initial` and `terminal` are no longer exported. `edge_vertices` is more appropriate for this.
- `geti`, `getj`, `getk` are no longer exported. `triangle_vertices` is more appropriate for this.
- `each_point` and `each_point_index` are no longer exported. You are encougared to rely on `each_solid_vertex` instead. Just note that `each_solid_vertex` is not an ordered iterator.
- `num_points` is no longer exported. You are encouraged to instead use `num_solid_vertices`.
- `remove_duplicate_triangles` has been removed. Just use `sort_triangles` and `unique!`.
- I used to refer to distinct parts of a boundary as `segments`. This is misleading, since it's the first name that I use for constrained edges. So, I now use the term `sections`. For example, `has_multiple_segments` is now `has_multiple_sections`.
- `has_multiple_sections` can no longer be used on the ghost vertex map.
- `getboundarynodes` has been removed. Just use and extend `get_boundary_nodes` instead.
- `BoundaryIndex` is now `GhostVertex`. Similar references to `boundary_index` have been changed to `ghost_vertex`.
- `boundary_map` is now `ghost_vertex_map`.
- `boundary_index_ranges` is now `ghost_vertex_ranges`.
- `integer_type` is removed (exceptfor `Triangulation`s). `number_type` should be used for this.
- Previously, `get_boundary_index` (now called `get_ghost_vertex`) errored if neither of the three arguments was a ghost vertex. Now, if none of the arguments is a ghost vertex, the function returns `k`. Similarly for the two-argument version.
- `rotate_ghost_triangle_to_standard_form` and `rotate_triangle_to_standard_form` are removed. Use `sort_triangle` instead.
- `num_outer_boundary_segments` has been removed.
- Iteration over `Adjacent` and `Adjacent2Vertex` has been removed.
- `Adjacent2Vertex` is now defined by `Adjacent2Vertex{IntType,EdgesType}` rather than `Adjacent2Vertex{IntType,EdgesType,EdgeType}`.
- `clear_empty_points!` is now `clear_empty_vertices!`.
- The fields of `ConvexHull` are now `points` and `vertices` rather than `points` and `indices`. Additionally, the type is now `ConvexHull{PointsType, IntegerType}` instead of `ConvexHull{PointsType, Vector{IntegerType}}`.
- `num_points` is no longer defined on `ConvexHull`s.
- The `boundary_map` (now `ghost_vertex_map`) and `boundary_index_ranges` (now `ghost_vertex_ranges`) are now `Dict`s instead of `OrderedDict`s.
- The constructor `Triangulation(points, triangles, boundary_nodes; kwargs...)` now uses the keyword argument `delete_ghosts` instead of `add_ghost_triangles`, with default `delete_ghosts=false` which is opposite to the previous default `add_ghost_triagnles=false`.
- The function `get_empty_representative_points` has been deleted.
- `get_convex_hull_indices` is now `get_convex_hull_vertices`.
- `min_max` is removed - just use `Base.minmax`.
- `nearest_power_of_two` has been removed.
- `intersection_of_ray_with_boundary` has been removed.
- `identify_side` has been removed.
- `intersection_of_ray_with_edge` has been removed.
- The field names of `TriangulationStatistics` have been renamed to match the new segment terminology, and the field `num_convex_hull_points` is now `num_convex_hull_vertices`.
- `get_total_area` is now `get_area`. The field name `total_area` in `TriangulationStatistics` is now `area`.
- `point_position_relative_to_box` has been removed.
- `is_boundary_edge` has had its definition changed: `is_boundary_edge(tri, i, j)` is now `is_boundary_edge(tri, j, i)`, so that we have consistency with `is_boundary_triangle`.
- `is_outer_boundary_index` is now `is_exterior_ghost_vertex`.
- `is_outer_ghost_triangle` is now `is_exterior_ghost_triangle`.
- `is_outer_ghost_edge` is now `is_exterior_ghost_edge`.
- `is_outer_boundary_node` is now `is_exterior_ghost_vertex`.
- The keyword arguments in `add_ghost_triangles!` and `delete_ghost_triangles!` have been removed.
- The keyword argument `exterior_curve_index` in `triangulate` has been removed. The exterior curves will be automatically computed.
- `peek_triangle_ρ` has been removed, and `triangle_dequeue!` now does a `dequeue_pair!`.
- The keyword `recompute_representative_point` is now `recompute_representative_points` in `triangulate`.
- The keyword `add_ghost_triangles` in `triangulate_rectangle` (and `triangulate_convex`) is now `delete_ghosts` to be consistent with `triangulate`. The default is `delete_ghosts=false`, consistent with the previous default of `add_ghost_triangles=true`.
- `voronoi` now has the signature `voronoi(tri; clip=false, smooth=false, kwargs...)` instead of `voronoi(tri, clip=false)` (note the difference in `clip` as a positional argument to a keyword argument). `centroidal_smooth` is still exported (the `kwargs...` are passed to it), but this can be a good alternative.
- Gmsh support has been removed.
- `get!` is no longer defined on `adjacent2vertex`.
- `Certificate` is no longer exported.
- The arguments to `construct_ghost_vertex_map` (previously `construct_boundary_map`), `construct_boundary_edge_map`, and `construct_ghost_vertex_ranges` (previously `construct_boundary_index_ranges`) are now positional.
- `get_vertices` is no longer exported.
- `map_boundary_index` is now `map_ghost_vertex`.
- `polylabel` is no longer an alias for `pole_of_inaccessibility`.
- `convert_to_boundary_edge` is now `convert_to_edge_adjoining_ghost_vertex`.
- References to a `point_order` are now an `insertion_order`.
- `add_boundary_information!` is no longer exported.
- `balanced_power_of_two_quarternary_split` has been removed.
- To provide a custom constraint for triangulation, you now need to provide a constraint of the form `(tri::Triangulation, T::Triangle) -> Bool` with the keyword argument `custom_constraint`, where the returned result is `true` if the triangle `T` should be refined and `false` otherwise.
- `is_triangle_seditious` no longer has the unused argument `ρ` at the end of its signature.
- It is no longer possible to use `delete_point!` on points that are on the boundary or adjoin segments. It is too finicky, and for constrained triangulations you shouldn't be doing those operations anyway. For unconstrained triangulations, deleting a point on the boundary might be reasonable but it is difficult and `delete_point!` will not update the convex hull field of the triangulation for you either. The previous implementation did work sometimes, but it required playing around with how ghost triangles are interpreted and was unfortunately not reliable. PRs to address this are welcome, but it would need some careful thought and highly extensive testing (e.g. thousands of random tests for random triangulations, especially triangulations with very skinny or very large triangles on the boundary; triangulations with many collinearities such as from `triangulate_rectangle`; triangulations with holes; and triangulations with disjoint domains).
- `is_positively_oriented` now only accepts integers for the second argument, to avoid issues with misleading results when testing triangles.
- Some methods like `get_right_boundary_node` have now been defined only for `Triangulation`s, instead of lowering to a more complicated method that uses the `adjacent` field and other fields like we used to.
