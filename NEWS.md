# Changelog

## 1.7.1

- Fixed a clipped Voronoi tessellation leaving a polygon unclipped when the polygon only crosses the convex hull through vertices that lie exactly on the hull. Circumcenters that lie exactly on the hull are also no longer duplicated when clipping. See [#223](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/issues/223).
- Fixed neighbouring polygons of a tessellation clipped to a `clip_polygon` getting slightly different copies of a shared vertex, which meant `get_adjacent` could not find the neighbour. The coordinates of an unbounded polygon from `get_polygon_coordinates` now include a point on each unbounded edge that only depends on that edge. Clipped polygons also no longer repeat a vertex when an edge passes through a corner of the clip polygon. See [#239](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/issues/239).
- `triangulate` now throws a `CollinearPointsError` when there are more than three points and they are all collinear, rather than an `UndefVarError` or hanging. Nearly collinear point sets, where every triangle is degenerate or very nearly degenerate, are now triangulated correctly. See [#216](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/issues/216).
- `pole_of_inaccessibility` no longer runs out of memory for very thin polygons, such as the convex hull of nearly collinear points. See [#216](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/issues/216).
- `convert_boundary_points_to_indices` now skips points that repeat the previous point, e.g. a repeated final point. `triangulate` now throws a `DuplicateBoundaryNodeError` if a boundary node is a duplicate point that would be skipped, rather than overflowing the stack. See [#220](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/issues/220).
- `copy` of a `Triangulation` or a `VoronoiTessellation` now returns a copy that is fully independent of the original. Previously, the copy shared some mutable containers with the original, e.g. the neighbours in the graph, the sections of the boundary nodes, the boundary enricher of a curve-bounded domain, and the polygons of a tessellation, so mutating the copy (e.g. with `add_point!`, `refine!`, or by clipping a tessellation) could change the original.
- `refine!` no longer depends on the scale of the domain. The tolerances used to detect precision issues were absolute, so for small domains (e.g. a unit square scaled by `1e-3`) refinement stopped early with far fewer vertices than for the same domain at a larger scale. For domains smaller than `1`, these tolerances are now relative to the size of the domain, and for larger domains they are unchanged. Also fixed `refine!` calling `pop!` on a `Set`, which removed an arbitrary vertex or threw an error, when a subsegment was too small to split.
- `triangulate` no longer overflows the stack for many points at a small scale, e.g. `1e-5 .* rand(2, 100_000)`. Whether a candidate initial triangle is nearly degenerate is now decided relative to the size of the triangle rather than by an absolute tolerance on its area.
- The default `precision` of `pole_of_inaccessibility` and `compute_representative_points!` is now `max(width, height) / 1000`, where `width` and `height` are the dimensions of the polygon's bounding box, instead of `1`. Previously, for polygons that are not much larger than `1`, the returned pole was just the initial guess. Distinct cells with the same priority were also treated as equal, so some were dropped from the search. The representative points of a triangulation may therefore differ slightly from those of previous versions.
- The representative point of the first boundary curve of a domain with multiple curves was computed using every curve, so for a domain with several disjoint parts it could lie inside another part. It now only uses the first curve and its holes. Previously, point location outside of such a domain (e.g. with `find_triangle`) could then fail to terminate.
- Points with coordinates `-0.0` and `0.0` are now treated as the same point when looking for duplicate points, e.g. in `triangulate` and `convex_hull`.
- Fixed `convex_hull` throwing a `BoundsError` when the points have duplicates and fewer than four unique points, and throwing an error when there are no points.
- Fixed `get_section_index(tri, ghost_vertex)` throwing an `UndefVarError`, and `is_interior_ghost_vertex(tri, i)` returning `true` for solid vertices `i`.
- Fixed `find_polygon` returning `0` for points inside a domain whose single boundary curve is split into sections, e.g. the default boundary from `triangulate_rectangle`.
- Fixed `get_all_stat` throwing an error for a triangulation without any triangles. The getters for the statistics of an individual triangle, e.g. `get_area(stats, T)`, now throw a `KeyError` rather than a `BoundsError` when `T` is not a triangle of the triangulation.
- `triangulate` now checks more of its arguments (unless `check_arguments = false`). It throws an `IntegerCoordinatesError` for integer coordinates, a `NonFiniteCoordinatesError` for a `NaN` or infinite coordinate, an `InconsistentWeightsError` if a vector of weights has fewer weights than there are points, an `InvalidSegmentError` for a segment `(i, i)` or a segment with a vertex that is not the index of a point, and an `InvalidBoundaryNodeError` for a boundary node that is not the index of a point. Previously, these gave errors such as an `InexactError`, `KeyError`, or `BoundsError`, and a segment `(i, i)` was either ignored or gave an `AssertionError` depending on the points.
- Invalid arguments to `BSpline`, `CatmullRomSpline`, `convert_boundary_points_to_indices`, `get_polygon_coordinates` (an invalid bounding box) and `voronoi` (`smooth = true` without `clip = true`) now throw an `ArgumentError` rather than an `AssertionError`. `triangulate_convex` now throws an `ArgumentError` when no three consecutive vertices of `S` are positively oriented (e.g. if they are all collinear or `S` is clockwise), rather than overflowing the stack. `voronoi` and `centroidal_smooth` now undo their temporary changes to the triangulation even if an error is thrown. `InconsistentConnectionError` messages now refer to sections of a boundary rather than segments.
- Triangulations use about 40% less memory. The two triangulations in the cache of a `Triangulation`, which are only used for small cavities, are no longer sized for all of the points. For 100,000 random points, a triangulation went from 233 MB to 143 MB.
- `validate_triangulation` is much faster for large triangulations, e.g. about 10 times faster for 8,000 vertices. Comparing the triangle iterators with the triangles took quadratic time.
- `validate_triangulation` no longer changes the triangulation. It rebuilt the ghost triangles and cleared the empty features of the triangulation itself, which changed the order of iteration over its triangles and edges, so that validating a triangulation could change the results of later operations. It now checks a copy instead. Its result also no longer depends on the global random number generator.
- `triangulate` is much faster for domains bounded by Bézier curves, B-splines or Catmull-Rom splines, e.g. about 1 ms instead of 84 ms for a single Bézier curve. To construct the polygon hierarchy, the curves were polygonised by repeatedly splitting them into pieces of equal total variation, which needs numerical integration for each split. Their points are now evenly spaced in the curve's parameter. The triangulation itself is unchanged.
- Fixed `voronoi(tri; clip = true, smooth = true)` throwing a `KeyError` for most refined triangulations `tri`. Smoothing locks the convex hull, and some of the points that `refine!` adds on the edges of the convex hull end up just inside of it. Constrained triangulations kept the thin triangles, outside of the domain, between such boundary nodes and the edge of the convex hull over them, since the computed distance from their centroid to the boundary did not have the right sign. Triangles whose vertices are all on the boundary are now classified as inside or outside of the domain using the boundary edges instead of their coordinates.
- `triangulate` now throws a `WeightedConstrainedTriangulationError` when weights are given together with segments or boundary nodes (unless `check_arguments = false`), since only unconstrained weighted triangulations are supported. Previously this was not rejected, and e.g. a segment with a submerged vertex gave a `KeyError`. `RefinementConstraints`, and so `refine!`, now throws an `ArgumentError` for a `min_angle` outside of `[0, 60]`, for which refinement could never finish since the smallest angle of a triangle is at most 60°, and for a `max_area` that is not positive.
- Clipping a Voronoi tessellation (`voronoi(tri; clip = true)`) no longer takes quadratic time when many of the points are on the convex hull, e.g. 45 ms instead of 260 ms for 8,000 points on a circle. Each intersection of a polygon with the boundary was looked for in the list of all such intersections.
- `refine!` no longer loops forever for curve-bounded domains smaller than about `1e-7`. Splitting the boundary curves still used an absolute tolerance to detect precision issues, so that every split was rejected. It now uses the same scale as the rest of the refinement, as does the enrichment of the boundary in `triangulate`.
- `triangulate` no longer overflows the stack when there are very many points outside of the domain, e.g. 300,000 random points around a small boundary. Finding the points to delete outside of the domain used recursion whose depth grew with the number of these points; it now uses a stack instead.
- Point location for a point outside of the triangulation (e.g. in `find_triangle`) could loop forever when the boundary is not convex as seen from its representative point, e.g. for a representative point outside of the domain. The search along the boundary now goes around it at most once, after which the search is restarted, eventually falling back to a brute-force search. `DelaunayTriangulation.exterior_find_triangle` now also returns whether it found the ghost triangle.
- `retriangulate` now copies the segments and boundary nodes it takes from the triangulation. Previously they were shared, so that, for example, the triangulation of a tessellation smoothed with `voronoi(tri; clip = true, smooth = true)` or `centroidal_smooth` lost its boundary nodes when the convex hull of the input triangulation was unlocked again, leaving it in an invalid state.
- The return types of `triangulate` and `triangulate_rectangle` are inferred again on Julia 1.14. Their types were passed to the `Triangulation` constructor as keyword arguments, which Julia 1.14.0-alpha1 no longer infers through.
- `Triangulation` and `VoronoiTessellation` now have a short two-argument `show`, e.g. `Triangulation with 4 vertices and 2 triangles`. Previously, printing one inside e.g. a `Vector` or a failed `@test` printed every field.
- The undocumented keyword arguments `full_polygon_hierarchy` and `boundary_curves` of `triangulate`, which were only used internally for curve-bounded domains, have been removed. They are now taken from `boundary_enricher`, which is documented as internal.
- The functions documented on the Primitive Interfaces page of the API reference (e.g. `geti`, `construct_edge`, `push_point!`, and `contains_edge`), which are the methods to overload when using custom primitives, are now marked `public`.
- Fixed several errors in the documentation. In particular, the mathematical details now use the correct definition of a Delaunay edge, the correct sign convention for the orientation predicate, the correct test for a point being inside a diametral circle or lens, and the correct second derivatives of circular and elliptical arcs, and they describe how unbounded Voronoi polygons are clipped since [#239](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/issues/239). Some docstrings that were attached to the wrong function (e.g. for `each_generator`) have also been fixed.
- `triangulate` is much faster for large numbers of points, e.g. about twice as fast for 100,000 random points, three times as fast for 300,000, and 4.5 times as fast for 1,000,000 (14 s instead of 66 s). With `randomise = true`, the points are now inserted in a biased randomised insertion order (see `DelaunayTriangulation.brio_order`) rather than in a uniformly random order. The order is still random, but within each of its rounds the points are sorted along a Hilbert curve, so that consecutive points are nearly always close together. The point location for each new point then starts from the closer of the last inserted point and one random vertex (the new default `num_sample_rule` with `randomise = true` is `DelaunayTriangulation.spatial_num_samples`), instead of the closest of `∛n` random vertices. For points in general position the triangulation is the same, but for degenerate point sets, e.g. points on a grid, `triangulate` may now return a different (still Delaunay) triangulation for the same `rng`, and the `rng` is used differently. To insert the points in the previous order, use `insertion_order = setdiff!(shuffle!(rng, collect(each_point_index(points))), skip_points)` and `num_sample_rule = DelaunayTriangulation.default_num_samples`.
- Fixed `refine!` throwing a `BoundsError` for some domains with narrow channels between boundaries, e.g. for about a quarter of the random point sets in the tight walls example in the tests. When the circumcenter of a triangle was behind the boundary at the vertex that its point location started from, the point location stopped before it found any triangle, and the empty triangle it returned was used. The centroid of the triangle is now inserted instead, as for other circumcenters that can not be seen from their triangle.
- `triangulate` is much faster for domains with many boundary nodes. To find the points outside of the domain, the distance from the boundary, which takes time proportional to the number of boundary nodes, was computed for every boundary node for each boundary edge, so that this took time proportional to the cube of the number of boundary nodes. For 10,000 random points in a square with 400 boundary nodes, `triangulate` now takes 63 ms instead of 225 ms (the same as without a boundary), and in a square with 4,000 boundary nodes it takes 0.09 s instead of minutes. The triangulations are the same as before.
- Fixed where `triangulate` and `refine!` split the boundary near small angles of curve-bounded domains. The point at the required distance from the vertex at the small angle was found by bisection from an interval that was on the wrong side of it when searching backwards along a `BezierCurve`, `BSpline` or `CatmullRomSpline`, so that it could be far from that distance. The bisection also stopped at an absolute rather than relative tolerance, so that for any curve the point was inaccurate at small scales. Similarly, `DelaunayTriangulation.get_circle_intersection` for a `LineSegment` now returns the parameter value on the segment rather than on the part of the segment that it searched.
- The `enrich` keyword argument of `triangulate`, for enriching a piecewise linear boundary as for a curve-bounded domain, is now documented.
- Added a Package Structure page to the documentation, in the Extended Reference, describing how the repository and the source code are organised and where each feature is implemented, tested, and documented, and the main steps of `triangulate`. The mathematical details of the Bowyer-Watson algorithm now describe the biased randomised insertion order and how the point location for each new point is started.
- Fixed `triangulate` with `segments` or `boundary_nodes`, and `add_segment!`, throwing a `MethodError` for triangulations with `IntegerType = Int32` (or any integer type other than `Int`).
- Fixed `voronoi(tri; clip = true, clip_polygon, smooth = true)` and `centroidal_smooth(vorn; clip_polygon)` losing generators whose centroids were outside of the convex hull of the points, e.g. about half of them for random points in a `clip_polygon` larger than their convex hull. The convex hull was still locked while smoothing, so these generators were deleted as if they were outside of the domain. Smoothing without a `clip_polygon` is unchanged.
- Fixed clipping a Voronoi tessellation of `Float32` points to a `clip_polygon` that intersects its unbounded polygons, and clipping a power diagram of `Float32` points, throwing a `TypeError`.
- Fixed the areas and centroids of polygons that are far from the origin compared with their size, as computed by `polygon_features` and used by e.g. `get_area` and `get_centroid` for Voronoi polygons. The sums used the coordinates themselves, so their precision was lost: the polygons of a clipped tessellation of points translated by `1e8` had a total area of `-1.0`, and centroidal smoothing of points translated by `1e4` moved the generators to the wrong places, and threw an `InexactError` for points translated by `1e7`. The coordinates are now taken relative to a point of the polygon.
- `refine!` now refines a domain far from the origin as finely as the same domain at the origin. A Steiner point was rejected as being too close to a vertex if it was within `sqrt(eps)` of it relative to the size of its coordinates, so that for a domain of size `1` at `1e6` the refinement stopped at triangles of size about `0.02`, with 26% fewer vertices than at the origin, and at `1e8` almost nothing was refined. This matters e.g. for projected map coordinates, which are often around `1e6`. The tolerance is now relative to the size of the coordinates or of the domain, whichever is smaller, which makes no difference for domains near the origin.
- Fixed weighted triangulations with weights of a wider type than the points, e.g. `Float64` weights with `Float32` points such as Makie's `Point2f`, throwing a `MethodError`, and their Voronoi tessellations throwing a `TypeError`. Centroidal smoothing also no longer moves a generator to infinity if the area of its polygon is lost to rounding.
- Adding points to a `Matrix` of points, e.g. with `refine!` or `add_point!(tri, x, y)`, now throws an `ArgumentError` saying to use a `Vector` of points or an `ElasticMatrix` instead, rather than a `MethodError` from `append!`.
- Added a Performance Tips page to the manual, together with a Performance Measurements page that plots measurements of the effects it describes. The measurements are made by `docs/performance_measurements.jl` rather than when the documentation is built. The tutorial on rectangular regions no longer uses `randomise = false` for the points of a lattice, which is much slower than the default for a lattice listed row by row. The table comparing other software on the introduction page of the documentation is now the same as in the README, which it was an older version of.

## 1.7.0

- Restructured the repository to reduce package size. The package now lives in the `DelaunayTriangulation/` subdirectory and all reference images live in the top-level `assets/` directory, which is not shipped with the package. See [#237](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/issues/237).
- The minimum supported Julia version is now the LTS, Julia 1.10.

## 1.6.7

- The `Triangulation` struct now stores a `has_ghosts::Bool` field and a `boundary_vertex_to_ghost` map for efficient lookup. See [#240](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/240).
- `has_ghost_vertices(tri)` no longer scans every vertex of the graph. The scan made each query O(n), and point location hit it on every insertion, so unconstrained `triangulate` was quadratic in the number of points (591 s for 250,000 points, now 14 s). See [#236](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/issues/236) and [#247](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/247).

## 1.6.5

- Clarified the counter-clockwise requirement for the `clip_polygon` argument in `voronoi`. See [#230](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/230).
- Reduced spam during package tests. `@info` statements are now only shown on a failed test, and are now shown in the `Context` output of the failed test. See [#231](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/231).
- Removed docstrings on some methods that extended `Base` functions for internal use only. You can still find documentation on them by going to the method location and viewing the commented-out docstring. See [#231](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/231). 

## 1.6.4 

- An error is no longer thrown for inputs with duplicate points. Instead, a warning is thrown and any duplicates are merged into the `skip_points` keyword argument. With this, `DuplicatePointsError` has been removed. To silence the new warning, use `DelaunayTriangulation.toggle_warn_on_dupes!()`.

## 1.6.1
- Fix issue with clipping Voronoi tessellation dual to a single right-angled triangle. See [#207](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/207)

## 1.6.0
- Define `reverse` for `AbstractParametricCurve`s, making it easier to reverse the orientation of a curve. See [#195](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/195).
- Fixed an issue with `LineSegment` not returning the exact endpoints at `t=1`, which can be problematic when joining boundary nodes. This has been fixed. See [#195](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/195).
- Introduced `is_linear` to fix issues with boundary enrichment of domains with `LineSegment`s. In particular, `LineSegment`s are no longer enriched. See [#195](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/195).
- `orientation_markers` now uses `uniquetol` instead of `unique` for the final set of markers (it already did it for the intermediate markers). See [#195](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/195).
- For large `Tuple`s, functions like `eval_fnc_at_het_tuple_two_elements` are problematic and allocate more than their non-type-stable counterparts. To get around this, for `Tuple`s of length `N > 32`, the non-type-stable version is used. See [#195](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/195).
- Fixed issue with `use_barriers` when a ghost edge is selected at random during point location. See [#196](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/196).
- Introduced the (currently internal) function `get_positive_curve_indices` for finding curves with positive orientation in a `Triangulation`. See [#196](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/196).
- `is_exterior_curve`, `is_interior_curve`, `num_exterior_curves`, and `is_disjoint` are now defined based on `get_positive_curve_indices` rather than `get_exterior_curve_indices`. See [#196](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/196).
- PrecompileTools.jl is now used. See [#200](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/200).
- Introduced the (currently internal) function `get_positive_curve_indices` for finding curves with positive orientation in a `Triangulation`. [#196](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/196).
- `PointLocationHistory` was not marked as public. This has been fixed. See [#198](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/198).
- Fixed an issue with missing docstrings and duplicate docstrings in the documentation. See [#198](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/198).
- `copy` and `deepcopy` are now correctly implemented for `PolygonTree`s and `PolygonHierarchy`s. See [#199](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/199)
- Implemented `copy` and `deepcopy` for `Triangulation` and `VoronoiTessellation`.  See [#201](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/201). 
- Fixed a bug with `Triangulation`s `polygon_hierarchy` not being correctly aliased with the `polygon_hierarchy` from the `BoundaryEnricher`, and similarly for the `boundary_edge_map`.  See [#201](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/201).
- Implemented `==` for `VoronoiTessellation`. See [#201](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/pull/201).

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
