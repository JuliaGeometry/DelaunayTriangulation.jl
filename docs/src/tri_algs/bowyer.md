```@meta
CurrentModule = DelaunayTriangulation
```

# Bowyer-Watson Algorithm

Here we will give a description of the Bowyer-Watson algorithm. This algorithm is the algorithm used by `triangulate`.

The main idea behind the algorithm is to insert points one at a time, deleting triangles at each step whose circumcircles contain the point to be inserted, then repairing the cavity. We give the procedure in steps below.

1. First, using `initialise_bowyer_watson`, we need to define the initial triangle. This triangle will be the first three points in the provided `point_order` (which is, by default, just a random permutation of the point indices - see the `point_order` keyword argument in `triangulate_bowyer_watson` and the `get_point_order` function). If these first three points are collinear, we do a circular shift of the point order until we get a non-degenerate triangle. See `get_initial_triangle`. At this step, we also reset the representative points field of the triangulation and initialise it with `DelaunayTriangulation.BoundaryIndex` mapping to the centroid of this initial triangle.
2. Once the initial triangle is selected, we move into `_triangulate_bowyer_watson`, where we loop over each point and add it in one at a time. In this loop, we start by selecting the initial point to start the jump and march algorithm at, making use of `select_initial_point`. Since points are added in one at a time, and a user's insertion order may often have points that are close together both in the order and in space, the `try_last_inserted_point` keyword argument is useful here in case we can start right next to the new point. With this point selected, we move into actually adding the point via `add_point_bowyer_watson!`.
3. The `add_point_bowyer_watson!` starts by using `jump_and_march` to find a triangle containing the point. The idea is to then find all triangles whose circumcircles, i.e. the circle through the three points of the triangle, contain this new point. These points need to be deleted since, by definition, these triangles are no longer Delaunay. This is done via a depth-first search, where we take the triangle we are currently in and step over its three edges into three new triangles, done via the recursive function `dig_cavity!`. If the new triangle also contains the point in its circumcircle, we delete it also, and we keep stepping. We stop at any triangles that don't contain the point in its circumcircle. Once we have stopped, we take the edge we did not step over and connect it with the new point, giving us a new triangle. 
4. Still in `add_point_bowyer_watson!`, an important case to consider is when the point we find is directly on the triangle we found. This does not cause any problems with `dig_cavity!`, but it may cause issues with how we update the boundary, so not only do we check if the point is on the triangle, but we also check that the triangle is either a boundary triangle or a ghost triangle (meaning the edge is on the boundary). If this is the case, then we find the edge of the triangle that the point is on with `find_edge`, and split the edge in half at the point, placing the new point correctly on the boundary and giving two new triangles.
5. Lastly, still in `add_point_bowyer_watson!`, we use `update_centroid_after_addition!` to update the centroid of the points with the new point.
6. Steps 2--5 are repeated for each new point, until we have finally added all points. Once this is done, we compute the convex hull of the points with `convex_hull!`, stepping over the boundary using the ghost triangles from the triangulation to get all the boundary nodes efficiently.
7. Next, if the keyword argument `recompute_representative_point` is true, we can give a better representative point for the central part of the domain than the centroid by computing the pole of inaccessibility. This is done with `compute_representative_points!`.
8. Finally, to clean up, we can delete all ghost triangles (if the keyword argument `delete_ghosts` is true) with `delete_ghost_triangles!`. Then, if the keyword argument `delete_empty_features` is true, we can delete all keys from the `Adjacent` map that map to empty values with `clear_empty_features!`, which would also clean up empty sets from the `Adjacent2Vertex` map and empty neighbourhoods from the `Graph`.

## Modifications for a constrained Delaunay triangulation

The Bowyer-Watson algorithm requires two modifications for adding points into a constrained Delaunay triangulation. The modifications are:

1. Avoid walking over any constrained edges when performing the depth-first search.
2. If a point is added onto a constrained segment, split the segment in two and perform the depth-first search on each side of the segment.

With just these two modifications, the algorithm works.