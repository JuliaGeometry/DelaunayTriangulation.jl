```@meta
CurrentModule = DelaunayTriangulation
```

# Point Location

The most expensive step for building a Delaunay triangulation is the point location step, wherein we need to find a triangle that contains a given point. The code for this turns out to be very complicated so that we can correctly handle points outside of the domain, inside interior holes, collinear with other points, on the corner, etc. The main function that handles all of this is `jump_and_march`, derived from the jump-and-march algorithm of Mücke, Saias, and Zhu (1999). We also provide a method that simply searches over all triangles, given below.

```@docs
brute_force_search
```

## The Main Method 

Below we list docstrings for the main jump and march algorithm.

```@docs 
jump_and_march 
```

The variant of the jump and march algorithm for points outside of the triangle is also accessed via the `jump_and_march`, calling into `exterior_jump_and_march`:

```@docs 
exterior_jump_and_march
```

You should not need to call into this method directly.

The core complexity of the algorithm comes from having to find the direction of the point from an initial search point. The docstrings for some of these initialisers are given below.

```@docs 
default_num_samples
select_initial_point 
select_initial_triangle_interior_node 
check_for_intersections_with_adjacent_boundary_edges
search_down_adjacent_boundary_edges
check_for_intersections_with_interior_edges_adjacent_to_boundary_node
```

## History 

If you need to, you can also store the history of the algorithm. This is primarily only used for (and only tested for) constrained triangulations, as it helps us locate which triangles to delete. The struct we use for storing history is given below.

```@docs 
PointLocationHistory 
```

## Basic Description of the Algorithm

Let us give a basic description of what `jump_and_march` does.

First, using `select_initial_point`, an initial point to start the algorithm is selected. This function samples some number $m$ of points, and then selects the point that is closest to the query point $q$.

Let the initially selected point be $p_k$. We break the discussion into two cases, where $p_k$ is an interior point and $p_k$ is a point on the boundary.

If $p_k$ is not a point on the boundary, then it is possible to completely rotate around the point searching for a triangle such that line $\overrightarrow{p_kq}$ intersects an edge of the triangle. This will give us an edge $e_{ij}$ that $\overrightarrow{p_kq}$ intersects, and we will put $p_i$ to the left of $\overrightarrow{p_kq}$ and $p_j$ to the right. The function that handles this selection is `select_initial_triangle_interior_node`, which starts by handling the case of collinear edges and then rotates around. 

Now suppose that $p_k$ is a point on the outer boundary. There are several possibilities in this case. First, $\overrightarrow{p_kq}$ might intersect a neighbouring boundary edge; secondly, $\overrightarrow{p_kq}$ could intersect a neighbouring interior edge; thirdly, $\overrightarrow{p_kq}$ could point away from the boundary, meaning $q$ is outside of the boundary. The function starts by checking the neighbouring boundary edges, done via `check_for_intersections_with_adjacent_boundary_edges`, making use of `get_right_boundary_node` and `get_left_boundary_node` to obtain the neighbouring boundary nodes. If we find that $\overrightarrow{p_kq}$ does intersect a neighbouring boundary edge, then we can search down adjacent boundary edges via `search_down_adjacent_boundary_edges` until we either find an edge that $q$ is on, or until we identify that $q$ is outside of the triangulation -- this function assumes the domain is convex, as triangulations being built are. If $q$ is outside of the triangulation, then we can use the exterior variant of the jump and march algorithm, `exterior_jump_and_march`, to find the ghost triangle containing $q$. This functon simply rotates around the boundary until we find two ghost edges enclosing the point $q$. Now let us assume we did not find a neighbouring boundary edge that intersects $\overrightarrow{p_kq}$. When this happens, we need to check the neighbouring interior edges, done via `check_for_intersections_with_interior_edges_adjacent_to_boundary_node`, which just rotates around the edges until we find an intersection. If we have still not found any intersection then, again assuming convexity, $q$ must be outside of the triangulation and so we use `exterior_jump_and_march`.

Now, if the algorithm is still going, then we need to start marching along the triangulation from $p_k$ towards $q$. The idea is to keep marching, keeping $p_i$ and $p_j$ to the left and right of $\overrightarrow{p_kq}$, respectively, until we find a case where $p_ip_jq$ is no longer a positively oriented triangle. When this happens, it must mean that we have passed $q$, and so we have found the triangle. 

Let us describe how this marching is done in more detail. First, we need to be careful of boundary indices, first checking for an outer boundary index. If we have found an outer boundary index, then we have marched into the boundary, and so $q$ will be outside of the domain, meaning we go into `exterior_jump_and_march`. If this is not the case, then we march using `get_adjacent` to step onto the next triangle from a given edge $e_{ij}$. Then, checking the positions of the new points relative to $\overrightarrow{p_kq}$ and rearranging accordingly, we can step forward. We keep doing this until we get a negatively oriented triangle $p_ip_jq$. Unfortunately, this loop could terminate even if $q$ is not in the found triangle, which can occasionally happen if $p_ip_jq$ is a degenerate triangle. In this case, we just restart the jump and march algorithm. This latter worry is a very rare concern and does not alter the runtime in any significant manner.