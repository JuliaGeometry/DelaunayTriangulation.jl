```@meta
CurrentModule = DelaunayTriangulation
```

# Chew's Algorithm for Triangulating Convex Polygons

Our algorithm for triangulating convex polygons comes from L. Paul Chew's "Building Voronoi diagrams for convex polygons in linear expected time" back in 1990, following the presentation of it in the book *Delaunay Mesh Generation* by Cheng, Dey, and Shewchuk (2013). The idea is that we can utilise the neighbour information already provided by the definition of a polygon to avoid the point location step, which turns out to be the most expensive step when computing triangulations normally.

A basic overview of the algorithm is as follows: Let us suppose that $S$ is a sequence that lists the $k$ vertices of some convex polygon in counterclockwise order, say $S = \{v_1, \ldots, v_k\}$. To improve the complexity of our algorithm, we should randomise the insertion order for the polygon vertices to reduce the work needed, noting that the average degree of a vertex in a triangulation should typically be around six. So, we construct some permutation $\pi$ of $\{1, \ldots, k\}$ that will define how we insert the points, taking the first three vertices $\pi_1, \pi_2, \pi_3$ as our initial triangle. Then, using the Bowyer-Watson algorithm, we can add in points one at a time, simply walking along the polygon to find triangles rather than the jump and march algorithm used in our main Bowyer-Watson algorithm. Note that when we add some vertex $u$, only one edge can separate it from the triangulation's interior, else the polygon wouldn't be convex. Hence, when $u$ is inserted, we just take the union of this disk with all other triangles whose open circumdisks contain $u$. This will give us our triangulation.

Let us now describe the algorithm in more detail. The main complexity in the algorithm comes in from devising a smart way for avoiding point location. This is done with backward analysis, where instead of imagining points in one at a time, we think about what happens if we delete vertices one at a time. Suppose we have a vertex $u$ with neighbours $v$ and $w$. When we delete $u$, the only way to reform the poylgon and keep it convex is to join $v$ and $w$, i.e. a deletion of a vertex also implies a new edge from its neighbours. This edge is exactly what we need for point location: the edge we need is simply the edge from a point's neighbours, thus allowing us to skip the jump and march algorithm. With this knowledge, we construct a circularly-, doubly-linked list of our polygon vertices, called $\mathcal N$ (next) and $\mathcal P$ (previous) for the two neighbours of each vertex, and then we imagine deleting vertices one at a time by simply editing this linked list. The steps are as follows:

1. First, take $S = \{v_1, \ldots, v_k\}$ and a permutation $\pi$ of $\{1,\ldots,k\}$.
2. Next, compute $\mathcal N = \{2, 3, \ldots, 1\}$ and $\mathcal P = \{k, 1, \ldots, k-1\}$.
3. Then, from $i=k$ down to $i=4$:
    a. To delete $v_{\pi_i}$ from the polygon, we can set $\mathcal N[\mathcal P(\pi_i)] \mapsfrom \mathcal N(\pi_i)$ and $\mathcal P[\mathcal N(\pi_i)] = \mathcal P(\pi_i)$, i.e. just connect the two neighbours together.
4. Now, start the triangulation with the initial triangle through $(v_{\pi_1}$, $v_{\mathcal N(\pi_1)}$, $v_{\mathcal P(\pi_1)}$. If this triangle is degenerate, just repeat the above steps until you find a non-degenerate triangular.
5. Now that we have our initial triangle, we can start adding points in. So, for $i = 4$ up to $i = k$:
    a. Call the routine $\mathcal C(v_{\pi_i}, v_{\mathcal N(\pi_i)}, v_{\mathcal P(\pi_i)})$ to insert $v_{\pi_i}$, defined below.

This routine $\mathcal C(u, v, w)$ is defined as follows:

1. First, take $x = \mathcal A(w, v)$ where $\mathcal A$ is the adjacent map. With this definition, the triangle $wvx$ is opposite the ede $vw$ from $u$.
2. It is possible that $x$ is not on the polygon, i.e. $x$ might be $\emptyset$ or, if this algorithm is being used for deleting a vertex so that the polygon is embedded inside a larger triangulation, $x$ could be a vertex away from the polygon, so we must check that $x \in S$. Moreover, we must check that $x$ is actually inside the triangle through $u$, $v$, and $w$. If these two conditions hold, then the triangles $uvw$ and $wvx$ must not be Delaunay, so delete the triangle $wvx$, and continue stepping through the triangulation by calling $\mathcal C(u, v, x)$ and $\mathcal C(u, x, w)$ so that more non-Delaunay triangles can be identified and deleted.
3. If the two conditions in the last step did not hold, we start this third step (otherwise, return in step 2 above). These two conditions holding mean that $vw$ is still a Delaunay edge, so we can edge the triangle $uvw$ into the triangle.

With this routine and the steps for the original algorithm, we can triangulate convex polygons.