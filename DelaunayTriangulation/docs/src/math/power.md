# Power Diagrams

Now we discuss power diagrams, the dual of a weighted Delaunay triangulation. The main difference between a power diagram and a Voronoi tessellation is that the power diagram is defined using a _power distance_ instead of a _Euclidean distance_. Recall that weighted Delaunay triangulations are derived from a point set $\mathcal P$ and an associated set of weights $\mathcal W$ so that each point $p_i$ has an associated weight $w_i$. If we let $p[w]$ denote a point $p$ with weight $w$, then the power distance between two points $p$ and $q$ is defined by
```math 
\pi(p[w_p], q[w_p]) = d(p, q)^2 - w_p - w_q,
```
where $d(p, q)$ is the Euclidean distance between $p$ and $q$. This notation is also used for unweighted points, which are assigned a weight of $0$.

We can use this power distance to define a power diagram. For a given site $u$, its _power cell_ is 
```math 
W_u = \{p \in \mathbb R^2 : \forall q \in \mathcal P^+, \pi(p, u[w_u]) \leq \pi(p, q[w_q])\},
```
where $\mathcal P^+$ is the set of points $\mathcal P$ with the assigned weights from $\mathcal W$. In particular, $W_u$ is the set of all points closer to $u$ than to any other point in $\mathcal P$, where the distance is measured using the power distance.

## Orthoballs and Orthocenters: Generalising Circumcircles

Recall that Voronoi tessellations are computed by simply connecting the circumcenters of the Delaunay triangles. For power diagrams, we need a generalisation of the circumcenter. Remember that the circumcenter can be derived by intersecting the perpendicular bisectors of a triangle's edges. A similar definition can be used for obtaining what is called the _orthocenter_.

To start, we need a geometric interpretion of a weighted point. We can relate it to the notion of a [power of a point](https://en.wikipedia.org/wiki/Power_of_a_point). Define the ball $B_p = \{x \in \mathbb R^2: d(p, x) \leq \sqrt{w_p}\}$, and note that the squared length of the line segment extending from a point $x$ to touch $B_p$ tangentially is given by $d(x, p)^2 - w_p$ (see the first equation on the linked Wikipedia page).  This distance is exactly $\pi(x, p[w_p])$. Thus:
- If $x \notin B_p$, then $\pi(x, p[w_p]) > 0$.
- If $x$ is exactly on the boundary of $B_p$, then $\pi(x, p[w_p]) = 0$.
- If $x \in B_p$, then $\pi(x, p[w_p]) < 0$.
Note that, since weights may be negative, $\sqrt{w_p}$ may be complex, in which case every point lies outside $B_p$. Thus, our geometric interpretation of a weighted point is exactly as the ball $B_p$.

We can use this interpretation to understand how the generators of the power diagram are related to the vertices of the power cells. Remember that the definition of a power cell implies that a vertex $v$ of the diagram if equidistant in power distance from the weighted generators that generate it, e.g. if $v$ neighbours generators $u$ and $w$, then $\pi(v, u[w_u]) = \pi(v, w[w_w])$. Let this constant power distance be $\omega_v$. Then, the power distance from $v$ to its weighted generators is 
```math 
\pi(v[\omega_v], u[w_u]) = d(v, u)^2 - \omega_v - w_u,
```
but we also know $\pi(v, u[w_u]) = \omega_v = d(v, u)^2 - w_u$. Thus,
```math 
\pi(v[\omega_v], u[w_u]) = 0.
```
We say that $v[\omega_v]$ is _orthogonal_ to $u[w_u]$. With this definition, $v[\omega_v]$ is _orthogonal_ to all of its generators. Letting $B_v$ be the ball centered at $v$ with radius $\sqrt{\omega_v}$ (which could be complex), we call $B_v$ the _orthoball_ centered at $v$. More generally, we say that two weighted points are
1. _Orthogonal_ if the power distance between them is zero;
1. _Farther than orthogonal_ if the power distance between them is positive;
1. _Closer than orthogonal_ if the power distance between them is negative.

The orthoball can be used to generalise circumcircles. The _orthoball_ of a triangle is defined to be the ball that is orthogonal to all three vertices of the triangle. In particular, if $c$ if the ball's center, called the _orthocenter_, and $\sqrt{\omega_c}$ is its radius, called the _orthoradius_, then $d(v, c)^2 - w_v - \omega_c = 0$ for each vertex $v$ of the triangle. The center and radius can be computed by solving this system of equations $\{d(v, c)^2 - w_v - \omega_c = 0\}$. It is possible to show these orthoballs are exactly what we need to define duality between the power diagram and the weighted Delaunay triangulation. 

## Computing the Power Diagram

The above discussion leads us to a method for computing the power diagram. Just as a Voronoi tessellation can be computed by joining up the circumcenters, simply join up the orthocenters of each triangle. Special care may be needed to accommodate the fact that not all generator's power cells are non-empty.