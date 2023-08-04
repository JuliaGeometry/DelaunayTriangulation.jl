# # Example Applications 
# ## Edge Flipping 

# In this application, we implement the edge flipping algorithm for computing an 
# (unconstrained) Delaunay triangulation. This algorithm, 
# as described by de Berg et al. in the book *Computational Geometry: Algorithms and Applications* (3 ed., 2008),
# is as follows: For a given point set $P$ and some initial triangulation $\mathcal T$ 
# of the point set:
#
# * while $\mathcal T$ contains an illegal edge $e_{ij}$
#   * do
#       * find the triangle $T_{ijk}$ that contains $e_{ij}$* Let $T_{ijk}$ and $T_{ji\ell}$ be the two triangles adjacent to $e_{ij}$.
#       * Remove $e_{ij}$ from $\mathcal T$, and add $e_{k\ell}$ instead. That is, flip the edge $e_{ij}$.
# * return $\mathcal T$
#
# Here, $e_{ij}$ denotes the edge $\overline{p_ip_j}$, where $p_i$ is the $i$th point, and $T_{ijk}$ is the triangle $\triangle p_ip_jp_k$. 
# The edge $e_{ij}$ is *illegal* if the circle through
# $p_i$, $p_j$, and $p_k$ contains $p_\ell$ in its interior, where 
# $p_k$ and $p_\ell$ are the other points on the triangles adjoining $e_{ij}$.
# The idea behind the algorithm is that flipping illegal edges _increases_ the minimum angle 
# of the triangulation, and thus, since the Delaunay triangulation is the triangulation with the
# largest minimum angle, the algorithm converges to the Delaunay triangulation and does 
# so in a finite number of steps. 
# 
# One complication with the algorithm is that we need to find a way to 
# construct an initial triangulation $\mathcal T$ for the given set of points. 
# One way to do this is to construct the convex hull, pick one point on the convex hull 
# and connect all other points on the convex hull to that point, and then for all remaining points 
# in the set, find the triangle that contains the point and add the point to the triangulation, connecting 
# that point to the triangle's vertices. Here is one implementation of this. 
using DelaunayTriangulation
using CairoMakie 

function get_initial_triangulation(points)
    tri 
end