```@meta
CurrentModule = DelaunayTriangulation
```

# Ruppert's Algorithm for Mesh Refinement 

One problem with triangulations is that just providing a boundary is not enough to give a good mesh - we would ideally want a good distribution of triangles. The way that this is done is through Delaunay refinement. Following Chapter 6 of the Delaunay Mesh Generation book by Cheng, Dey, and Shewchuk, a general procedure for Delaunay refinement is:

1. Compute $\mathcal D\mathcal T(\mathcal S)$, where $\mathcal S = \{\mathcal P, \mathcal E, \mathcal B\}$ is a given set of points, edges, and boundary edges. 
2. If $\mathcal D\mathcal T(\mathcal S)$ does not meet some desired geometric or topological property, insert some point near the violation into $\mathcal D\mathcal T(\mathcal S)$ to meet the property. Update $\mathcal D\mathcal T(\mathcal S)$.
3. If some element $\mathcal T \in \mathcal D \mathcal T(\mathcal S)$ is too large or too small, put a point the element and update $\mathcal D\mathcal T(\mathcal S)$ accordingly.
4. Do the above until all elements are of good quality.

