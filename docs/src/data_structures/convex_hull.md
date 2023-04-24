```@meta
CurrentModule = DelaunayTriangulation
```

# Convex Hull 

We represent convex hulls using a `ConvexHull` type, which is simply a type containing points and indices:

```@docs 
ConvexHull
```

For a triangulation, convex hulls are obtained from the unconstrained form, but if they need to be reconstructed then we can do so with the monotone chain algorithm.

```@docs 
convex_hull(::Any)
```

