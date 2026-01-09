```@meta 
CurrentModule = DelaunayTriangulation
```

# Representing Boundaries

Here we will give a description of how we represent boundaries when considering constrained triangulations. There are three possible types of boundaries that can be considered:

1. A contiguous boundary with no holes. For example, a circle.
2. A single boundary with no holes, but the boundary is split into multiple sections that can be identified separately via ghost vertices (ghost vertices are described in the [next section](ghost_triangles.md)). For example, a square with a different ghost vertex for each side.
3. A boundary comprising multiple disjoint curves and possibly disjoint domains. For example, a square with a circular hole inside and a separate circular domain.

The way to represent boundaries can be customised as needed, but any such interface for the boundaries must conform to the following specifications; these specifications are what we use in the `convert_boundary_points_to_indices` function.


## Contiguous boundary

In the case of a contiguous boundary, the boundary should be defined in a counter-clockwise orientation (meaning that, if you were to walk along the boundary, the interior of the domain would be on your left). Moreover, the first and last points of the boundary must be the same. For example, if 

```@example contgex
points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
```

Then the boundary nodes should be defined as 

```@example contgex 
boundary_nodes = [1, 2, 3, 4, 1]
```

The vector `[1, 2, 3, 4]` would not be correct as the first and last points are 
not correct, and `[1, 4, 3, 2, 1]` would not be correct as the orientation is
clockwise.

## Sectioned boundary 

In the case of a single boundary that is split into sections, the boundary should be defined in a counter-clockwise orientation, and the first and last points of each section must be the same. For example, suppose we have 

```@example sectex 
points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
```

and we want to define each side of the square to be its own section of the boundary. We should then define the boundary nodes as 

```@example sectex 
boundary_nodes = [[1, 2], [2, 3], [3, 4], [4, 1]]
```

The vector `[[1, 2], [2, 3], [3, 4]]` would not be correct as the first and last points of the boundary are not the same, and `[[1, 4], [4, 3], [3, 2], [2, 1]]` would not be correct as the orientation is clockwise.

## Multiply-connected domain

The specification of a multiply-connected domain is a bit more cumbersome. For each individual boundary, the specification is the same as that of a sectioned boundary, but more care is needed for the orientations. The orientation of the entire domain's boundary must be positive, but individual boundaries may need to be defined in a clockwise orientation so that the domain's interior remains on the left when it is traversed. For example, if a circle is to be defined inside of a square boundary, then the circle would have to be defined clockwise.