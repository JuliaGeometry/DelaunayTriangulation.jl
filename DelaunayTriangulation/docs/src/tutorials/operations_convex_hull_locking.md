```@meta
EditURL = "https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/operations_convex_hull_locking.jl"
```

# Triangulation Operations
## Locking and Unlocking the Convex Hull

For unconstrained triangulations, the only boundary
is the convex hull of the point set, but the
`boundary_nodes` field is empty because it is reserved
for a constrained boundary. There may be cases where
you want to treat the convex hull as if it were a
constrained boundary. For example, this is done internally
inside `refine!` when providing an unconstrained triangulation
for mesh refinement. Let us give an example of how this can be done,
in case you want to do this for your own application.

````@example operations_convex_hull_locking
using DelaunayTriangulation
using CairoMakie

points = rand(2, 50)
tri = triangulate(points)
get_boundary_nodes(tri)
````

As you can see, the boundary nodes field is empty.
We can lock the convex hull using [`lock_convex_hull!`](@ref):

````@example operations_convex_hull_locking
lock_convex_hull!(tri)
get_boundary_nodes(tri)
````

Now the boundary nodes field is not empty. Note that if you try
and lock the convex hull again, you will get an error because
`DelaunayTriangulation.has_boundary_nodes(tri)` is now true.
To now unlock the convex hull, we use [`unlock_convex_hull!`](@ref):

````@example operations_convex_hull_locking
unlock_convex_hull!(tri)
get_boundary_nodes(tri)
````

This function will error if it detects that the existing boundary
isn't actually equal to the convex hull.

Note that this locking/unlocking doesn't actually change anything about the triangulation,
it just adds information into `tri` to treat it as if you had provided
the convex hull as a constrained boundary to start with.

## Just the code
An uncommented version of this example is given below.
You can view the source code for this file [here](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/operations_convex_hull_locking.jl).

```julia
using DelaunayTriangulation
using CairoMakie

points = rand(2, 50)
tri = triangulate(points)
get_boundary_nodes(tri)

lock_convex_hull!(tri)
get_boundary_nodes(tri)

unlock_convex_hull!(tri)
get_boundary_nodes(tri)
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

