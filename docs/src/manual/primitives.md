# Representing Primitives 
In this package, we need methods for representing points, edges, triangles, and boundaries. By default, we provide a simple interface for defining these primitives, typically treating these primitives as:

- Points: `Vector{NTuple{2, Float64}}` or `Matrix{Float64}`, where in the latter case points are defined by the columns of the matrix. Where mutability is needed, typically the former is preferred as an `ElasticMatrix` from `ElasticMatrix.jl` is needed if a matrix is used instead.
- Edges: `Set{NTuple{2, Int}}`.
- Triangles: `Set{NTuple{3, Int}}`.
- Boundaries: The recommended way to represent a boundary depends on the type, as described in the [Representing Boundaries](boundaries.md) section. For a contiguous boundary, `Vector{Int}` is preferred; for a boundary with a single sectioned curve, `Vector{Vector{Int}}` is preferred; for a boundary with multiple disjoint sections, `Vector{Vector{Vector{Int}}}` is preferred.

For defining new methods, as discussed in [this tutorial](../tutorials/custom_primitive.md) and this [API reference](../api/primitives.md), there are many methods that could be overloaded. When defining these methods, though, there are some important things to consider:

- Points are constantly being accessed via index, so a vector-like structure is recommended, and fixed-size tuples (like `Tuple` or a `SVector` from StaticArrays.jl) are recommended for individual points. For mutability, e.g. in the case of refinement or translating points with centroidal Voronoi tessellations, you of course want a mutable structure.
- Edges and triangles are never accessed via index, so a vector-like structure is not needed. Instead, these collections are always being used for finding specific objects and for deleting/pushing objects, so a set-like structure is highly recommended. For individual edges and triangles, fixed-size tuples of length 2 and 3, respectively, are recommended.
- `Int` is recommended to be used for vertices, like in an edge or a triangle, and for coordinates `Float64` is highly recommended. As discussed in the [Geometrical Predicates](predicates.md) section, the use of `Float64` is important for the robustness of the algorithms in this package - using `Float32` is highly likely to lead to errors or infinite loops in the algorithms; as of v1.1.0, the use of `AdaptiveKernel()` predicate kernel makes issues from `Float32` less likely.

Of course, how you implement the methods for your custom types comes down to what you need. 

For using these methods inside [`triangulate`](@ref), you need to specify all types other than the point used for points, as described in its docstring.