```@meta
CurrentModule = DelaunayTriangulation
```

# Interface 

The package makes its simple for customing the interface used for defining points, edges, and triangles, as described in the docstring for `Interfaces` (see the end of this section). Without any customisation, the default forms are:

- Edges: `NTuple{2, Int}`.
- Collections of edges: `Set{NTuple{2, Int}}`.
- Triangles: `NTuple{3, Int}`.
- Collections of triangles: `Set{NTuple{3, Int}}`

We also give support for customing how points are represented, and by default we support collections of points given as matrices (with each point its own column), or vectors of vectors. The number type used for representing coordinates has to be Float64 to support ExactPredicates.jl, although if you like you could customise `orient_predicate` and `incircle_predicate`, even circumventing ExactPredicates.jl if you like. See the predicates section in the sidebar for a further discussion of changing these predicate definitions.

We also provide a customisable interface for representing boundary nodes, although for unconstrained triangulations this is relevant. By default, we support boundary nodes represented according to the following, where we let `BN` refer to the collection of boundary nodes:

- `Vector{Int}`: In this case, there is only one fixed boundary and it is represented as a contiguous set of nodes. We must have `BN[begin] == BN[end]`, and the nodes must be listed in counter-clockwise order.
- `Vector{Vector{Int}}`: In this case, there is only one fixed boundary, but it is made up of separate segments, with `BN[n]` the nodes for the `n`th segment. This makes it possible to more easily support, for example, a domain with different boundary conditions on separate parts of the boundary. We must have `BN[n][end] == BN[n+1][begin]` and `BN[end][end] == BN[begin][begin]`, and each segment must be listed in counter-clockwise order.
- `Vector{Vector{Vector{Int}}}`: In this case, there are multiple fixed boundaries, each of which are assumed to be made up of separate segments. This makes it possible to support multiply-connected domains, e.g. an annulus with each circle split into its lower and upper halves. Here, `BN[m][n]` is the set of nodes for the `n`th segment of the `m`th boundary curve, and `BN[begin]` the outer-most boundary curve and `BN[m]`, `m > 1`, nodes for curves contained within `BN[begin]`. As in the previous case, `BN[m][n][end] == BN[m][n+1][end]` and `BN[m][end][end] == BN[m][begin][begin]` for each `m`. Moreover, `BN[begin]` should be a counter-clockwise list of nodes while `BN[m]` is a clockwise list of nodes for `m > 1`.

For more information about how we handle boundaries, and how they are handled in our data structures, see the boundary handling section in the sidebar. 