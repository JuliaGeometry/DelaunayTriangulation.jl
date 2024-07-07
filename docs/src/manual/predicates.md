```@meta 
CurrentModule = DelaunayTriangulation
```

# Geometrical Predicates

This section discusses how geometrical predicates are defined in this package. The predicates in this package are primarily derived from those implemented in [ExactPredicates.jl](https://github.com/lairez/ExactPredicates.jl). The choice of exact predicates is important for the robustness of the algorithms in this package. Without using exact predicates, you may quickly find issues such as infinite loops or errors in the algorithms, as discussed for example at the start of p.3 of [these notes](https://perso.uclouvain.be/jean-francois.remacle/LMECA2170/robnotes.pdf) by Shewchuk. If you do want to disable exact predicates, see [here](disabling_ea.md).

## Certificates 

All predicates defined in this package return a [`Certificate`](@ref) which simply specifies the result of the predicate. This is easier than working with `Bool`s only as (1) not all predicates have only two outcomes and (2) it is easier to see the certificate than to remember exactly what outcome is represented by a `Bool`. For example:

```@example certex 
using DelaunayTriangulation
p, q, r = (0.0, 0.0), (1.0, 0.0), (0.0, 1.0)
flag = DelaunayTriangulation.triangle_orientation(p, q, r)
```

We could then inspect the result using e.g.

```@example certex 
DelaunayTriangulation.is_positively_oriented(flag)
```

## Predicates

The predicates we define, as mentioned, are primarily derived from those in [ExactPredicates.jl](https://github.com/lairez/ExactPredicates.jl), where we simply extend the [`orient`](@ref orient_predicate), [`incircle`](@ref incircle_predicate), [`parallelorder`](@ref parallelorder_predicate), [`sameside`](@ref sameside_predicate), and [`meet`](@ref meet_predicate) predicates, allowing us to define predicates for, for example, [`triangle_orientation`](@ref) and [`point_position_relative_to_line`](@ref). Predicates for working with the boundary and ghost vertices are also implemented, for example, [`is_ghost_edge`](@ref is_ghost_edge) and [`is_boundary_node`](@ref).