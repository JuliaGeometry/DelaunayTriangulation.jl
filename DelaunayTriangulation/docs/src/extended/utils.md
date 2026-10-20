```@meta 
CurrentModule = DelaunayTriangulation
```

# Utility Functions

This section lists some of the internal utility functions, or other miscellaneous functions, used in this package.

```@docs; canonical=false
number_type
```

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["src/utils/utils.jl"]
Filter = t -> !(t in (DelaunayTriangulation.number_type,))
```

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["src/geometry/geometry_utils.jl"]
```

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["src/predicates/predicate_definitions.jl"]
```

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["src/predicates/boundaries_and_ghosts.jl"]
```

```@docs
convert_certificate
DefaultAdjacentValue
𝒢
GhostVertex
ε
∅
fix_orient3_cache 
fix_incircle_cache 
validate_orient3_cache
validate_incircle_cache
```