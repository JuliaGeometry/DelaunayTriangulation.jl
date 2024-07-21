"""
    DefaultAdjacentValue = 0 

Default value used for representing an empty result 
from an adjacency query.
"""
const DefaultAdjacentValue = 0

"""
    ∅ = DefaultAdjacentValue

Alias for [`DefaultAdjacentValue`](@ref).
"""
const ∅ = DefaultAdjacentValue

"""
    GhostVertex = -1

Number used for representing initial ghost vertices. 
All other ghost vertices are derived from subtracting from 
this number. See https://juliageometry.github.io/DelaunayTriangulation.jl/stable/manual/ghost_triangles/.
"""
const GhostVertex = -1

"""
    𝒢 = GhostVertex

Alias for [`𝒢`](@ref).
"""
const 𝒢 = GhostVertex

"""
    ε = sqrt(eps(Float64))

Number used as a tolerance in certain functions, e.g. 
for mesh refinement when using [`check_precision`](@ref) to 
avoid degenerate circumcenters.
"""
const ε = sqrt(eps(Float64))

const INF_WARN = Ref(true)
"""
    toggle_inf_warn!()

Toggle the warning for infinite circumcenters in the Voronoi tessellation.
By default, this warning is enabled.
"""
toggle_inf_warn!() = (INF_WARN[] = !INF_WARN[])

@static if VERSION ≥ v"1.6"
    using Preferences
end

@static if VERSION ≥ v"1.6"
    const PREDICATES = @load_preference("PREDICATES", "EXACT")::String # This default is not guaranteed to be consistent between versions
else
    const PREDICATES = "EXACT"
end
@static if PREDICATES ∉ ("EXACT", "INEXACT")
    throw("You have set the PREDICATES option to PREDICATES = $PREDICATES. This is not allowed, only EXACT or INEXACT are possible choices.")
end

@doc """
    PREDICATES 

Type of predicates to use. This can be either 

- `PREDICATES = "EXACT"`: Use ExactPredicates.jl.
- `PREDICATES = "INEXACT"`: Compute the predicates numerically without any extra checks.

(In the future, this may include `"ADAPTIVE"`.) By default this is assumed to be 
`"EXACT"`, but this default is not guaranteed to be the same across versions. 
You can change this setting using Preferences.jl. For example, to 
use inexact predicates, do
```julia-repl
julia> using Preferences: set_preferences!

julia> set_preferences!("DelaunayTriangulation", "PREDICATES" => "INEXACT")

julia> using DelaunayTriangulation # load only after setting the preference

julia> DelaunayTriangulation.PREDICATES
"INEXACT"
```

!!! note "Precision"

    Regardless of the setting, all coordinates are computed to `Float64` before 
    computing any predicates.
"""
PREDICATES

const USE_EXACTPREDICATES = PREDICATES == "EXACT"
const USE_INEXACTPREDICATES = PREDICATES == "INEXACT"
