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
    const USE_EXACTPREDICATES = @load_preference("USE_EXACTPREDICATES", true)::Bool
else 
    const USE_EXACTPREDICATES = true 
end
@doc """
    USE_EXACTPREDICATES

Whether to use ExactPredicates.jl for computing predicates. By default, 
this is true, but a user can change this by defining a preference with Preferences.jl, i.e. 
you could do the following 

```julia-repl
julia> using Preferences: set_preferences!

julia> set_preferences!("DelaunayTriangulation", "USE_EXACTPREDICATES" => false)

julia> using DelaunayTriangulation # load only after setting the preference

julia> DelaunayTriangulation.USE_EXACTPREDICATES
false
```

You have set `USE_EXACTPREDICATES = $USE_EXACTPREDICATES`. 

!!! note "Precision"

    Even if you have disabled ExactPredicates.jl, the predicates 
    are still computed in Float64 precision.
"""
USE_EXACTPREDICATES