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

Alias for [`GhostVertex`](@ref).
"""
const 𝒢 = GhostVertex

"""
    ε(x) = sqrt(eps(number_type(x)))

Number used as a tolerance in certain functions, e.g. 
for mesh refinement when using [`check_precision`](@ref) to 
avoid degenerate circumcenters.
"""
ε(::Type{T}) where {T} = sqrt(eps(T))
ε(x) = ε(number_type(x))

const INF_WARN = Ref(true)
"""
    toggle_inf_warn!()

Toggle the warning for infinite circumcenters in the Voronoi tessellation.
By default, this warning is enabled.
"""
toggle_inf_warn!() = (INF_WARN[] = !INF_WARN[])

const WARN_ON_DUPES = Ref(true)
"""
    toggle_warn_on_dupes!()

Toggle the warning for duplicate points in the input data.
By default, this warning is enabled.
"""
toggle_warn_on_dupes!() = (WARN_ON_DUPES[] = !WARN_ON_DUPES[])

