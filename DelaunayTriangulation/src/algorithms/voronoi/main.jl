"""
    voronoi(tri::Triangulation; clip=false, smooth=false, kwargs...) -> VoronoiTessellation

Computes the Voronoi tessellation dual to a triangulation. If the triangulation is weighted, this computes the power diagram.

# Arguments
- `tri`: The triangulation.

# Keyword Arguments
- `clip=false`: If `true`, then the Voronoi tessellation is clipped to the convex hull of the triangulation. Otherwise, the Voronoi tessellation is unbounded.
- `clip_polygon=nothing`: If `clip=true`, then this is the polygon to clip the Voronoi tessellation to. If `nothing`, the convex hull of the triangulation is used. The polygon should be defined as a `Tuple` of the form `(points, boundary_nodes)` where the `boundary_nodes` are vertices mapping to coordinates in `points`, adhering to the usual conventions for defining boundaries (in particular, must be counter-clockwise - see https://juliageometry.github.io/DelaunayTriangulation.jl/stable/manual/boundaries/). 

!!! warning "Convex"

    The clipping polygon must be convex. If it is not, unexpected results may occur.

- `smooth=false`: If `true`, then the Voronoi tessellation is smoothed into a centroidal tessellation. Otherwise, the Voronoi tessellation is not smoothed. Must have `clip=true` if `smooth=true`.
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.
- `kwargs...`: Additional keyword arguments passed to [`centroidal_smooth`](@ref) if `smooth=true`.

# Output 
- `vorn`: The [`VoronoiTessellation`](@ref).
"""
function voronoi(tri::Triangulation; clip = false, smooth = false, clip_polygon = nothing, rng = Random.default_rng(), predicates::AbstractPredicateKernel = AdaptiveKernel(), kwargs...)
    if smooth
        @assert clip "Can't pass clip = false when smooth = true. Pass clip = true if you need smoothing."
    end
    if !clip && !isnothing(clip_polygon)
        @warn "Clip polygon provided but clip = false. Ignoring clip polygon."
    end
    has_ghost = has_ghost_triangles(tri)
    !has_ghost && add_ghost_triangles!(tri)
    vorn = initialise_voronoi_tessellation(tri)
    has_bnds = has_boundary_nodes(tri)
    is_true(clip) && !has_bnds && lock_convex_hull!(tri; rng, predicates)
    for i in each_generator(vorn)
        add_voronoi_polygon!(vorn, i)
    end
    if is_true(clip)
        clip_voronoi_tessellation!(vorn; clip_polygon, rng, predicates)
    end
    if is_true(smooth)
        vorn = centroidal_smooth(vorn; clip_polygon, rng, predicates, kwargs...)
    end
    is_true(clip) && !has_bnds && unlock_convex_hull!(tri)
    !has_ghost && delete_ghost_triangles!(tri)
    return vorn
end
