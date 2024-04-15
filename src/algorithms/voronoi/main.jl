"""
    voronoi(tri::Triangulation; clip=false, smooth=false, kwargs...) -> VoronoiTessellation

Computes the Voronoi tessellation dual to a triangulation.

# Arguments
- `tri`: The triangulation.

# Keyword Arguments
- `clip=false`: If `true`, then the Voronoi tessellation is clipped to the convex hull of the triangulation. Otherwise, the Voronoi tessellation is unbounded.
- `smooth=false`: If `true`, then the Voronoi tessellation is smoothed into a centroidal tessellation. Otherwise, the Voronoi tessellation is not smoothed. Must have `clip=true` if `smooth=true`.

# Output 
- `vorn`: The [`VoronoiTessellation`](@ref).
"""
function voronoi(tri::Triangulation; clip=false, smooth=false, kwargs...)
    if smooth
        @assert clip "Smoothing requires clipping."
    end
    has_ghost = has_ghost_triangles(tri)
    !has_ghost && add_ghost_triangles!(tri)
    vorn = initialise_voronoi_tessellation(tri)
    has_bnds = has_boundary_nodes(tri)
    is_convex = has_bnds
    is_true(clip) && !has_bnds && lock_convex_hull!(tri)
    for i in each_generator(vorn)
        add_voronoi_polygon!(vorn, i)
    end
    is_true(clip) && clip_voronoi_tessellation!(vorn, is_convex)
    is_true(smooth) && (vorn = centroidal_smooth(vorn; kwargs...))
    is_true(clip) && !has_bnds && unlock_convex_hull!(tri)
    !has_ghost && delete_ghost_triangles!(tri)
    return vorn
end
