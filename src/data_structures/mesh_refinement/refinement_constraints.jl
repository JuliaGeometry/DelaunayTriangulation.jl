"""
    RefinementConstraints{F}

A struct for storing constraints for mesh refinement.

# Fields 
- `min_angle=0.0`: The minimum angle of a triangle.
- `max_angle=180.0`: The maximum angle of a triangle. 
- `min_area=0.0`: The minimum area of a triangle.
- `max_area=Inf`: The maximum area of a triangle.
- `max_radius_edge_ratio=csd(min_angle) / 2`: The maximum radius-edge ratio of a triangle. This is computed from `min_angle` - you cannot provide a value for it yourself. 
- `max_points=typemax(Int)`: The maximum number of vertices allowed in the triangulation. Note that this refers to [`num_solid_vertices`](@ref), not the amount returned by [`num_points`](@ref).
- `seiditous_angle=20.0`: The inter-segment angle used to define seditious edges in degrees. Should not be substantially smaller than 20.0° or any greater than 60.0°.
- `custom_constraint::F=(tri, triangle) -> false`: A custom constraint function. This should take a [`Triangulation`](@ref) and a `triangle` as arguments, and return `true` if the `triangle` violates the constraints and `false` otherwise.
"""
struct RefinementConstraints{F}
    min_angle::Float64
    max_angle::Float64
    min_area::Float64
    max_area::Float64
    max_radius_edge_ratio::Float64
    max_points::Int
    seditious_angle::Float64
    custom_constraint::F
    function RefinementConstraints(;
        min_angle=0.0,
        max_angle=180.0,
        min_area=0.0,
        max_area=Inf,
        max_points=typemax(Int),
        seditious_angle=20.0,
        custom_constraint=(tri, triangle) -> false)
        max_radius_edge_ratio = cscd(min_angle) / 2
        min_angle, max_angle, min_area, max_area, max_radius_edge_ratio, seditious_angle = convert.(Float64, (min_angle, max_angle, min_area, max_area, max_radius_edge_ratio, seditious_angle))
        F = typeof(custom_constraint)
        return new{F}(min_angle, max_angle, min_area, max_area, max_radius_edge_ratio, Int(max_points), seditious_angle, custom_constraint)
    end
end
function Base.show(io::IO, ::MIME"text/plain", constraints::RefinementConstraints)
    println(io, "RefinementConstraints")
    println(io, "   θₘᵢₙ°: ", constraints.min_angle)
    println(io, "   θₘₐₓ°: ", constraints.max_angle)
    println(io, "   Aₘᵢₙ: ", constraints.min_area)
    println(io, "   Aₘₐₓ: ", constraints.max_area)
    println(io, "   ρₘₐₓ: ", constraints.max_radius_edge_ratio)
    println(io, "   max_points: ", constraints.max_points == typemax(Int) ? "∞" : constraints.max_points)
    println(io, "   θ_seditious°: ", constraints.seditious_angle)
    print(io, "   custom_constraint: ", constraints.custom_constraint)
end

"""
    has_max_angle_constraint(constraints::RefinementConstraints) -> Bool

Return `true` if `constraints` has a maximum angle constraint, `false` otherwise.
"""
has_max_angle_constraint(constraints::RefinementConstraints) = constraints.max_angle < 180.0

"""
    violates_custom_constraint(constraints::RefinementConstraints{F}, tri::Triangulation, T) where {F} -> Bool

Return `true` if `T` violates the custom constraint in `constraints`, `false` otherwise.
"""
violates_custom_constraint(constraints::RefinementConstraints{F}, tri::Triangulation, T) where {F} = constraints.custom_constraint(tri, T)
