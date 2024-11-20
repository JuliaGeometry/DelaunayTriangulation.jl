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
            min_angle = 0.0,
            max_angle = 180.0,
            min_area = 0.0,
            max_area = Inf,
            max_points = typemax(Int),
            seditious_angle = 20.0,
            custom_constraint = (tri, triangle) -> false,
        )
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

@inline has_max_angle_constraint(constraints::RefinementConstraints) = constraints.max_angle < 180.0
@inline violates_custom_constraint(constraints::RefinementConstraints{F}, tri::Triangulation, T) where {F} = constraints.custom_constraint(tri, T)
