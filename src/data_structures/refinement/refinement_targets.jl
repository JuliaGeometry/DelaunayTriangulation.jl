"""
    RefinementTargets{A,R,P}

A struct containing the user-specified refinement targets. 

# Fields

- `max_area::A`

The maximum area of a triangle. This can also be a function of the form `f(T, p, q, r, A)`, where `T` is the triangle` with area and coordinates `p`, `q`, `r`, returning `true` if the triangle should be refined.

- `max_radius_edge_ratio::R`

The maximum permitted radius-edge ratio. This can also be a function of the form `f(T, p, q, r, ρ)`, where `T` is the triangle` with ratio `ρ` and coordinates `p`, `q`, `r`, returning `true` if the triangle should be refined. Defaults to 1.0, corresponding to a minimum angle of 30°.
- `max_points::P`

The maximum number of points in the mesh. Defaults to Inf.
# Constructors 

The constructor is 

    RefinementTargets(;
            max_area=typemax(Float64),
            max_radius_edge_ratio=nothing,
            max_points=typemax(Int64),
            min_angle = nothing
            )

which allows for a user to specify either the maximum radius-edge ratio `ρ` or the minimum angle `min_angle` in degrees, using the relationship `min_angle = asin[1/(2ρ)]`. If both are provided, the minimum angle is ignored. If neither are provided, the default value of `ρ = 1.0` is used, corresponding to a minimum angle of 30°.

Note that you cannot use `min_angle` as a function, unlike `max_radius_edge_ratio`. If you want to use a function, use `max_radius_edge_ratio` instead.
"""
struct RefinementTargets{A,R,P}
    max_area::A
    max_radius_edge_ratio::R
    max_points::P
    function RefinementTargets(;
        max_area=typemax(Float64),
        max_radius_edge_ratio=nothing,
        max_points=typemax(Int64),
        min_angle=nothing
    )
        if min_angle isa Function
            throw(ArgumentError("Cannot provide min_angle as a function."))
        end
        if max_radius_edge_ratio isa Number && max_radius_edge_ratio < sqrt(3) / 3
            @warn "The provided maximum radius-edge ratio, $max_radius_edge_ratio, is below the theoretical limit of 1/√3. Replacing it with 1/sqrt(3)."
            max_radius_edge_ratio = sqrt(3) / 3
        end
        if max_radius_edge_ratio === nothing && min_angle === nothing
            max_radius_edge_ratio = 1.0
        elseif max_radius_edge_ratio === nothing && min_angle !== nothing
            max_radius_edge_ratio = cscd(min_angle) / 2
        elseif max_radius_edge_ratio !== nothing && min_angle !== nothing
            @warn "Both max_radius_edge_ratio and min_angle are provided. Ignoring min_angle."
        end
        min_angle = max_radius_edge_ratio isa Number ? asind(1 / (2max_radius_edge_ratio)) : nothing
        if !isnothing(min_angle) && (33.9 < min_angle < 60)
            @warn "The provided max_radius_edge_ratio, ρ = $max_radius_edge_ratio, corresponds to a minimum angle of $(min_angle)°. The algorithm may fail to halt with a minimum angle constraint exceeding 33.9° (corresponding to ρ <
             0.9) or higher (meaning lower ρ). Consider using a smaller minimum angle constraint or a larger maximum radius-edge ratio. Will proceed with the provided maximum radius-edge ratio."
        end
        if max_area isa Number && max_area < 0.0
            @warn "The provided maximum area constraint, $max_area, is negative. Replacing with Inf."
            max_area = typemax(Float64)
        end
        if !isnothing(min_angle) && (min_angle < 0.0 || min_angle > 60.00005) # 60.00005 for the =sqrt(3)/3 case
            @warn "The provided max_radius_edge_ratio, ρ = $max_radius_edge_ratio, corresponds to a minimum angle of $(min_angle)°, which is outside the range [0, 60]. Replacing ρ with 1.0, corresponding to a minimum angle of 30°."
            max_radius_edge_ratio = 1.0
        end
        if max_points isa Number && max_points < 0
            @warn "The provided maximum number of points, $max_points, is negative. Replacing with Inf."
            max_points = typemax(Int64)
        end
        return new{typeof(max_area),typeof(max_radius_edge_ratio),typeof(max_points)}(max_area, max_radius_edge_ratio, max_points)
    end
end