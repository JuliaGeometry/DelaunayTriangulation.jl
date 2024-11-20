@inline function get_plane_through_three_points(a::Tuple, b::Tuple, c::Tuple)
    ax, ay, az = a
    bx, by, bz = b
    cx, cy, cz = c
    α = ay * bz - az * by - ay * cz + az * cy + by * cz - bz * cy
    β = az * bx - ax * bz + ax * cz - az * cx - bx * cz + bz * cx
    γ = ax * by - ay * bx - ax * cy + ay * cx + bx * cy - by * cx
    δ = ax * bz * cy - ax * by * cz + ay * bx * cz - ay * bz * cx - az * bx * cy + az * by * cx
    return α, β, γ, δ
end

function get_steepest_descent_direction(a::Tuple, b::Tuple, c::Tuple)
    α, β, γ, _ = get_plane_through_three_points(a, b, c)
    sγ = sign(γ)
    return sγ * α, sγ * β
end

function get_vertical_distance_to_plane(a::Tuple, b::Tuple, c::Tuple, p::Tuple)
    α, β, γ, δ = get_plane_through_three_points(a, b, c)
    x, y, z = p
    plane_z = -(α * x + β * y + δ) / γ
    return z - plane_z
end

function get_distance_to_plane(a::Tuple, b::Tuple, c::Tuple, p::Tuple)
    α, β, γ, δ = get_plane_through_three_points(a, b, c)
    if γ < 0 # so that nz > 0
        α, β, γ, δ = -α, -β, -γ, -δ
    end
    nmag = sqrt(α^2 + β^2 + γ^2)
    nx, ny, nz, Jd = α / nmag, β / nmag, γ / nmag, δ / nmag
    x, y, z = p
    return x * nx + y * ny + z * nz + Jd
end