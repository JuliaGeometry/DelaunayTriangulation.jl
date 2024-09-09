using DelaunayTriangulation, Random, GeometryBasics, GLMakie
GLMakie.activate!()
const DT = DelaunayTriangulation
struct UnitSphere{T} end
UnitSphere() = UnitSphere{Float64}()
struct SphericalPoint{T} # assumed to be on the unit sphere
    x::T # could also parametrise using latitude and longitude
    y::T
    z::T
end
Random.eltype(::Type{UnitSphere{T}}) where {T} = SphericalPoint{T}
function Random.rand(rng::AbstractRNG, ::Random.SamplerTrivial{UnitSphere{T}}) where {T}
    x, y, z = ntuple(_ -> randn(rng, T), Val(3))
    r⁻¹ = 1 / sqrt(x^2 + y^2 + z^2)
    return SphericalPoint(x * r⁻¹, y * r⁻¹, z * r⁻¹)
end
GeometryBasics.Point3(p::SphericalPoint) = Point3(p.x, p.y, p.z)
project(p::SphericalPoint) = (projectx(p), projecty(p))
function invproject(p)
    X, Y = getxy(p)
    den = 1 + X^2 + Y^2
    return SphericalPoint(2X / den, 2Y / den, (X^2 + Y^2 - 1) / den)
end
projectx(p::SphericalPoint) = p.x / (1 - p.z)
projecty(p::SphericalPoint) = p.y / (1 - p.z)
DT.getx(p::SphericalPoint) = projectx(p)
DT.gety(p::SphericalPoint) = projecty(p)
DT.number_type(::Type{SphericalPoint{T}}) where {T} = T

function DT.triangle_circumcenter(tri::Triangulation{Vector{SphericalPoint{T}}}, V) where {T} # If we used a SphericalTriangle type instead, we'd dispatch on that
    points = get_points(tri)
    u, v, w = triangle_vertices(V)
    p, q, r = points[u], points[v], points[w] # can't use get_point because it will return the projection instead
    return spherical_circumcenter(p, q, r)
end
function spherical_circumcenter(p::SphericalPoint, q::SphericalPoint, r::SphericalPoint)
    _p, _q, _r = Point3(p), Point3(q), Point3(r)
    cross_pq = Point3(_p[2] * _q[3] - _q[2] * _p[3], _q[1] * _p[3] - _p[1] * _q[3], _p[1] * _q[2] - _q[1] * _p[2])
    cross_qr = Point3(_q[2] * _r[3] - _r[2] * _q[3], _r[1] * _q[3] - _q[1] * _r[3], _q[1] * _r[2] - _r[1] * _q[2])
    cross_rp = Point3(_r[2] * _p[3] - _p[2] * _r[3], _p[1] * _r[3] - _r[1] * _p[3], _r[1] * _p[2] - _p[1] * _r[2])
    num = cross_pq + cross_qr + cross_rp
    den = sqrt(num[1]^2 + num[2]^2 + num[3]^2)
    return project(SphericalPoint(num[1] / den, num[2] / den, num[3] / den)) # project because DelaunayTriangulation uses NTuple{2}
end

function spherical_triangulate(points; kwargs...)
    idx = findfirst(==(SphericalPoint(0.0, 0.0, 1.0)), DT.each_point(points))
    skip_points = isnothing(idx) ? 0 : idx
    tri = triangulate(points; skip_points = skip_points, kwargs...) # assuming skip_points not in kwargs
    # We could just note that the missing part of the triangulation near the north pole 
    # can be represented by each_ghost_triangle(tri), mapping -1 to the north pole. To use Mesh, 
    # though, this won't work. So, let's put the triangles into each_solid_triangle and add the 
    # north pole to points. 
    if isnothing(idx)
        push!(points, SphericalPoint(0.0, 0.0, 1.0))
        n = length(points)
    else
        n = idx 
    end
    for T in each_ghost_triangle(tri)
        u, v, w = triangle_vertices(DT.sort_triangle(T))
        u′, v′, w′ = DT.sort_triangle(u, v, w)
        add_triangle!(tri, v′, u′, n)
    end
    # delete_ghost_triangles!(tri)
    return GeometryBasics.Mesh(
        map(Point3, points),
        map(TriangleFace, each_solid_triangle(tri))
    )
end

function spherical_voronoi(points; kwargs...)
    tri = triangulate(points; kwargs...)
    vorn = voronoi(tri)
    # Now we need to put in the 
end

function line_between_points_on_sphere(p, q)
    pqorth = Point3(
        -p[2] * (p[1] * q[2] - q[1] * p[2]) - p[3] * (p[1] * q[3] - q[1] * p[3]),
        p[1] * (p[1] * q[2] - q[1] * p[2]) - p[3] * (p[2] * q[3] - q[2] * p[3]),
        p[1] * (p[1] * q[3] - q[1] * p[3]) + p[2] * (p[2] * q[3] - q[2] * p[3])
    )
    pqorthnorm = sqrt(pqorth[1]^2 + pqorth[2]^2 + pqorth[3]^2)
    pqorth = pqorth / pqorthnorm
    pqcross = Point3(
        p[2] * q[3] - q[2] * p[3],
        q[1] * p[3] - p[1] * q[3],
        p[1] * q[2] - q[1] * p[2]
    )
    pqcrossnorm = sqrt(pqcross[1]^2 + pqcross[2]^2 + pqcross[3]^2)
    pqdot = p[1] * q[1] + p[2] * q[2] + p[3] * q[3]
    ϕ = mod2pi(atan(pqcrossnorm, pqdot))
    t = LinRange(0, ϕ, 50)
    points = map(t) do τ
        s, c = sincos(τ)
        return p * c + pqorth * s
    end
    return points
end
function to_spherical_triangles(tri)
    points = tri.position
    triangles = faces(tri)
    spherical_triangles = Point3d[]
    foreach(triangles) do T
        i, j, k = triangle_vertices(T)
        p, q, r = points[i], points[j], points[k]
        L1 = line_between_points_on_sphere(p, q)
        L2 = line_between_points_on_sphere(q, r)
        L3 = line_between_points_on_sphere(r, p)
        for L in (L1, L2, L3)
            append!(spherical_triangles, L)
            push!(spherical_triangles, Point3d(NaN, NaN, NaN))
        end
    end
    return spherical_triangles
end

points = rand(UnitSphere(), 100)
tri = spherical_triangulate(points)
fig, ax, sc = poly(tri, color=:blue, strokewidth=1, strokecolor=:black)
scatter!(ax, [Point3(0.0, 0.0, 1.0)], color=:red, markersize=14)
spherical_triangles = to_spherical_triangles(tri)
fig, ax, sc = poly(Sphere(Point3(0.0, 0.0, 0.0), 1.0))
lines!(ax, spherical_triangles, color=:red, linewidth=1)
fig