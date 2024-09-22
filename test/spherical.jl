using DelaunayTriangulation
using DelaunayTriangulation.SphericalDelaunay
using Test
using Makie, CairoMakie, GLMakie
CairoMakie.activate!()
using Random
using LinearAlgebra

const DT = DelaunayTriangulation
const SD = SphericalDelaunay

Makie.to_vertices(p::Vector{<:SphericalPoint}) = Makie.Point3.(SD.__getxyz.(p))
Makie.convert_arguments(p::PointBased, points::Vector{<:SphericalPoint}) = Makie.convert_arguments(p, Makie.to_vertices(points))
Makie.convert_arguments(p::PointBased, L::SphericalLine) = Makie.convert_arguments(p, SD._to_coords(L))
Makie.convert_arguments(p::PointBased, T::SphericalTriangle) = Makie.convert_arguments(p, SD._to_coords(T))
Makie.convert_arguments(p::PointBased, P::SphericalPolygon) = Makie.convert_arguments(p, SD._to_coords(P))
function Makie.convert_arguments(p::Type{<:Lines}, tri::SphericalTriangulation)
    points = SD._to_lines(tri)
    return Makie.convert_arguments(p, points)
end
function Makie.convert_arguments(p::Type{<:Lines}, vorn::SphericalTessellation)
    points = SD._to_lines(vorn)
    return Makie.convert_arguments(p, points)
end
# These mesh recipes render SO poorly...
function Makie.convert_arguments(p::Type{<:Mesh}, T::SphericalTriangle)
    points = SD._to_coords(T, false)
    faces = SD._to_mesh(points)
    return Makie.convert_arguments(p, points, faces)
end
function Makie.convert_arguments(p::Type{<:Mesh}, tri::SphericalTriangulation)
    points, vertices = SD._to_triangles(tri)
    return Makie.convert_arguments(p, points, vertices)
end
function Makie.convert_arguments(p::Type{<:Poly}, tri::SphericalTriangulation)
    res = Makie.convert_arguments(Mesh, tri)
    return Makie.convert_arguments(p, res...)
end

function uniform_rand_sphere(N)
    # http://web.archive.org/web/20120421191837/http://www.cgafaq.info/wiki/Evenly_distributed_points_on_sphere
    coords = Vector{SphericalPoint{NTuple{3,Float64}}}(undef, N)
    s = 3.6 / sqrt(N)
    dz = 2 / N
    z = 1 - dz / 2
    long = 0.0
    for k in 1:N
        r = sqrt(1 - z^2)
        latitude = asind(z)
        longitude = rad2deg(long)
        coords[k] = SphericalPoint(latitude, longitude % 360)
        long += s / r
        z -= dz
    end
    return coords
end

@testset "Basic" begin
    @test SD.UnitSphere() == SD.UnitSphere{Float64}()
    p = (1.0, 0.0, 0.0)
    @test SphericalPoint(p) == SphericalPoint{typeof(p)}(p)
    @test SphericalPoint(p) == SphericalPoint(p)
    @test SphericalPoint(0.3, 0.2, sqrt(1 - 0.3^2 - 0.2^2)) == SphericalPoint((0.3, 0.2, sqrt(1 - 0.3^2 - 0.2^2)))
    q = (0.0, 1.0, 0.0)
    @test SphericalPoint(p) != SphericalPoint(q)
    p = rand(UnitSphere())
    @test p isa SphericalPoint
    @test sum(p.p .^ 2) ≈ 1.0
    @test SD.getp(p) == p.p
    @test SD.__getxyz(p) == p.p
    x, y, z = p.p
    @test SD.projectx(p) == getx(p) == x / (1 - z)
    @test SD.projecty(p) == gety(p) == y / (1 - z)
    @test SD.project(p) == getxy(p)
    X, Y = SD.project(p)
    @test collect(SD.invproject((X, Y)).p) ≈ collect(p.p)
    @test DT.number_type(p) == Float64
    @test DT.number_type(rand(UnitSphere{Float32}())) == Float32
end

@testset "SphericalPoints" begin
    points = rand(UnitSphere(), 25)
    spoints = SD.SphericalPoints(points)
    @test SD.SphericalPoints(spoints) == spoints
    @test eltype(points) == SphericalPoint{NTuple{3,Float64}}
    @test eachindex(spoints) == eachindex(points)
    @test iterate(spoints) == iterate(points)
    @test collect(spoints) == collect(points)
    @test length(spoints) == length(points)
    @test size(spoints) == size(points)
    @test points[2] == spoints[2]
    @test DT.number_type(spoints) == Float64
    DT.set_point!(spoints, 7, SD.project(SphericalPoint(0.3, 0.2, sqrt(1 - 0.3^2 - 0.2^2)))...)
    @test collect(spoints[7].p) ≈ collect(SphericalPoint(0.3, 0.2, sqrt(1 - 0.3^2 - 0.2^2)).p)
    val2 = spoints[end]
    val = pop!(spoints)
    @test val == val2
    @test length(spoints) == length(points) == 24
    DT.push_point!(spoints, SD.project(SphericalPoint(0.3, 0.2, sqrt(1 - 0.3^2 - 0.2^2)))...)
    @test collect(spoints[end].p) ≈ collect(SphericalPoint(0.3, 0.2, sqrt(1 - 0.3^2 - 0.2^2)).p)
    @test SD._getindex(spoints, -1) == SD.north_pole(Float64)
end

function spherical_distance(p, q)
    a1, a2, a3 = p
    b1, b2, b3 = q
    return acos(a1 * b1 + a2 * b2 + a3 * b3)
end
_cross(p, q) = cross(collect(p), collect(q))
@testset "Spherical Triangles" begin
    for _ in 1:10000
        p, q = rand(UnitSphere(), 2)
        c = cross(collect(p.p), collect(q.p))
        @test collect(c) ≈ collect(SD.cross(p, q))
        p, q, r = rand(UnitSphere(), 3)
        sc = SD.spherical_circumcenter(p, q, r)
        dists = spherical_distance.((sc.p,), (p.p, q.p, r.p))
        @test all(dists .≈ dists[1])
        p, q = rand(UnitSphere(), 2)
        d = dot(collect(p.p), collect(q.p))
        @test d ≈ SD.dot(p, q)
        u, v, w = rand(3)
        @test SD._hypot(u, v, w) ≈ sqrt(u^2 + v^2 + w^2)
        p, q, r = SphericalPoint(30.2672, 97.7431), SphericalPoint(42.3584, 71.0598), SphericalPoint(51.0501, 114.0853)
        A = SD.spherical_triangle_area(p, q, r)
        @test A ≈ 0.090030685
        q, r = SphericalPoint(29.7947, 98.7320), SphericalPoint(32.9746, 96.8899)
        A = SD.spherical_triangle_area(p, q, r)
        @test A ≈ 0.0003032270459
        p, q, r = rand(UnitSphere(), 3)
        s = SD.spherical_triangle_circumradius(p, q, r)
        sc = SD.spherical_circumcenter(p, q, r)
        dists = spherical_distance.((sc.p,), (p.p, q.p, r.p))
        @test all(dists .≈ s)
        p, q = rand(UnitSphere(), 2)
        @test SD.spherical_distance(p, q) ≈ spherical_distance(p.p, q.p)
        p, q, r = rand(UnitSphere(), 3)
        stp = SD.scalar_triple_product(p, q, r)
        @test stp ≈ dot(collect(p.p), cross(collect(q.p), collect(r.p)))

        p, q, r = rand(UnitSphere(), 3)
        if SD.scalar_triple_product(p, q, r) < 0
            p, q, r = q, p, r
        end
        A, B, C = SD.spherical_angles(p, q, r)
        @test sin(A) / norm(cross(collect(q.p), collect(r.p))) ≈
              sin(B) / norm(cross(collect(r.p), collect(p.p))) ≈
              sin(C) / norm(cross(collect(p.p), collect(q.p))) ≈
              SD.scalar_triple_product(p, q, r) / (norm(cross(collect(p.p), collect(q.p))) * norm(cross(collect(q.p), collect(r.p))) * norm(cross(collect(r.p), collect(p.p))))
        excess = A + B + C - π
        @test SD.spherical_triangle_area(p, q, r) ≈ excess

        p, q, r = rand(UnitSphere(), 3)
        if SD.scalar_triple_product(p, q, r) < 0
            p, q, r = q, p, r
        end
        c = SD.spherical_triangle_centroid(p, q, r)
        @test norm(collect(c.p)) ≈ 1.0
        A1 = SD.spherical_triangle_area(p, q, c)
        A2 = SD.spherical_triangle_area(q, r, c)
        A3 = SD.spherical_triangle_area(r, p, c)
        @test A1 ≈ A2 ≈ A3 ≈ SD.spherical_triangle_area(p, q, r) / 3
    end
    p, q, r = rand(UnitSphere(), 3)
    @inferred SD.spherical_triangle_circumradius(p, q, r)
    @inferred SD.spherical_triangle_area(p, q, r)
    @inferred SD.spherical_circumcenter(p, q, r)
    @inferred SD.cross(p, q)
    @inferred SD.dot(p, q)
    @inferred SD._hypot(rand(3)...)
    @inferred SD.projectx(p)
    @inferred SD.projecty(p)
    @inferred SD.project(p)
    X, Y = SD.project(rand(UnitSphere()))
    @inferred SD.invproject((X, Y))
    @inferred SD.getp(p)
    @inferred SD.__getxyz(p)
    @test SD.north_pole(Float64) == SphericalPoint(0.0, 0.0, 1.0)
    @test SD.north_pole(Float32) == SphericalPoint(0.0, 0.0, 1.0f0)
    @test SD.invproject((NaN, NaN)) == SphericalPoint(0.0, 0.0, 1.0)
end

@testset "SphericalLine/Triangle" begin
    for _ in 1:100
        p, q = rand(UnitSphere(), 2)
        L = SD.SphericalLine(p, q)
        t = LinRange(0.0, L.upper, 500)
        pts = L.(t)
        _pts = [collect(x.p) for x in pts]
        dist = sum(norm.(diff(_pts)))
        @test dist ≈ spherical_distance(p.p, q.p) rtol = 1e-4

        p, q, r = rand(UnitSphere(), 3)
        L1, L2, L3 = SphericalLine(p, q), SphericalLine(q, r), SphericalLine(r, p)
        T = SD.SphericalTriangle(L1, L2, L3)
        T2 = SD.SphericalTriangle(p, q, r)
        @test T == T2
    end
end

points = uniform_rand_sphere(35)
tri = spherical_triangulate(points)