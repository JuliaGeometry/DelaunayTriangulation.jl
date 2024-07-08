using Test, DataStructures
using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using LinearAlgebra
EllipticalArc = DT.EllipticalArc # shadow
using CairoMakie
using StableRNGs
using Preferences
using StructEquality
using ForwardDiff
using ReferenceTests

@testset "LineSegment" begin
    # Construction 
    p = rand(2) |> Tuple
    q = rand(2) |> Tuple
    L = LineSegment(p, q)
    @test !DT.is_interpolating(L)
    @test DT.is_curve_bounded(L)
    @test !DT.is_piecewise_linear(L)
    @test L.first == p
    @test L.last == q
    @test L.length ≈ sqrt(sum((p .- q) .^ 2))

    # Evaluation 
    @test collect(L(0.0)) ≈ collect(p)
    @test collect(L(1.0)) ≈ collect(q)
    @test collect(L(0.398881)) ≈ collect(p .+ 0.398881 .* (q .- p))

    # Sideof 
    p = (1.0, 1.0)
    q = (6.0, 1.0)
    L = LineSegment(p, q)
    revL = LineSegment(q, p)
    d = (3.0, 3.0)
    cert1 = DT.point_position_relative_to_curve(L, d)
    cert2 = DT.point_position_relative_to_curve(revL, d)
    @test DT.is_left(cert1) && DT.is_right(cert2)
    d = (3.0, 1.0)
    cert1 = DT.point_position_relative_to_curve(L, d)
    cert2 = DT.point_position_relative_to_curve(revL, d)
    @test DT.is_on(cert1) && DT.is_on(cert2)
    d = (3.0, 0.0)
    cert1 = DT.point_position_relative_to_curve(L, d)
    cert2 = DT.point_position_relative_to_curve(revL, d)

    ## Arclength 
    p = rand(2) |> Tuple
    q = rand(2) |> Tuple
    L = LineSegment(p, q)
    s = DT.arc_length(L)
    @test s ≈ norm(p .- q)
    s = DT.arc_length(L, 0.3, 0.7)
    @test s ≈ norm(p .- q) * 0.4
    @test s ≈ norm(p .+ 0.3 .* (q .- p) .- (p .+ 0.7 .* (q .- p)))
    @test s ≈ slow_arc_length(L, 0.3, 0.7)

    ## Differentiate 
    p = rand(2) |> Tuple
    q = rand(2) |> Tuple
    L = LineSegment(p, q)
    for t in LinRange(0, 1, 100)
        der1 = DT.differentiate(L, t)
        h = 1e-8
        der2 = (L(t + h) .- L(t - h)) ./ (2h)
        @test der1 ⪧ der2 rtol = 1e-5 atol = 1e-5
    end
    @test DT.twice_differentiate(L, rand()) == (0.0, 0.0)
    @inferred DT.differentiate(L, rand())
    @inferred DT.twice_differentiate(L, rand())

    ## Curvature 
    @test DT.curvature(L, rand()) ≈ 0.0

    ## Total variation 
    TV = DT.total_variation(L)
    @test TV ≈ 0.0 atol = 1e-6
    TV = DT.total_variation(L, 0.2, 0.5)
    @test TV ≈ 0.0 atol = 1e-6
    @test TV ≈ slow_total_absolute_curvature(L, 0.2, 0.5) atol = 1e-3

    ## Equidistant split 
    t = DT.get_equidistant_split(L, 0, 1)
    @test t ≈ 1 / 2
    t = DT.get_equidistant_split(L, 0.38, 0.95)
    @test t ≈ (0.38 + 0.95) / 2

    ## Equivariation split 
    t, T = DT.get_equivariation_split(L, 0, 1)
    @test t ≈ 1 / 2 && T == 0.0
    t, T = DT.get_equivariation_split(L, 0.38, 0.95)
    @test t ≈ (0.38 + 0.95) / 2 && T == 0.0

    ## Inverse 
    L1 = LineSegment(Tuple(rand(2)), Tuple(rand(2)))
    L2 = LineSegment((2.0, 5.0), (2.0, 10.0))
    L3 = LineSegment((2.0, 5.0), (-5.0, 5.0))
    for t in LinRange(0, 1, 2500)
        for L in (L1, L2, L3)
            @test t ≈ DT.get_inverse(L, L(t))
            @inferred DT.get_inverse(L, L(t))
        end
    end

    ## get_circle_intersection
    p, q = (0.0, 0.0), (1.0, 1.0)
    L = LineSegment(p, q)
    r = 2L.length / 3
    t′, q′ = DT.get_circle_intersection(L, 0.0, 1.0, r)
    @test t′ ≈ 2 / 3 && q′ ⪧ L(t′) ⪧ (2 / 3, 2 / 3)
    @inferred DT.get_circle_intersection(L, 0.0, 1.0, r)
    @test DT.dist(L(0.0), q′) ≈ r

    r = 0.1
    t′, q′ = DT.get_circle_intersection(L, 0.2, 1.0, r)
    @test t′ ≈ 0.08838834764831845 && q′ ⪧ (0.2707106781186548, 0.2707106781186548)
    @test DT.dist(L(0.2), q′) ≈ r

    r = 0.5
    t′, q′ = DT.get_circle_intersection(L, 1.0, 0.0, r)
    @test t′ ≈ 0.6464466094067263 && q′ ⪧ (0.6464466094067263, 0.6464466094067263)

    # == 
    p = rand(2) |> Tuple
    q = rand(2) |> Tuple
    L1 = LineSegment(p, q)
    L2 = LineSegment(p, q)
    @test L1 == L2
    L2 = LineSegment(q, p)
    @test L1 ≠ L2
    L2 = LineSegment(p, p)
    @test L1 ≠ L2
    L2 = LineSegment(q, q)
    @test L1 ≠ L2
end

@testset "PiecewiseLinear" begin
    L = DT.PiecewiseLinear([(1.0, 0.0), (1.0, 1.0)], [1, 1, 1])
    @test !DT.is_curve_bounded(L)
    @test DT.is_piecewise_linear(L)
    @test !DT.is_interpolating(L)
    @test get_points(L) == [(1.0, 0.0), (1.0, 1.0)]
    @test get_boundary_nodes(L) == [1, 1, 1]

    L = DT.PiecewiseLinear(rand(2, 50), [5, 17, 23, 19, 25])
    @test L(0) == L.points[:, 5] |> Tuple
    @test L(1) == L.points[:, 25] |> Tuple
    @inferred L(0)
    @inferred L(1)
    @inferred L(1 / 2)
    ∂ = DT.differentiate(L, 0.0)
    @test ∂ ⪧ DT.differentiate(LineSegment(L.points[:, 5] |> Tuple, L.points[:, 17] |> Tuple), 0.0)
    @inferred DT.differentiate(L, 0.0)
    ∂ = DT.differentiate(L, 1.0)
    @test ∂ ⪧ DT.differentiate(LineSegment(L.points[:, 19] |> Tuple, L.points[:, 25] |> Tuple), 0.0)

    points = [(-1.0, 1.0), (0.0, 0.0), (0.5, 0.2), (0.7, 1.0), (1.0, 2.0), (1.2, 0.5)]
    boundary_nodes = [1, 2, 3, 4, 5, 6]
    L = DT.PiecewiseLinear(points, boundary_nodes)
    r = 0.5
    t′, q′ = DT.get_circle_intersection(L, 0.0, 1.0, r)
    @test q′ ⪧ (-0.6464466094067263, 0.6464466094067263)
    @inferred DT.get_circle_intersection(L, 0.0, 1.0, r)
    @test DT.dist(L(0.0), q′) ≈ r
    r = 1.3
    t′, q′ = DT.get_circle_intersection(L, 0.0, 1.0, r)
    @test q′ ⪧ (-0.08076118445748826, 0.08076118445748826)
    @test DT.dist(L(0.0), q′) ≈ r
    r = 0.3
    t′, q′ = DT.get_circle_intersection(L, 1.0, 0.0, r)
    @test q′ ⪧ (1.1603508839726946, 0.7973683702047905)
    @test DT.dist(L(1.0), q′) ≈ r
    @test DT.PiecewiseLinear(points, boundary_nodes) == DT.PiecewiseLinear(points, boundary_nodes)
    @test DT.PiecewiseLinear(points, boundary_nodes) ≠ DT.PiecewiseLinear(points, [1, 2, 3, 4, 5, 6, 7])
    @test DT.PiecewiseLinear(points, boundary_nodes) ≠ DT.PiecewiseLinear(points, boundary_nodes[2:end])
    @test DT.PiecewiseLinear(points, boundary_nodes) ≠ DT.PiecewiseLinear(points[2:end], boundary_nodes)
end

@testset "CircularArc" begin
    ## Construction 
    center1 = (1.0, 0.0)
    p1 = (5.0, 5.0)
    q1 = (-4.0, -4.0)
    arc = CircularArc(p1, q1, center1)
    @test DT.is_curve_bounded(arc)
    @test !DT.is_piecewise_linear(arc)
    @test !DT.is_interpolating(arc)
    @inferred CircularArc(p1, q1, center1)
    negarc = CircularArc(p1, q1, center1, positive=false)
    revarc = CircularArc(q1, p1, center1)
    revnegarc = CircularArc(q1, p1, center1, positive=false)
    center2 = (3.0, -5.0)
    p2 = (6.0, 5.0)
    circ = CircularArc(p2, p2, center2)
    revcirc = CircularArc(p2, p2, center2, positive=false)
    @test arc.center == negarc.center == revarc.center == revnegarc.center == center1
    @test circ.center == revcirc.center == center2
    @test arc.radius == negarc.radius == revarc.radius == revnegarc.radius ≈ sqrt(41)
    @test circ.radius == revcirc.radius ≈ sqrt(109)
    @test arc.start_angle ≈ 0.8960554793197029 rtol = 1e-6
    @test negarc.start_angle ≈ 0.8960554793197029 rtol = 1e-6
    @test arc.sector_angle ≈ deg2rad(167.3196155081802) rtol = 1e-6
    @test negarc.sector_angle ≈ deg2rad(167.3196155081802) - 2π rtol = 1e-6
    @test arc.first == negarc.first == revarc.last == revnegarc.last == p1
    @test arc.last == negarc.last == revarc.first == revnegarc.first == q1
    @test circ.first == revcirc.first == circ.last == revcirc.last == p2
    @test arc.pqr ⪧ ((1 + sqrt(41), 0.0), (1.0, sqrt(41.0)), (1.0 - sqrt(41.0), 0.0))
    @test negarc.pqr ⪧ ((1 - sqrt(41), 0.0), (1.0, sqrt(41.0)), (1.0 + sqrt(41.0), 0.0))
    @test revarc.pqr ⪧ ((1 + sqrt(41), 0.0), (1.0, sqrt(41.0)), (1.0 - sqrt(41.0), 0.0))
    @test revnegarc.pqr ⪧ ((1 - sqrt(41), 0.0), (1.0, sqrt(41.0)), (1.0 + sqrt(41.0), 0.0))

    ## Evaluation 
    @test arc(0.0) ⪧ negarc(0.0) ⪧ revarc(1.0) ⪧ revnegarc(1.0) ⪧ p1
    @test arc(1.0) ⪧ negarc(1.0) ⪧ revarc(0.0) ⪧ revnegarc(0.0) ⪧ q1
    @test arc(1e-16) ⪧ negarc(1e-16) ⪧ revarc(1 - 1e-16) ⪧ revnegarc(1 - 1e-16) ⪧ p1
    @test arc(1 - 1e-16) ⪧ negarc(1 - 1e-16) ⪧ revarc(1e-16) ⪧ revnegarc(1e-16) ⪧ q1

    t = LinRange(0, 1, 500)
    θ₀_arc, θ₁_arc = arc.start_angle, arc.start_angle + arc.sector_angle
    θ₀_negarc, θ₁_negarc = negarc.start_angle, negarc.start_angle + negarc.sector_angle
    θ₀_revarc, θ₁_revarc = revarc.start_angle, revarc.start_angle + revarc.sector_angle
    θ₀_revnegarc, θ₁_revnegarc = revnegarc.start_angle, revnegarc.start_angle + revnegarc.sector_angle
    θ₀_circ, θ₁_circ = circ.start_angle, circ.start_angle + circ.sector_angle
    θ₀_revcirc, θ₁_revcirc = revcirc.start_angle, revcirc.start_angle + revcirc.sector_angle
    for t in t
        θarc = (θ₁_arc - θ₀_arc) * t + θ₀_arc
        θnegarc = (θ₁_negarc - θ₀_negarc) * t + θ₀_negarc
        θrevarc = (θ₁_revarc - θ₀_revarc) * t + θ₀_revarc
        θrevnegarc = (θ₁_revnegarc - θ₀_revnegarc) * t + θ₀_revnegarc
        θcirc = (θ₁_circ - θ₀_circ) * t + θ₀_circ
        θrevcirc = (θ₁_revcirc - θ₀_revcirc) * t + θ₀_revcirc
        @test arc(t) ⪧ (1.0 + sqrt(41.0) * cos(θarc), sqrt(41.0) * sin(θarc))
        @test negarc(t) ⪧ (1.0 + sqrt(41.0) * cos(θnegarc), sqrt(41.0) * sin(θnegarc))
        @test revarc(t) ⪧ (1.0 + sqrt(41.0) * cos(θrevarc), sqrt(41.0) * sin(θrevarc))
        @test revnegarc(t) ⪧ (1.0 + sqrt(41.0) * cos(θrevnegarc), sqrt(41.0) * sin(θrevnegarc))
        @test circ(t) ⪧ (3.0 + sqrt(109.0) * cos(θcirc), -5.0 + sqrt(109.0) * sin(θcirc))
        @test revcirc(t) ⪧ (3.0 + sqrt(109.0) * cos(θrevcirc), -5.0 + sqrt(109.0) * sin(θrevcirc))
        @inferred arc(t)
    end

    ## Sideof 
    d = (-7.0, 5.0) # outside
    certarc = DT.point_position_relative_to_curve(arc, d)
    certnegarc = DT.point_position_relative_to_curve(negarc, d)
    certrevarc = DT.point_position_relative_to_curve(revarc, d)
    certrevnegarc = DT.point_position_relative_to_curve(revnegarc, d)
    @test DT.is_right(certarc) && DT.is_left(certnegarc) && DT.is_right(certrevarc) && DT.is_left(certrevnegarc)

    d = (-1.0, 3.0) # inside
    certarc = DT.point_position_relative_to_curve(arc, d)
    certnegarc = DT.point_position_relative_to_curve(negarc, d)
    certrevarc = DT.point_position_relative_to_curve(revarc, d)
    certrevnegarc = DT.point_position_relative_to_curve(revnegarc, d)
    @test DT.is_left(certarc) && DT.is_right(certnegarc) && DT.is_left(certrevarc) && DT.is_right(certrevnegarc)

    d = (1.0 + sqrt(41), 0.0) # on 
    certarc = DT.point_position_relative_to_curve(arc, d)
    certnegarc = DT.point_position_relative_to_curve(negarc, d)
    certrevarc = DT.point_position_relative_to_curve(revarc, d)
    certrevnegarc = DT.point_position_relative_to_curve(revnegarc, d)
    @test DT.is_on(certarc) && DT.is_on(certnegarc) && DT.is_on(certrevarc) && DT.is_on(certrevnegarc)

    d = (-11.073, 5.1257) # outside
    certcirc = DT.point_position_relative_to_curve(circ, d)
    certrevcirc = DT.point_position_relative_to_curve(revcirc, d)
    @test DT.is_right(certcirc) && DT.is_left(certrevcirc)

    d = (-4.0, 0.0) # inside 
    certcirc = DT.point_position_relative_to_curve(circ, d)
    certrevcirc = DT.point_position_relative_to_curve(revcirc, d)
    @test DT.is_left(certcirc) && DT.is_right(certrevcirc)

    d = (3.0 + sqrt(109.0), -5.0) # on
    certcirc = DT.point_position_relative_to_curve(circ, d)
    certrevcirc = DT.point_position_relative_to_curve(revcirc, d)
    @test DT.is_on(certcirc) && DT.is_on(certrevcirc)

    ## Arclength
    for c in (arc, negarc, revarc, revnegarc, circ, revcirc)
        @test DT.arc_length(c) ≈ slow_arc_length(c, 0, 1)
        @test DT.arc_length(c, 0.1, 0.15) ≈ slow_arc_length(c, 0.1, 0.15)
        @test DT.arc_length(c, 0.0, 0.0) ≈ 0.0
        @test DT.arc_length(c, 1.0, 1.0) ≈ 0.0
        @test DT.arc_length(c, 0.0, 1.0) ≈ slow_arc_length(c, 0.0, 1.0)
        @test DT.arc_length(c, 0.3, 0.9) ≈ slow_arc_length(c, 0.3, 0.9)
        @test DT.arc_length(c, 0.2, 0.21) ≈ slow_arc_length(c, 0.2, 0.21)
        for _ in 1:1000
            t₁, t₂ = rand(2)
            t₁, t₂ = minmax(t₁, t₂)
            @test DT.arc_length(c, t₁, t₂) ≈ slow_arc_length(c, t₁, t₂)
        end
    end

    ## Differentiate 
    for c in (arc, negarc, revarc, revnegarc, circ, revcirc)
        for t in LinRange(0, 1, 100)
            der1 = DT.differentiate(c, t)
            h = 1e-8
            der2 = (c(t + h) .- c(t - h)) ./ (2h)
            @test der1 ⪧ der2 rtol = 1e-5 atol = 1e-5
        end
    end

    ## Twice differentiate 
    for c in (arc, negarc, revarc, revnegarc, circ, revcirc)
        for t in LinRange(0, 1, 100)
            der1 = DT.twice_differentiate(c, t)
            h = 1e-4
            der2 = (c(t + h) .- 2 .* c(t) .+ c(t - h)) ./ (h^2)
            @test der1 ⪧ der2 rtol = 1e-5 atol = 1e-5
        end
    end

    ## Curvature 
    for c in (arc, negarc, revarc, revnegarc)
        for t in LinRange(0, 1, 100)
            ∂x, ∂y = DT.differentiate(c, t)
            ∂²x, ∂²y = DT.twice_differentiate(c, t)
            vec1 = [∂x, ∂y, 0.0]
            vec2 = [∂²x, ∂²y, 0.0]
            cur1 = DT.curvature(c, t)
            cur2 = cross(vec1, vec2)[3] / norm(vec1)^3
            @test cur1 ≈ cur2 rtol = 1e-5 atol = 1e-5
            @test cur1 ≈ sign(c.sector_angle) / c.radius
        end
    end

    ## Total variation
    for c in (arc, negarc, revarc, revnegarc)
        TV = DT.total_variation(c)
        @inferred DT.total_variation(c)
        TVslow = slow_total_absolute_curvature(c, 0, 1)
        @test TV ≈ TVslow ≈ abs(c.sector_angle)
        for _ in 1:1000
            t1, t2 = minmax(rand(2)...)
            TV = DT.total_variation(c, t1, t2)
            @inferred DT.total_variation(c, t1, t2)
            TVslow = slow_total_absolute_curvature(c, t1, t2)
            @test TV ≈ TVslow rtol = 1e-1 atol = 1e-1
        end
    end

    ## Equidistant split 
    for c in (arc, negarc, revarc, revnegarc)
        for _ in 1:1000
            t1, t2 = minmax(rand(2)...)
            @inferred DT.get_equidistant_split(c, t1, t2)
            t = DT.get_equidistant_split(c, t1, t2)
            s1 = DT.arc_length(c, t1, t)
            s2 = DT.arc_length(c, t, t2)
            @test s1 ≈ s2 rtol = 1e-1 atol = 1e-1
            @test t1 ≤ t ≤ t2
            @test DT.arc_length(c, t1, t2) ≈ 2s1
        end
    end

    ## Equivariation split 
    for c in (arc, negarc, revarc, revnegarc)
        for _ in 1:1000
            t1, t2 = minmax(rand(2)...)
            @inferred DT.get_equivariation_split(c, t1, t2)
            t, T = DT.get_equivariation_split(c, t1, t2)
            T1 = DT.total_variation(c, t1, t)
            T2 = DT.total_variation(c, t, t2)
            @test T1 ≈ T2 rtol = 1e-1 atol = 1e-1
            @test t1 ≤ t ≤ t2
            @test T ≈ T1
        end
    end

    ## Inverse 
    for c in (arc, negarc, revarc, revnegarc, circ, revcirc)
        for t in LinRange(0, 1, 2500)
            t′ = DT.get_inverse(c, c(t))
            @inferred DT.get_inverse(c, c(t))
            if c ∉ (circ, revcirc) || (t ≠ 0.0 && t ≠ 1.0)
                @test t ≈ t′
            else
                @test t′ ≈ 0.0 || t′ ≈ 1.0
            end
        end
    end

    ## get_circle_intersection
    c = CircularArc((0.0, 0.0), (1.0, 1.0), (0.5, 0.5))
    t₁, t₂, r = 0.0, 0.5, 0.5
    t′, q′ = DT.get_circle_intersection(c, t₁, t₂, r)
    @test t′ ≈ 0.22995991983967937 && q′ ⪧ (0.45551153816025436, -0.2057058712828833)
    @test DT.dist(c(t₁), q′) ≈ r rtol = 1e-3
    c = CircularArc((0.0, 0.0), (0.0, 0.0), (0.5, 0.5))
    t₁, t₂, r = 0.0, 0.5, 1.3
    t′, q′ = DT.get_circle_intersection(c, t₁, t₂, r)
    @test t′ ≈ 0.37124248496993983 && q′ ⪧ (1.2069097223005656, 0.48330735141035414)
    @test DT.dist(c(t₁), q′) ≈ r rtol = 1e-3

    ## == 
    @test arc ≠ negarc
    @test CircularArc((0.0, 0.0), (1.0, 1.0), (0.5, 0.5)) == CircularArc((0.0, 0.0), (1.0, 1.0), (0.5, 0.5))
    @test arc ≠ revarc
    @test arc ≠ revnegarc
    @test circ ≠ revcirc
    @test circ ≠ arc
end

@testset "EllipticalArc" begin
    ## Construction
    Rx, Ry, Cx, Cy, s = 2.2, 2.75, 0.6, 0.65, 60.0
    A, B = (-1.9827767129992, 1.4963893939715), (2.5063071261053, -1.2608915244961)
    arc = EllipticalArc(A, B, (Cx, Cy), Rx, Ry, s)
    @test DT.is_curve_bounded(arc)
    @test !DT.is_piecewise_linear(arc)
    @test !DT.is_interpolating(arc)
    @inferred EllipticalArc(A, B, (Cx, Cy), Rx, Ry, s)
    negarc = EllipticalArc(A, B, (Cx, Cy), Rx, Ry, s, positive=false)
    revarc = EllipticalArc(B, A, (Cx, Cy), Rx, Ry, s)
    revnegarc = EllipticalArc(B, A, (Cx, Cy), Rx, Ry, s, positive=false)
    @test arc.center == negarc.center == revarc.center == revnegarc.center == (Cx, Cy)
    @test arc.first == negarc.first == revarc.last == revnegarc.last == A
    @test arc.last == negarc.last == revarc.first == revnegarc.first == B
    @test arc.vert_radius == negarc.vert_radius == revarc.vert_radius == revnegarc.vert_radius == Ry
    @test arc.horz_radius == negarc.horz_radius == revarc.horz_radius == revnegarc.horz_radius == Rx
    @test arc.sector_angle ≈ 2.5603307829868798
    @test negarc.sector_angle ≈ arc.sector_angle - 2π
    @test revarc.sector_angle ≈ 2π - arc.sector_angle
    @test revnegarc.sector_angle ≈ -arc.sector_angle
    @test arc.rotation_scales ⪧ negarc.rotation_scales ⪧ revarc.rotation_scales ⪧ revnegarc.rotation_scales ⪧ sincosd(s)
    @test arc.sector_angle > 0
    @test negarc.sector_angle < 0
    @test revarc.sector_angle > 0
    @test revnegarc.sector_angle < 0
    closed_arc = EllipticalArc(A, A, (Cx, Cy), Rx, Ry, s)
    revclosed_arc = EllipticalArc(A, A, (Cx, Cy), Rx, Ry, s, positive=false)
    @test closed_arc.sector_angle ≈ 2π
    @test revclosed_arc.sector_angle ≈ -2π
    @test closed_arc.center == revclosed_arc.center == (Cx, Cy)
    @test closed_arc.first == revclosed_arc.first == closed_arc.last == revclosed_arc.last == A
    @test closed_arc.vert_radius == revclosed_arc.vert_radius == Ry
    @test closed_arc.horz_radius == revclosed_arc.horz_radius == Rx

    ## Evaluation
    @test arc(0.0) ⪧ negarc(0.0) ⪧ revarc(1.0) ⪧ revnegarc(1.0) ⪧ A
    @test arc(1.0) ⪧ negarc(1.0) ⪧ revarc(0.0) ⪧ revnegarc(0.0) ⪧ B
    @test arc(1e-16) ⪧ negarc(1e-16) ⪧ revarc(1 - 1e-16) ⪧ revnegarc(1 - 1e-16) ⪧ A
    @test arc(1 - 1e-16) ⪧ negarc(1 - 1e-16) ⪧ revarc(1e-16) ⪧ revnegarc(1e-16) ⪧ B
    @test closed_arc(0.0) ⪧ revclosed_arc(0.0) ⪧ A
    @test closed_arc(1.0) ⪧ revclosed_arc(1.0) ⪧ A
    @test closed_arc(1e-16) ⪧ revclosed_arc(1e-16) ⪧ A
    @test closed_arc(1 - 1e-16) ⪧ revclosed_arc(1 - 1e-16) ⪧ A
    arc1, arc2, arc3 = arc(0), arc(1 / 2), arc(1)
    negarc1, negarc2, negarc3 = negarc(0), negarc(1 / 2), negarc(1)
    revarc1, revarc2, revarc3 = revarc(0), revarc(1 / 2), revarc(1)
    revnegarc1, revnegarc2, revnegarc3 = revnegarc(0), revnegarc(1 / 2), revnegarc(1)
    closed_arc1, closed_arc2, closed_arc3 = closed_arc(0), closed_arc(1 / 2), closed_arc(0.95)
    revclosed_arc1, revclosed_arc2, revclosed_arc3 = revclosed_arc(0), revclosed_arc(1 / 2), revclosed_arc(0.95)
    @test DT.is_positively_oriented(DT.triangle_orientation(arc1, arc2, arc3))
    @test DT.is_negatively_oriented(DT.triangle_orientation(negarc1, negarc2, negarc3))
    @test DT.is_positively_oriented(DT.triangle_orientation(revarc1, revarc2, revarc3))
    @test DT.is_negatively_oriented(DT.triangle_orientation(revnegarc1, revnegarc2, revnegarc3))
    @test DT.is_positively_oriented(DT.triangle_orientation(closed_arc1, closed_arc2, closed_arc3))
    @test DT.is_negatively_oriented(DT.triangle_orientation(revclosed_arc1, revclosed_arc2, revclosed_arc3))
    t = LinRange(0, 1, 1000)
    for t in t
        arc_eval = arc(t)
        negarc_eval = negarc(t)
        revarc_eval = revarc(t)
        revnegarc_eval = revnegarc(t)
        closedarc_eval = closed_arc(t)
        revclosedarc_eval = revclosed_arc(t)
        arc_θ = arc.sector_angle * t + arc.start_angle
        negarc_θ = negarc.sector_angle * t + negarc.start_angle
        revarc_θ = revarc.sector_angle * t + revarc.start_angle
        revnegarc_θ = revnegarc.sector_angle * t + revnegarc.start_angle
        closedarc_θ = closed_arc.sector_angle * t + closed_arc.start_angle
        revclosedarc_θ = revclosed_arc.sector_angle * t + revclosed_arc.start_angle
        arc_eval′ = (Rx * cos(arc_θ) * cosd(s) - Ry * sin(arc_θ) * sind(s) + Cx, Rx * cos(arc_θ) * sind(s) + Ry * sin(arc_θ) * cosd(s) + Cy)
        negarc_eval′ = (Rx * cos(negarc_θ) * cosd(s) - Ry * sin(negarc_θ) * sind(s) + Cx, Rx * cos(negarc_θ) * sind(s) + Ry * sin(negarc_θ) * cosd(s) + Cy)
        revarc_eval′ = (Rx * cos(revarc_θ) * cosd(s) - Ry * sin(revarc_θ) * sind(s) + Cx, Rx * cos(revarc_θ) * sind(s) + Ry * sin(revarc_θ) * cosd(s) + Cy)
        revnegarc_eval′ = (Rx * cos(revnegarc_θ) * cosd(s) - Ry * sin(revnegarc_θ) * sind(s) + Cx, Rx * cos(revnegarc_θ) * sind(s) + Ry * sin(revnegarc_θ) * cosd(s) + Cy)
        closedarc_eval′ = (Rx * cos(closedarc_θ) * cosd(s) - Ry * sin(closedarc_θ) * sind(s) + Cx, Rx * cos(closedarc_θ) * sind(s) + Ry * sin(closedarc_θ) * cosd(s) + Cy)
        revclosedarc_eval′ = (Rx * cos(revclosedarc_θ) * cosd(s) - Ry * sin(revclosedarc_θ) * sind(s) + Cx, Rx * cos(revclosedarc_θ) * sind(s) + Ry * sin(revclosedarc_θ) * cosd(s) + Cy)
        @test arc_eval ⪧ arc_eval′
        @test negarc_eval ⪧ negarc_eval′
        @test revarc_eval ⪧ revarc_eval′
        @test revnegarc_eval ⪧ revnegarc_eval′
        @test closedarc_eval ⪧ closedarc_eval′
        @test revclosedarc_eval ⪧ revclosedarc_eval′
    end

    ## Sideof 
    d = (3.5, 3.0) # outside
    arc_cert = DT.point_position_relative_to_curve(arc, d)
    negarc_cert = DT.point_position_relative_to_curve(negarc, d)
    revarc_cert = DT.point_position_relative_to_curve(revarc, d)
    revnegarc_cert = DT.point_position_relative_to_curve(revnegarc, d)
    @test DT.is_right(arc_cert) && DT.is_left(negarc_cert) && DT.is_right(revarc_cert) && DT.is_left(revnegarc_cert)
    d = (0.0, 2.0) # inside
    arc_cert = DT.point_position_relative_to_curve(arc, d)
    negarc_cert = DT.point_position_relative_to_curve(negarc, d)
    revarc_cert = DT.point_position_relative_to_curve(revarc, d)
    revnegarc_cert = DT.point_position_relative_to_curve(revnegarc, d)
    @test DT.is_left(arc_cert) && DT.is_right(negarc_cert) && DT.is_left(revarc_cert) && DT.is_right(revnegarc_cert)
    d = (-1.0373942027943, -0.869928168413) # inside 
    arc_cert = DT.point_position_relative_to_curve(arc, d)
    negarc_cert = DT.point_position_relative_to_curve(negarc, d)
    revarc_cert = DT.point_position_relative_to_curve(revarc, d)
    revnegarc_cert = DT.point_position_relative_to_curve(revnegarc, d)
    @test DT.is_left(arc_cert) && DT.is_right(negarc_cert) && DT.is_left(revarc_cert) && DT.is_right(revnegarc_cert)
    d = (-1.0385707473021, -0.8709576448573) # outside 
    arc_cert = DT.point_position_relative_to_curve(arc, d)
    negarc_cert = DT.point_position_relative_to_curve(negarc, d)
    revarc_cert = DT.point_position_relative_to_curve(revarc, d)
    revnegarc_cert = DT.point_position_relative_to_curve(revnegarc, d)
    @test DT.is_right(arc_cert) && DT.is_left(negarc_cert) && DT.is_right(revarc_cert) && DT.is_left(revnegarc_cert)

    t = LinRange(0, 1, 15000) |> collect
    pop!(t)
    points = closed_arc.(t)
    boundary_nodes = [eachindex(t); 1]
    reverse_points = revclosed_arc.(t)
    reverse_boundary_nodes = [eachindex(t); 1]
    for _ in 1:10000
        p = (2rand() - 1, 2rand() - 1)
        δ = DT.distance_to_polygon(p, points, boundary_nodes)
        δ2 = DT.distance_to_polygon(p, reverse_points, reverse_boundary_nodes)
        cert = DT.point_position_relative_to_curve(closed_arc, p)
        cert2 = DT.point_position_relative_to_curve(revclosed_arc, p)
        if δ > 0
            @test DT.is_left(cert)
        elseif δ < 0
            @test DT.is_right(cert)
        else
            @test DT.is_on(cert)
        end
        if δ2 < 0
            @test DT.is_left(cert2)
        elseif δ2 > 0
            @test DT.is_right(cert2)
        else
            @test DT.is_on(cert2)
        end
    end

    ## Differentiate 
    for c in (arc, negarc, revarc, revnegarc, closed_arc, revclosed_arc)
        for t in LinRange(0, 1, 100)
            der1 = DT.differentiate(c, t)
            h = 1e-8
            der2 = (c(t + h) .- c(t - h)) ./ (2h)
            @test der1 ⪧ der2 rtol = 1e-5 atol = 1e-5
        end
    end

    ## Arclength 
    for c in (arc, negarc, revarc, revnegarc, closed_arc, revclosed_arc)
        @test DT.arc_length(c) ≈ slow_arc_length(c, 0, 1)
        @inferred DT.arc_length(c)
        @test DT.arc_length(c, 0.1, 0.15) ≈ slow_arc_length(c, 0.1, 0.15)
        @test DT.arc_length(c, 0.0, 0.0) ≈ 0.0
        @test DT.arc_length(c, 1.0, 1.0) ≈ 0.0
        @test DT.arc_length(c, 0.0, 1.0) ≈ slow_arc_length(c, 0.0, 1.0)
        @test DT.arc_length(c, 0.3, 0.9) ≈ slow_arc_length(c, 0.3, 0.9)
        @test DT.arc_length(c, 0.2, 0.21) ≈ slow_arc_length(c, 0.2, 0.21)
        @inferred DT.arc_length(c, 0.2, 0.21)
        for _ in 1:1000
            t₁, t₂ = rand(2)
            t₁, t₂ = minmax(t₁, t₂)
            @test DT.arc_length(c, t₁, t₂) ≈ slow_arc_length(c, t₁, t₂)
        end
    end

    ## Twice differentiate 
    for c in (arc, negarc, revarc, revnegarc, closed_arc, revclosed_arc)
        for t in LinRange(0, 1, 100)
            der1 = DT.twice_differentiate(c, t)
            @inferred DT.twice_differentiate(c, t)
            h = 1e-4
            der2 = (c(t + h) .- 2 .* c(t) .+ c(t - h)) ./ (h^2)
            @test der1 ⪧ der2 rtol = 1e-5 atol = 1e-5
        end
    end

    ## Curvature 
    for c in (arc, negarc, revarc, revnegarc, closed_arc, revclosed_arc)
        for t in LinRange(0, 1, 100)
            ∂x, ∂y = DT.differentiate(c, t)
            ∂²x, ∂²y = DT.twice_differentiate(c, t)
            vec1 = [∂x, ∂y, 0.0]
            vec2 = [∂²x, ∂²y, 0.0]
            cur1 = DT.curvature(c, t)
            cur2 = cross(vec1, vec2)[3] / norm(vec1)^3
            @test cur1 ≈ cur2 rtol = 1e-5 atol = 1e-5
        end
    end

    ## Total variation
    for c in (arc, negarc, revarc, revnegarc, closed_arc, revclosed_arc)
        TV = DT.total_variation(c)
        @inferred DT.total_variation(c)
        TVslow = slow_total_absolute_curvature(c, 0, 1)
        @test TV ≈ TVslow
        for _ in 1:1000
            t1, t2 = minmax(rand(2)...)
            TV = DT.total_variation(c, t1, t2)
            @inferred DT.total_variation(c, t1, t2)
            TVslow = slow_total_absolute_curvature(c, t1, t2)
            @test TV ≈ TVslow rtol = 1e-1 atol = 1e-1
        end
    end

    ## Equidistant split
    for c in (arc, negarc, revarc, revnegarc, closed_arc, revclosed_arc)
        for _ in 1:1000
            t1, t2 = minmax(rand(2)...)
            @inferred DT.get_equidistant_split(c, t1, t2)
            t = DT.get_equidistant_split(c, t1, t2)
            s1 = DT.arc_length(c, t1, t)
            s2 = DT.arc_length(c, t, t2)
            @test s1 ≈ s2 rtol = 1e-1 atol = 1e-1
            @test t1 ≤ t ≤ t2
            @test DT.arc_length(c, t1, t2) ≈ 2s1 rtol = 1e-1
        end
    end

    ## Equivariation split 
    for c in (arc, negarc, revarc, revnegarc, closed_arc, revclosed_arc)
        for _ in 1:1000
            t1, t2 = minmax(rand(2)...)
            @inferred DT.get_equivariation_split(c, t1, t2)
            t, T = DT.get_equivariation_split(c, t1, t2)
            T1 = DT.total_variation(c, t1, t)
            T2 = DT.total_variation(c, t, t2)
            @test T1 ≈ T2 rtol = 1e-1 atol = 1e-1
            @test t1 ≤ t ≤ t2
            @test T ≈ T1 rtol = 1e-1 atol = 1e-1
        end
    end

    ## Inverse 
    for c in (arc, negarc, revarc, revnegarc, closed_arc, revclosed_arc)
        for t in LinRange(0, 1, 2500)
            t′ = DT.get_inverse(c, c(t))
            @inferred DT.get_inverse(c, c(t))
            if c ∉ (closed_arc, revclosed_arc) || (t ≠ 0.0 && t ≠ 1.0)
                @test t ≈ t′
            else
                @test t′ ≈ 0.0 || t′ ≈ 1.0
            end
        end
    end

    ## get_circle_intersection
    Rx, Ry, Cx, Cy, s = 2.2, 2.75, 0.6, 0.65, 60.0
    A, B = (-1.9827767129992, 1.4963893939715), (2.5063071261053, -1.2608915244961)
    c = EllipticalArc(A, B, (Cx, Cy), Rx, Ry, s)
    t₁, t₂, r, = 0.2, 0.5, 0.8
    t′, q′ = DT.get_circle_intersection(c, t₁, t₂, r)
    @test t′ ≈ 0.32234468937857575 && q′ ⪧ (-1.489595166177109, -0.3863541400710244)
    @test DT.dist(c(t₁), q′) ≈ r rtol = 1e-3
    t₁, t₂, r, = 1.0, 0.5, 0.18889
    t′, q′ = DT.get_circle_intersection(c, t₁, t₂, r)
    @test t′ ≈ 0.967434869739479 && q′ ⪧ (2.349597231970133, -1.368111418490932)
    @test DT.dist(c(t₁), q′) ≈ r rtol = 1e-2

    ## == 
    @test arc ≠ negarc
    @test EllipticalArc(A, B, (Cx, Cy), Rx, Ry, s) == EllipticalArc(A, B, (Cx, Cy), Rx, Ry, s)
    @test arc ≠ revarc
    @test arc ≠ revnegarc
    @test closed_arc ≠ revclosed_arc
    @test closed_arc ≠ arc
end

@testset "BezierCurve" begin
    ## Construction
    control_points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (1 / 2, 1 / 2)]
    bezier = BezierCurve(control_points)
    @test DT.is_curve_bounded(bezier)
    @test !DT.is_piecewise_linear(bezier)
    @test !DT.is_interpolating(bezier)
    @test bezier.control_points == control_points
    @test typeof(bezier.cache) == Vector{NTuple{2,Float64}} && length(bezier.cache) == 5
    @test length(bezier.lookup_table) == 5000
    for i in 1:5000
        t = (i - 1) / 4999
        @test bezier(t) ⪧ bezier.lookup_table[i]
    end

    ## Evaluation
    @test bezier(0.0) ⪧ (0.0, 0.0)
    @test bezier(1e-16) ⪧ (0.0, 0.0) atol = 1e-9
    @test bezier(1.0) ⪧ bezier(1 - 1e-16) ⪧ (1 / 2, 1 / 2)
    t = LinRange(0, 1, 1500)
    for t in t
        @test bezier(t) ⪧ slow_bezier_eval(control_points, t)
        @inferred bezier(t)
    end

    ## Differentiation  
    t = LinRange(0, 1, 1500)
    for t in t
        der1 = ForwardDiff.derivative(t -> slow_bezier_eval(control_points, t), t)
        der2 = (bezier(t + 1e-6) .- bezier(t - 1e-6)) ./ (2 * 1e-6)
        der3 = DT.differentiate(bezier, t)
        @test der1 ⪧ der2 ⪧ der3
        @inferred DT.differentiate(bezier, t)
    end

    ## Closest point 
    invert(points) =
        let y = maximum(last.(points))
            return [(x, y - y′) for (x, y′) in points]
        end # just to match my reference figure
    control_points = invert([(207.0, 65.0), (84.0, 97.0), (58.0, 196.0), (261.0, 130.0), (216.0, 217.0), (54.0, 122.0), (52.0, 276.0), (85.0, 330.0), (258.0, 334.0), (209.0, 286.0)])
    for lookup_steps in (100, 500, 1000)
        if lookup_steps ≠ 500
            bezier = BezierCurve(control_points; lookup_steps)
        else
            bezier = BezierCurve(control_points)
        end
        p = (100.0, 100.0)
        t′, q = DT.get_closest_point(bezier, p)
        @test q ⪧ bezier(t′)
        for _ in 1:50
            p = rand(2) * 300 |> Tuple
            t′, q = DT.get_closest_point(bezier, p)
            @test q ⪧ bezier(t′)
            _t, _q = closest_point_on_curve(bezier, p)
            @test _t ≈ t′ rtol = 1e-1 atol = 1e-1
            @test q ⪧ _q rtol = 1e-1 atol = 1e-1
            @test DT.dist(p, _q) ≈ DT.dist(p, q) rtol = 1e-1 atol = 1e-1
        end
        for t in LinRange(0, 1, 150)
            p = bezier(t)
            t′, q = DT.get_closest_point(bezier, p)
            @test q ⪧ bezier(t′)
            @test t ≈ t′ rtol = 1e-1 atol = 1e-1
            @test p ⪧ q rtol = 1e-1 atol = 1e-1
        end
        @test DT.get_closest_point(bezier, bezier(0.0))[1] == 0.0
        @test DT.get_closest_point(bezier, bezier(1.0))[1] == 1.0
    end

    ## Sideof 
    control_points = invert([(207.0, 65.0), (84.0, 97.0), (58.0, 196.0), (261.0, 130.0), (216.0, 217.0), (54.0, 122.0), (52.0, 276.0), (85.0, 330.0), (258.0, 334.0), (209.0, 286.0)])
    bezier = BezierCurve(control_points)
    t = LinRange(0, 1, 15000)
    d1 = (150.0, 200.0)
    d2 = (125.0, 200.0)
    d3 = (115.0, 100.0)
    d4 = (112.0, 100.0)
    d5 = (200.0, 40.0)
    d6 = (225.0, 50.0)
    d7 = (200.0, 250.0)
    d8 = (250.0, 270.0)
    d9 = (207.0, 269.0)
    d10 = (209.0, 48.0)
    cert1 = DT.point_position_relative_to_curve(bezier, d1) # left 
    cert2 = DT.point_position_relative_to_curve(bezier, d2) # right
    cert3 = DT.point_position_relative_to_curve(bezier, d3) # left
    cert4 = DT.point_position_relative_to_curve(bezier, d4) # right
    cert5 = DT.point_position_relative_to_curve(bezier, d5) # left
    cert6 = DT.point_position_relative_to_curve(bezier, d6) # right
    cert7 = DT.point_position_relative_to_curve(bezier, d7) # left
    cert8 = DT.point_position_relative_to_curve(bezier, d8) # left
    cert9 = DT.point_position_relative_to_curve(bezier, d9) # on
    cert10 = DT.point_position_relative_to_curve(bezier, d10) # on
    @test DT.is_left(cert1) && DT.is_right(cert2) && DT.is_left(cert3) && DT.is_right(cert4) && DT.is_left(cert5) && DT.is_right(cert6) && DT.is_left(cert7) && DT.is_left(cert8) && DT.is_on(cert9) && DT.is_on(cert10)

    control_points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]
    bezier = BezierCurve(control_points)
    t = LinRange(0, 1, 15000)
    d1 = (0.25, 0.1)
    d2 = (0.01, 0.01)
    d3 = (0.5, 0.5)
    d4 = (0.25, 0.0)
    d5 = (0.0, 0.0)
    d6 = (0.75, 0.5)
    d7 = (1.0, 1.0)
    cert1 = DT.point_position_relative_to_curve(bezier, d1) # left
    cert2 = DT.point_position_relative_to_curve(bezier, d2) # left 
    cert3 = DT.point_position_relative_to_curve(bezier, d3) # Left
    cert4 = DT.point_position_relative_to_curve(bezier, d4) # right
    cert5 = DT.point_position_relative_to_curve(bezier, d5) # on
    cert6 = DT.point_position_relative_to_curve(bezier, d6) # right
    cert7 = DT.point_position_relative_to_curve(bezier, d7) # right
    @test DT.is_left(cert1) && DT.is_left(cert2) && DT.is_left(cert3) && DT.is_right(cert4) && DT.is_on(cert5) && DT.is_right(cert6) && DT.is_right(cert7)

    control_points = [(0.0, 0.0), (1.0, 1.0)]
    bezier = BezierCurve(control_points)
    for _ in 1:100
        p = randn(2) |> Tuple
        cert = DT.point_position_relative_to_curve(bezier, p)
        cert2 = DT.point_position_relative_to_curve(LineSegment((0.0, 0.0), (1.0, 1.0)), p)
        @test cert == cert2
    end

    control_points = [(0.0, 0.0), (1 / 2, 1 / 2), (0.0, 1.0)]
    bezier = BezierCurve(control_points)
    d1 = (0.1, 0.5)
    d2 = (0.2, 0.1)
    d3 = (0.2, 0.5)
    d4 = (0.2, 1.0)
    d5 = (0.0, 0.0)
    d6 = (0.0, 1.0)
    cert1 = DT.point_position_relative_to_curve(bezier, d1) # left
    cert2 = DT.point_position_relative_to_curve(bezier, d2) # right
    cert3 = DT.point_position_relative_to_curve(bezier, d3) # left 
    cert4 = DT.point_position_relative_to_curve(bezier, d4) # right
    cert5 = DT.point_position_relative_to_curve(bezier, d5) # on
    cert6 = DT.point_position_relative_to_curve(bezier, d6) # on
    @test DT.is_left(cert1) && DT.is_right(cert2) && DT.is_left(cert3) && DT.is_right(cert4) && DT.is_on(cert5) && DT.is_on(cert6)

    ## Arclength 
    control_points1 = invert([(207.0, 65.0), (84.0, 97.0), (58.0, 196.0), (261.0, 130.0), (216.0, 217.0), (54.0, 122.0), (52.0, 276.0), (85.0, 330.0), (258.0, 334.0), (209.0, 286.0)])
    control_points2 = [(0.0, 0.0), (1.0, 1.0)]
    control_points3 = [(0.0, 0.0), (1 / 2, 1 / 2), (0.0, 1.0)]
    control_points4 = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]
    bezier1 = BezierCurve(control_points1)
    bezier2 = BezierCurve(control_points2)
    bezier3 = BezierCurve(control_points3)
    bezier4 = BezierCurve(control_points4)
    for bezier in (bezier1, bezier2, bezier3, bezier4)
        @test DT.arc_length(bezier) ≈ slow_arc_length(bezier, 0, 1) rtol = 1e-4
        @inferred DT.arc_length(bezier)
        @test DT.arc_length(bezier, 0.1, 0.15) ≈ slow_arc_length(bezier, 0.1, 0.15)
        @test DT.arc_length(bezier, 0.0, 0.0) ≈ 0.0
        @test DT.arc_length(bezier, 1.0, 1.0) ≈ 0.0
        @test DT.arc_length(bezier, 0.0, 1.0) ≈ slow_arc_length(bezier, 0.0, 1.0) rtol = 1e-4
        @test DT.arc_length(bezier, 0.3, 0.9) ≈ slow_arc_length(bezier, 0.3, 0.9)
        @test DT.arc_length(bezier, 0.2, 0.21) ≈ slow_arc_length(bezier, 0.2, 0.21)
        for _ in 1:1000
            t₁, t₂ = rand(2)
            t₁, t₂ = minmax(t₁, t₂)
            @test DT.arc_length(bezier, t₁, t₂) ≈ slow_arc_length(bezier, t₁, t₂) rtol = 1e-2
        end
    end

    ## Twice differentiate 
    for bezier in (bezier1, bezier2, bezier3, bezier4)
        for t in LinRange(0, 1, 1000)
            der1 = DT.twice_differentiate(bezier, t)
            h = 1e-4
            der2 = (bezier(t + h) .- 2 .* bezier(t) .+ bezier(t - h)) ./ (h^2)
            @test der1 ⪧ der2 rtol = 1e-5 atol = 1e-5
        end
    end

    ## Thrice differentiate
    for bezier in (bezier1, bezier2, bezier3, bezier4)
        for t in LinRange(0, 1, 1000)
            der1 = DT.thrice_differentiate(bezier, t)
            h = 1e-6
            f = t -> DT.twice_differentiate(bezier, t)
            der2 = (f(t + h) .- f(t - h)) ./ (2h)
            @test der1 ⪧ der2 rtol = 1e-3 atol = 1e-3
        end
    end

    ## Curvature 
    for c in (bezier1, bezier2, bezier3, bezier4)
        for t in LinRange(0, 1, 100)
            ∂x, ∂y = DT.differentiate(c, t)
            ∂²x, ∂²y = DT.twice_differentiate(c, t)
            vec1 = [∂x, ∂y, 0.0]
            vec2 = [∂²x, ∂²y, 0.0]
            cur1 = DT.curvature(c, t)
            cur2 = cross(vec1, vec2)[3] / norm(vec1)^3
            @test cur1 ≈ cur2 rtol = 1e-5 atol = 1e-5
        end
    end

    ## Orientation markers
    control_points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (1 / 2, 1 / 2)]
    bezier = BezierCurve(control_points)
    tx = DT.horizontal_turning_points(bezier)
    ty = DT.vertical_turning_points(bezier)
    @test tx ≈ [0.3785321505748616, 0.8559779432249344]
    @test ty ≈ [0.0, 0.7101020514433644]
    κx = DT.horizontal_inflection_points(bezier)
    κy = DT.vertical_inflection_points(bezier)
    @test κx ≈ [0.6324555320336759]
    @test κy ≈ [0.3101020514433644]
    κ = DT.inflection_points(bezier)
    @test isempty(κ)
    allt = DT.orientation_markers(bezier)
    @test allt ≈ sort([tx; ty; κx; κy; κ; 1.0])
    @test allt == bezier.orientation_markers

    control_points = invert([(207.0, 65.0), (84.0, 97.0), (58.0, 196.0), (261.0, 130.0), (216.0, 217.0), (54.0, 122.0), (52.0, 276.0), (85.0, 330.0), (258.0, 334.0), (209.0, 286.0)])
    bezier = BezierCurve(control_points)
    tx = DT.horizontal_turning_points(bezier)
    ty = DT.vertical_turning_points(bezier)
    @test tx ≈ [0.15407958022481555, 0.36232751887181364, 0.6486424410086364, 0.9665387508202394]
    @test ty ≈ [0.8839329480604027]
    κx = DT.horizontal_inflection_points(bezier)
    κy = DT.vertical_inflection_points(bezier)
    @test κx ≈ [0.24362760326392038, 0.511382356213115, 0.8499089636970797]
    @test κy ≈ [0.06532566474143786, 0.3462184144341068, 0.6727307742621746]
    κ = DT.inflection_points(bezier)
    @test κ ≈ [0.2645455199588649, 0.47387783454147714]
    allt = DT.orientation_markers(bezier)
    @test allt ≈ sort([tx; ty; κx; κy; κ; 0.0; 1.0])
    @test allt == bezier.orientation_markers

    control_points = [(0.0, 0.0), (1.0, 1.0)]
    bezier = BezierCurve(control_points)
    tx = DT.horizontal_turning_points(bezier)
    ty = DT.vertical_turning_points(bezier)
    @test isempty(tx) && isempty(ty)
    κx = DT.horizontal_inflection_points(bezier)
    κy = DT.vertical_inflection_points(bezier)
    @test isempty(κx) && isempty(κy)
    κ = DT.inflection_points(bezier)
    @test isempty(κ)
    allt = DT.orientation_markers(bezier)
    @test allt ≈ [0.0, 1.0]
    @test allt == bezier.orientation_markers

    control_points = [(0.0, 0.0), (1 / 2, 1 / 2), (0.0, 1.0)]
    bezier = BezierCurve(control_points)
    tx = DT.horizontal_turning_points(bezier)
    ty = DT.vertical_turning_points(bezier)
    @test tx ≈ [0.5]
    @test isempty(ty)
    κx = DT.horizontal_inflection_points(bezier)
    κy = DT.vertical_inflection_points(bezier)
    @test isempty(κx) && isempty(κy)
    κ = DT.inflection_points(bezier)
    @test isempty(κ)
    allt = DT.orientation_markers(bezier)
    @test allt ≈ [0.0, 0.5, 1.0]
    @test allt == bezier.orientation_markers

    control_points = [(-12.0, 5.0), (-9.77, 6.29),
        (-8.11, 4.55), (-7.47, 1.49),
        (-7.61, -1.29), (-9.63, -3.69),
        (-13.65, -4.37), (-15.65, -1.25),
        (-15.39, 0.93), (-14.17, 1.63),
        (-12.37, -0.93), (-13.51, -1.17),
        (-12.59, -2.39), (-10.6, -2.47),
        (-9.19, 0.11), (-9.95, 2.79)]
    bezier = BezierCurve(control_points)
    tx = DT.horizontal_turning_points(bezier)
    ty = DT.vertical_turning_points(bezier)
    @test tx ≈ [0.18075623805341032, 0.5604462043127618, 0.970975589403759]
    @test ty ≈ [0.034765208488561516, 0.3911416826751549, 0.6084649566889487, 0.8025198442401662]
    κx = DT.horizontal_inflection_points(bezier)
    κy = DT.vertical_inflection_points(bezier)
    @test κx ≈ [0.3486922619973087, 0.8601849568818399]
    @test κy ≈ [0.17686163858006337, 0.4891637919747223, 0.7138800043385305]
    κ = DT.inflection_points(bezier)
    @test κ ≈ [0.696315552183285]
    allt = DT.orientation_markers(bezier)
    @test allt ≈ sort([tx; ty; κx; κy; κ; 0.0; 1.0])
    @test allt == bezier.orientation_markers

    control_points = [(-12.0, -4.0), (-8.0, -8.0), (-4.0, -6.0), (-2.0, -4.0),
        (-6.0, 0.0), (-10.0, 4.0), (-14.0, 8.0), (-10.0, 12.0), (-4.0, 14.0), (4.0, 10.0),
        (0.0, 8.0), (-2.0, 6.0), (2.0, 4.0), (6.0, -2.0), (6.0, -8.0), (0.0, -12.0), (-10.0, -12.0),
        (-18.0, -12.0), (-18.0, -2.0), (-18.0, -2.0), (-12.0, -4.0)]
    bezier = BezierCurve(control_points)
    tx = DT.horizontal_turning_points(bezier)
    ty = DT.vertical_turning_points(bezier)
    @test tx ≈ [0.13161961540843561, 0.2704114277681997, 0.6092229033971368, 0.9251598879841408]
    @test ty ≈ [0.0525819800261663, 0.40273720300019505, 0.8007966486184426, 0.9669119698489039]
    κx = DT.horizontal_inflection_points(bezier)
    κy = DT.vertical_inflection_points(bezier)
    @test κx ≈ [0.0, 0.19546372693114597, 0.4126806384017861, 0.8024455147185351]
    @test κy ≈ [0.20849602194605005, 0.634052074772975, 0.8968171648064026]
    κ = DT.inflection_points(bezier)
    @test κ ≈ [0.19306516592293574, 1.0]
    allt = DT.orientation_markers(bezier)
    @test allt ≈ sort([tx; ty; κx; κy; κ])
    @test allt == bezier.orientation_markers

    ## Total variation
    control_points1 = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (1 / 2, 1 / 2)]
    control_points2 = invert([(207.0, 65.0), (84.0, 97.0), (58.0, 196.0), (261.0, 130.0), (216.0, 217.0), (54.0, 122.0), (52.0, 276.0), (85.0, 330.0), (258.0, 334.0), (209.0, 286.0)])
    control_points3 = [(0.0, 0.0), (1.0, 1.0)]
    control_points4 = [(0.0, 0.0), (1 / 2, 1 / 2), (0.0, 1.0)]
    control_points5 = [(-12.0, 5.0), (-9.77, 6.29),
        (-8.11, 4.55), (-7.47, 1.49),
        (-7.61, -1.29), (-9.63, -3.69),
        (-13.65, -4.37), (-15.65, -1.25),
        (-15.39, 0.93), (-14.17, 1.63),
        (-12.37, -0.93), (-13.51, -1.17),
        (-12.59, -2.39), (-10.6, -2.47),
        (-9.19, 0.11), (-9.95, 2.79)]
    control_points6 = [(-12.0, -4.0), (-8.0, -8.0), (-4.0, -6.0), (-2.0, -4.0),
        (-6.0, 0.0), (-10.0, 4.0), (-14.0, 8.0), (-10.0, 12.0), (-4.0, 14.0), (4.0, 10.0),
        (0.0, 8.0), (-2.0, 6.0), (2.0, 4.0), (6.0, -2.0), (6.0, -8.0), (0.0, -12.0), (-10.0, -12.0),
        (-18.0, -12.0), (-18.0, -2.0), (-18.0, -2.0), (-12.0, -4.0)] # periodic 
    control_points7 = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)] # periodic
    bezier1 = BezierCurve(control_points1)
    bezier2 = BezierCurve(control_points2)
    bezier3 = BezierCurve(control_points3)
    bezier4 = BezierCurve(control_points4)
    bezier5 = BezierCurve(control_points5)
    bezier6 = BezierCurve(control_points6)
    bezier7 = BezierCurve(control_points7)

    for c in (bezier1, bezier2, bezier3, bezier4, bezier5, bezier6, bezier7)
        TV = DT.total_variation(c)
        @inferred DT.total_variation(c)
        TVslow = slow_total_absolute_curvature(c, 0, 1)
        @test TV ≈ TVslow rtol = 1e-3 atol = 1e-3
        for _ in 1:1000
            t1, t2 = minmax(rand(2)...)
            TV = DT.total_variation(c, t1, t2)
            @inferred DT.total_variation(c, t1, t2)
            TVslow = slow_total_absolute_curvature(c, t1, t2)
            @test TV ≈ TVslow rtol = 1e-1 atol = 1e-1
        end
    end

    ## Equidistant split
    for c in (bezier1, bezier2, bezier3, bezier4, bezier5, bezier6, bezier7)
        for _ in 1:1000
            t1, t2 = minmax(rand(2)...)
            @inferred DT.get_equidistant_split(c, t1, t2)
            t = DT.get_equidistant_split(c, t1, t2)
            s1 = DT.arc_length(c, t1, t)
            s2 = DT.arc_length(c, t, t2)
            @test s1 ≈ s2 rtol = 1e-2 atol = 1e-2
            @test t1 ≤ t ≤ t2
            @test DT.arc_length(c, t1, t2) ≈ 2s1 rtol = 1e-1
        end
    end

    ## Equivariation split 
    for c in (bezier1, bezier2, bezier3, bezier4, bezier5, bezier6, bezier7)
        for _ in 1:1000
            t1, t2 = minmax(rand(2)...)
            @inferred DT.get_equivariation_split(c, t1, t2)
            t, T = DT.get_equivariation_split(c, t1, t2)
            T1 = DT.total_variation(c, t1, t)
            T2 = DT.total_variation(c, t, t2)
            @test T1 ≈ T2 rtol = 1e-1 atol = 1e-1
            @test t1 ≤ t ≤ t2
            @test T ≈ T1
        end
    end

    ## Inverse 
    for c in (bezier1, bezier2, bezier3, bezier4, bezier5, bezier6, bezier7)
        for t in LinRange(0, 1, 1000)
            t′ = DT.get_inverse(c, c(t))
            @inferred DT.get_inverse(c, c(t))
            if c ∉ (bezier7, bezier6) || (t ≠ 0.0 && t ≠ 1.0)
                @test t ≈ t′
            else
                @test t′ ≈ 0.0 || t′ ≈ 1.0
            end
        end
    end

    ## get_circle_intersection
    control_points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]
    c = BezierCurve(control_points)
    t₁, t₂, r, = 0.2, 0.5, 0.2
    t′, q′ = DT.get_circle_intersection(c, t₁, t₂, r)
    @test t′ ≈ 0.3019603920784157 && q′ ⪧ (0.6773887113970392, 0.34344590594423263)
    @test DT.dist(c(t₁), q′) ≈ r rtol = 1e-3
    t₁, t₂, r, = 1.0, 0.5, 0.6
    t′, q′ = DT.get_circle_intersection(c, t₁, t₂, r)
    @test t′ ≈ 0.7954590918183637 && q′ ⪧ (0.186063568848561, 0.5706424003210613)
    @test DT.dist(c(t₁), q′) ≈ r rtol = 1e-2

    ## == 
    ctrl1 = [(0.0, 0.0), (1.3, 1.5), (2.7, 5.3), (17.3, 5.2), (23.5, -0.5)]
    ctrl2 = control_points
    @test BezierCurve(ctrl1) == BezierCurve(ctrl1)
    @test BezierCurve(ctrl1) ≠ BezierCurve(ctrl2)
end

@testset "BSpline" begin
    ## Construction 
    control_points = [(1.3, 5.0), (10.0, 2.0), (17.3, 25.0), (17.3, 17.3)]
    spline = BSpline(control_points)
    @test DT.is_curve_bounded(spline)
    @test !DT.is_piecewise_linear(spline)
    @test !DT.is_interpolating(spline)
    @test spline.control_points ⪧ control_points
    @test typeof(spline.cache) == Vector{NTuple{2,Float64}} && length(spline.cache) == 4
    @test spline.knots == [0, 0, 0, 0, 1, 1, 1, 1]
    @test length(spline.lookup_table) == 5000
    for i in 1:5000
        @test spline.lookup_table[i] ⪧ spline((i - 1) / 4999)
    end

    periodic_control_points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (1 / 2, 1 / 2), (0.0, 0.0)]
    periodic_spline = BSpline(periodic_control_points)
    @test periodic_spline.control_points ⪧ [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (1 / 2, 1 / 2), (0.0, 0.0)]
    @test typeof(periodic_spline.cache) == Vector{NTuple{2,Float64}} && length(periodic_spline.cache) == 6
    @test periodic_spline.knots == [0, 0, 0, 0, 1, 2, 3, 3, 3, 3]
    @test length(periodic_spline.lookup_table) == 5000
    for i in 1:5000
        @test periodic_spline.lookup_table[i] ⪧ periodic_spline((i - 1) / 4999)
    end

    longer_control_points = [(0.3, 0.3), (0.5, -1.0), (2.0, 0.0), (2.5, 3.2), (-10.0, 10.0)]
    longer_spline = BSpline(longer_control_points; lookup_steps=2500)
    @test longer_spline.control_points ⪧ longer_control_points
    @test typeof(longer_spline.cache) == Vector{NTuple{2,Float64}} && length(longer_spline.cache) == 5
    @test longer_spline.knots == [0, 0, 0, 0, 1, 2, 2, 2, 2]
    @test length(longer_spline.lookup_table) == 2500
    for i in 1:2500
        @test longer_spline.lookup_table[i] ⪧ longer_spline((i - 1) / 2499)
    end

    quadratic_control_points = [(0.1, 0.1), (0.2, 0.2), (0.5, 0.8), (1.0, 2.0), (0.0, 15.0), (-5.0, 10.0), (-10.0, 0.0)]
    quadratic_spline = BSpline(quadratic_control_points; degree=2)
    @test quadratic_spline.control_points ⪧ quadratic_control_points
    @test typeof(quadratic_spline.cache) == Vector{NTuple{2,Float64}} && length(quadratic_spline.cache) == 7
    @test quadratic_spline.knots == [0, 0, 0, 1, 2, 3, 4, 5, 5, 5]
    @test length(quadratic_spline.lookup_table) == 5000
    for i in 1:5000
        @test quadratic_spline.lookup_table[i] ⪧ quadratic_spline((i - 1) / 4999)
    end

    sextic_control_points = [(cos(t), sin(t)) for t in LinRange(0, π - π / 3, 10)]
    sextic_spline = BSpline(sextic_control_points; degree=6)
    @test sextic_spline.control_points ⪧ sextic_control_points
    @test typeof(sextic_spline.cache) == Vector{NTuple{2,Float64}} && length(sextic_spline.cache) == 10
    @test sextic_spline.knots == [0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4]
    @test length(sextic_spline.lookup_table) == 5000
    for i in 1:5000
        @test sextic_spline.lookup_table[i] ⪧ sextic_spline((i - 1) / 4999)
    end

    quadratic_periodic_control_points = [(cos(t), sin(t)) for t in LinRange(0, 2π, 5)]
    quadratic_periodic_control_points[end] = quadratic_periodic_control_points[1]
    quadratic_periodic_spline = BSpline(quadratic_periodic_control_points; degree=2)
    @test quadratic_periodic_spline.control_points ⪧ quadratic_periodic_control_points
    @test typeof(quadratic_periodic_spline.cache) == Vector{NTuple{2,Float64}} && length(quadratic_periodic_spline.cache) == 5
    @test quadratic_periodic_spline.knots == [0, 0, 0, 1, 2, 3, 3, 3]
    @test length(quadratic_periodic_spline.lookup_table) == 5000
    for i in 1:5000
        @test quadratic_periodic_spline.lookup_table[i] ⪧ quadratic_periodic_spline((i - 1) / 4999)
    end
    @test length(quadratic_periodic_control_points) == 5 # test that we don't mutate

    ## Evaluation 
    for spl in (spline, periodic_spline, longer_spline, quadratic_spline, sextic_spline, quadratic_periodic_spline)
        local control_points, knots
        control_points = spl.control_points
        knots = spl.knots
        for t in LinRange(0, 1, 15000)
            @test slow_eval_bspline(control_points, knots, t) ⪧ spl(t)
        end
    end

    @test spline(0.0) ⪧ control_points[begin]
    @test spline(1.0) ⪧ control_points[end]
    @test spline(1e-16) ⪧ control_points[begin] atol = 1e-9
    @test spline(1 - 1e-16) ⪧ control_points[end] atol = 1e-9
    @test periodic_spline(0.0) ⪧ periodic_control_points[begin]
    @test periodic_spline(1.0) ⪧ periodic_control_points[end]
    @test periodic_spline(1e-16) ⪧ periodic_control_points[begin] atol = 1e-9
    @test periodic_spline(1 - 1e-16) ⪧ periodic_control_points[end] atol = 1e-9
    @test longer_spline(0.0) ⪧ longer_control_points[begin]
    @test longer_spline(1.0) ⪧ longer_control_points[end]
    @test longer_spline(1e-16) ⪧ longer_control_points[begin] atol = 1e-9
    @test longer_spline(1 - 1e-16) ⪧ longer_control_points[end] atol = 1e-9
    @test quadratic_spline(0.0) ⪧ quadratic_control_points[begin]
    @test quadratic_spline(1.0) ⪧ quadratic_control_points[end]
    @test quadratic_spline(1e-16) ⪧ quadratic_control_points[begin] atol = 1e-9
    @test quadratic_spline(1 - 1e-16) ⪧ quadratic_control_points[end] atol = 1e-9
    @test sextic_spline(0.0) ⪧ sextic_control_points[begin]
    @test sextic_spline(1.0) ⪧ sextic_control_points[end]
    @test sextic_spline(1e-16) ⪧ sextic_control_points[begin] atol = 1e-9
    @test sextic_spline(1 - 1e-16) ⪧ sextic_control_points[end] atol = 1e-9
    @test quadratic_periodic_spline(0.0) ⪧ quadratic_periodic_control_points[begin]
    @test quadratic_periodic_spline(1.0) ⪧ quadratic_periodic_control_points[end]
    @test quadratic_periodic_spline(1e-16) ⪧ quadratic_periodic_control_points[begin] atol = 1e-9
    @test quadratic_periodic_spline(1 - 1e-16) ⪧ quadratic_periodic_control_points[end] atol = 1e-9

    t = 0:0.1:1
    spline_t = [
        [1.3, 5],
        [3.8621000000000003, 4.8233],
        [6.304800000000001, 5.866400000000001],
        [8.5927, 7.7891],
        [10.6904, 10.2512],
        [12.5625, 12.9125],
        [14.1736, 15.4328],
        [15.4883, 17.471899999999998],
        [16.471200000000003, 18.6896],
        [17.0869, 18.745700000000003],
        [17.3, 17.3]
    ]
    periodic_spline_t = [
        [0, 0],
        [0.6525000000000001, 0.11475000000000002],
        [0.9, 0.37800000000000006],
        [0.8775, 0.6682500000000001],
        [0.7169999999999999, 0.8710000000000001],
        [0.515625, 0.953125],
        [0.34800000000000003, 0.934],
        [0.286875, 0.8336250000000002],
        [0.315, 0.657],
        [0.2756250000000001, 0.3858750000000003],
        [6.661338147750935e-16, 6.661338147750937e-16]
    ]
    longer_spline_t = [
        [0.3, 0.3],
        [0.4796, -0.274],
        [0.7528000000000001, -0.4760000000000001],
        [1.0812000000000002, -0.3659999999999999],
        [1.4264000000000001, -0.003999999999999851],
        [1.75, 0.55],
        [1.92, 1.2832],
        [1.4300000000000004, 2.3716],
        [-0.3199999999999982, 4.038399999999999],
        [-3.929999999999996, 6.506799999999997],
        [-9.999999999999995, 9.999999999999995]
    ]
    quadratic_spline_t = [
        [0.1, 0.1],
        [0.21250000000000002, 0.25],
        [0.35, 0.5],
        [0.525, 0.8750000000000001],
        [0.75, 1.4],
        [0.8125, 3.475],
        [0.5, 8.5],
        [-0.5, 12.75],
        [-2.4999999999999982, 12.500000000000004],
        [-5.625, 8.125],
        [-9.999999999999991, 1.7763568394002498e-14]
    ]
    sextic_spline_t = [
        [1, 0],
        [0.9157231187374586, 0.3586302031184757],
        [0.8162235486801132, 0.543152180420363],
        [0.7124296177871137, 0.6707273332765804],
        [0.6037119762372597, 0.7679083660495618],
        [0.48804323753567636, 0.8453156837021978],
        [0.3631721646588902, 0.9067840910151551],
        [0.2246521007365532, 0.9523458140503687],
        [0.0622718120248865, 0.9784464185342435],
        [-0.1472786929037561, 0.9723545852185906],
        [-0.49999999999999956, 0.8660254037844393]
    ]
    quadratic_periodic_spline_t = [
        [1, 0],
        [0.44499999999999995, 0.465],
        [-0.020000000000000073, 0.66],
        [-0.39500000000000013, 0.5849999999999999],
        [-0.6600000000000001, 0.2999999999999999],
        [-0.75, 9.71445146547012e-17],
        [-0.6600000000000001, -0.29999999999999977],
        [-0.3950000000000005, -0.5849999999999997],
        [-0.02000000000000024, -0.66],
        [0.44499999999999945, -0.4650000000000003],
        [0.9999999999999991, -8.881784197001249e-16]
    ]
    @test spline.(t) ⪧ spline_t
    @test periodic_spline.(t) ⪧ periodic_spline_t
    @test longer_spline.(t) ⪧ longer_spline_t
    @test quadratic_spline.(t) ⪧ quadratic_spline_t
    @test sextic_spline.(t) ⪧ sextic_spline_t
    @test quadratic_periodic_spline.(t) ⪧ quadratic_periodic_spline_t
    @inferred spline(0.0)
    @inferred spline(1 / 2)
    @inferred spline(1)

    ## Differentiation  
    t = LinRange(0, 1, 25000)
    f = tuple.(sin.(π .* t), cos.(π .* t))
    spl = BSpline(f)
    for t in t[3:end-2]
        @test spl(t) ⪧ (sin(π * t), cos(π * t)) rtol = 1e-2
        der = DT.differentiate(spl, t)
        @test ⪧(der, (π * cos(π * t), -π * sin(π * t)); atol=1e-2)
    end

    t = LinRange(0, 1, 25000)
    f = tuple.(sin.(2π .* t), cos.(2π .* t))
    f[end] = f[begin]
    spl = BSpline(f)
    ctr = 0
    for t in t[3:end-2]
        @test spl(t) ⪧ (sin(2π * t), cos(2π * t)) rtol = 1e-2
        der = DT.differentiate(spl, t)
        @test ⪧(der, (2π * cos(2π * t), -2π * sin(2π * t)); atol=1e-2)
    end

    t = LinRange(0, 1, 1500)
    for spl in (spline, sextic_spline, periodic_spline, quadratic_periodic_spline, longer_spline, quadratic_spline)
        for t in t
            der1 = ForwardDiff.derivative(t -> slow_eval_bspline(spl.control_points, spl.knots, t), t)
            der2 = DT.differentiate(spl, t)
            @inferred DT.differentiate(spl, t)
            h = 1e-4
            if h < t < 1 - h
                der3 = (spl(t + h) .- spl(t - h)) ./ (2h)
            elseif t < h
                der3 = (spl(t + h) .- spl(t)) ./ h
            else
                der3 = (spl(t) .- spl(t - h)) ./ h
            end
            flag1 = ⪧(der1, der2, rtol=1e-6, atol=1e-6)
            flag2 = ⪧(der2, der3, rtol=1e-3, atol=1e-6)
            @test flag1 || flag2
        end
    end

    ## Twice differentiate 
    for spl in (spline, sextic_spline, periodic_spline, quadratic_periodic_spline, longer_spline, quadratic_spline)
        for t in LinRange(0, 1, 1500)
            der1 = DT.twice_differentiate(spl, t)
            @inferred DT.twice_differentiate(spl, t)
            h = 1e-6
            if 3h < t < 1 - 3h
                der2 = (spl(t + h) .- 2 .* spl(t) .+ spl(t - h)) ./ (h^2)
            elseif t < 3h
                # h = 1e-2
                der2 = (2.0 .* spl(t) .- 5.0 .* spl(t + h) .+ 4.0 .* spl(t + 2h) .- spl(t + 3h)) ./ h^2
            else
                #  h = 1e-2
                der2 = (2.0 .* spl(t) .- 5.0 .* spl(t - h) .+ 4.0 .* spl(t - 2h) .- spl(t - 3h)) ./ h^2
            end
            @test der1 ⪧ der2 rtol = 1e-1 atol = 1e-1
        end
    end

    t = LinRange(0, 1, 25000)
    f = tuple.(sin.(π .* t), cos.(π .* t))
    spl = BSpline(f)
    for t in t[4:end-4]
        der = DT.twice_differentiate(spl, t)
        @test ⪧(der, (-π^2 * sin(π * t), -π^2 * cos(π * t)); atol=1e-2)
    end

    ## Thrice differentiate 
    t = LinRange(0, 1, 15000)
    for spl in (spline, sextic_spline, periodic_spline, quadratic_periodic_spline, longer_spline, quadratic_spline)
        _der1 = t -> ForwardDiff.derivative(t -> ForwardDiff.derivative(t -> ForwardDiff.derivative(t -> slow_eval_bspline(spl.control_points, spl.knots, t), t), t), t)
        for t in t
            der1 = _der1(t)
            der2 = DT.thrice_differentiate(spl, t)
            @inferred DT.thrice_differentiate(spl, t)
            h = 1e-4
            if h < t < 1 - h
                der3 = (DT.twice_differentiate(spl, t + h) .- DT.twice_differentiate(spl, t - h)) ./ (2h)
            elseif t < h
                der3 = (DT.twice_differentiate(spl, t + h) .- DT.twice_differentiate(spl, t)) ./ h
            else
                der3 = (DT.twice_differentiate(spl, t) .- DT.twice_differentiate(spl, t - h)) ./ h
            end
            flag1 = ⪧(der1, der2, rtol=1e-6, atol=1e-6)
            flag2 = ⪧(der2, der3, rtol=1e-3, atol=1e-6)
            @test flag1 || flag2
        end
    end

    t = LinRange(0, 1, 2500)
    f = tuple.(sin.(π .* t), cos.(π .* t))
    spl = BSpline(f)
    for t in t[4:end-4]
        der = DT.thrice_differentiate(spl, t)
        @test ⪧(der, (-π^3 * cos(π * t), π^3 * sin(π * t)); atol=1e-2, rtol=1e-2)
    end

    ## Closest point 
    for spl in (spline, periodic_spline, longer_spline, quadratic_spline, sextic_spline, quadratic_periodic_spline)
        for _ in 1:500
            p = rand(2) * 30 |> Tuple
            t′, q = DT.get_closest_point(spl, p)
            @test q ⪧ spl(t′)
            _t, _q = closest_point_on_curve(spl, p)
            if !(spl == periodic_spline || spl == quadratic_periodic_spline)
                @test _t ≈ t′ rtol = 1e-1 atol = 1e-1
            elseif t ≠ 0 && t ≠ 1
                @test _t ≈ t′ rtol = 1e-1 atol = 1e-1
            end
            @test q ⪧ _q rtol = 1e-1 atol = 1e-1
            @test DT.dist(p, _q) ≈ DT.dist(p, q) rtol = 1e-1 atol = 1e-1
        end
        for t in LinRange(0, 1, 150)
            p = spl(t)
            t′, q = DT.get_closest_point(spl, p)
            @test q ⪧ spl(t′)
            if !(spl == periodic_spline || spl == quadratic_periodic_spline)
                @test t ≈ t′ rtol = 1e-1 atol = 1e-1
            elseif t ≠ 0 && t ≠ 1
                @test t ≈ t′ rtol = 1e-1 atol = 1e-1
            end
            @test p ⪧ q rtol = 1e-1 atol = 1e-1
        end
        @test DT.get_closest_point(spl, spl(0.0)) == (0.0, spl(0.0))
        if spl ∉ (periodic_spline, quadratic_periodic_spline)
            @test DT.get_closest_point(spl, spl(1.0)) == (1.0, spl(1.0))
        else
            @test DT.get_closest_point(spl, spl(1.0)) == (0.0, spl(0.0))
        end
    end

    ## Sideof 
    d1 = (5.0, 6.0) # left 
    d2 = (1.3, 5.1) # left 
    d3 = (1.3, 4.9) # right 
    d4 = (10.0, 15.0) # left
    d5 = (1.3, 5.0) # on 
    d6 = (15.0, 15.0) # right 
    d7 = (17.0, 17.0) # right
    d8 = (17.5, 17.0) # left
    d9 = (10.0, 10.0) # left 
    d10 = (10.0, 9.0) # right
    cert1 = DT.point_position_relative_to_curve(spline, d1)
    cert2 = DT.point_position_relative_to_curve(spline, d2)
    cert3 = DT.point_position_relative_to_curve(spline, d3)
    cert4 = DT.point_position_relative_to_curve(spline, d4)
    cert5 = DT.point_position_relative_to_curve(spline, d5)
    cert6 = DT.point_position_relative_to_curve(spline, d6)
    cert7 = DT.point_position_relative_to_curve(spline, d7)
    cert8 = DT.point_position_relative_to_curve(spline, d8)
    cert9 = DT.point_position_relative_to_curve(spline, d9)
    cert10 = DT.point_position_relative_to_curve(spline, d10)
    @test DT.is_left(cert1) && DT.is_left(cert2) && DT.is_right(cert3) && DT.is_left(cert4) && DT.is_on(cert5) && DT.is_right(cert6) && DT.is_right(cert7) && DT.is_left(cert8) && DT.is_left(cert9) && DT.is_right(cert10)

    d1 = (0.5, 0.5) # left 
    d2 = (0.2, 0.2) # left
    d3 = (0.4, 0.0) # right 
    d4 = (0.0, 0.5) # right
    d5 = (0.0, 0.0) # on 
    cert1 = DT.point_position_relative_to_curve(periodic_spline, d1)
    cert2 = DT.point_position_relative_to_curve(periodic_spline, d2)
    cert3 = DT.point_position_relative_to_curve(periodic_spline, d3)
    cert4 = DT.point_position_relative_to_curve(periodic_spline, d4)
    cert5 = DT.point_position_relative_to_curve(periodic_spline, d5)
    @test DT.is_left(cert1) && DT.is_left(cert2) && DT.is_right(cert3) && DT.is_right(cert4) && DT.is_on(cert5)
    t = LinRange(0, 1, 15000) |> collect
    pop!(t)
    points = periodic_spline.(t)
    boundary_nodes = [eachindex(t); 1]
    reverse_periodic_spline = BSpline(reverse(periodic_spline.control_points))
    for _ in 1:10000
        p = (rand(), rand())
        δ = DT.distance_to_polygon(p, points, boundary_nodes)
        cert = DT.point_position_relative_to_curve(periodic_spline, p)
        cert2 = DT.point_position_relative_to_curve(reverse_periodic_spline, p)
        if δ > 0
            @test DT.is_left(cert)
            @test DT.is_right(cert2)
        elseif δ < 0
            @test DT.is_right(cert)
            @test DT.is_left(cert2)
        else
            @test DT.is_on(cert)
            @test DT.is_on(cert2)
        end
    end

    d1 = (0.0, 0.0) # right 
    d2 = (1.0, 0.0) # left
    d3 = (1.0, 1.0) # left 
    d4 = (5.0, 3.0) # right
    d5 = (0.0, 5.0) # right
    d6 = (-5.0, 5.0) # left
    d7 = (-3.0, 6.0) # right
    d8 = (-3.2, 6.0) # left
    d9 = (-10.0, 10.0) # on 
    d10 = (0.3, 0.3) # on
    cert1 = DT.point_position_relative_to_curve(longer_spline, d1)
    cert2 = DT.point_position_relative_to_curve(longer_spline, d2)
    cert3 = DT.point_position_relative_to_curve(longer_spline, d3)
    cert4 = DT.point_position_relative_to_curve(longer_spline, d4)
    cert5 = DT.point_position_relative_to_curve(longer_spline, d5)
    cert6 = DT.point_position_relative_to_curve(longer_spline, d6)
    cert7 = DT.point_position_relative_to_curve(longer_spline, d7)
    cert8 = DT.point_position_relative_to_curve(longer_spline, d8)
    cert9 = DT.point_position_relative_to_curve(longer_spline, d9)
    cert10 = DT.point_position_relative_to_curve(longer_spline, d10)
    @test DT.is_right(cert1) && DT.is_left(cert2) && DT.is_left(cert3) && DT.is_right(cert4) && DT.is_right(cert5) && DT.is_left(cert6) && DT.is_right(cert7) && DT.is_left(cert8) && DT.is_on(cert9) && DT.is_on(cert10)

    d1 = (0.0, 5.0) # left 
    d2 = (1.0, 5.0) # right 
    d3 = (-2.0, 15.0) # right
    d4 = (0.1, 0.1) # on 
    d5 = (-10.0, 0.0) # on 
    d6 = (-5.0, 5.0) # left 
    d7 = (-5.0, 10.0) # right 
    d8 = (-10.0, 5.0) # right 
    d9 = (-10.0, -0.01) # left 
    d10 = (-10.0, 0.01) # right
    cert1 = DT.point_position_relative_to_curve(quadratic_spline, d1)
    cert2 = DT.point_position_relative_to_curve(quadratic_spline, d2)
    cert3 = DT.point_position_relative_to_curve(quadratic_spline, d3)
    cert4 = DT.point_position_relative_to_curve(quadratic_spline, d4)
    cert5 = DT.point_position_relative_to_curve(quadratic_spline, d5)
    cert6 = DT.point_position_relative_to_curve(quadratic_spline, d6)
    cert7 = DT.point_position_relative_to_curve(quadratic_spline, d7)
    cert8 = DT.point_position_relative_to_curve(quadratic_spline, d8)
    cert9 = DT.point_position_relative_to_curve(quadratic_spline, d9)
    cert10 = DT.point_position_relative_to_curve(quadratic_spline, d10)
    @test DT.is_left(cert1) && DT.is_right(cert2) && DT.is_right(cert3) && DT.is_on(cert4) && DT.is_on(cert5) && DT.is_left(cert6) && DT.is_right(cert7) && DT.is_right(cert8) && DT.is_left(cert9) && DT.is_right(cert10)

    d1 = (0.9, 0.0) # left 
    d2 = (0.6, 0.5) # left 
    d3 = (0.5, 1.0) # right 
    d4 = (1.0, 0.0) # on 
    d5 = (-0.3, 0.5) # left 
    d6 = (-0.3, 1.0) # right
    d7 = (0.3, 0.5) # left 
    d8 = (0.6, 0.8) # right
    cert1 = DT.point_position_relative_to_curve(sextic_spline, d1)
    cert2 = DT.point_position_relative_to_curve(sextic_spline, d2)
    cert3 = DT.point_position_relative_to_curve(sextic_spline, d3)
    cert4 = DT.point_position_relative_to_curve(sextic_spline, d4)
    cert5 = DT.point_position_relative_to_curve(sextic_spline, d5)
    cert6 = DT.point_position_relative_to_curve(sextic_spline, d6)
    cert7 = DT.point_position_relative_to_curve(sextic_spline, d7)
    cert8 = DT.point_position_relative_to_curve(sextic_spline, d8)
    @test DT.is_left(cert1) && DT.is_left(cert2) && DT.is_right(cert3) && DT.is_on(cert4) && DT.is_left(cert5) && DT.is_right(cert6) && DT.is_left(cert7) && DT.is_right(cert8)

    t = LinRange(0, 1, 15000) |> collect
    pop!(t)
    points = quadratic_periodic_spline.(t)
    boundary_nodes = [eachindex(t); 1]
    reverse_quadratic_periodic_spline = BSpline(reverse(quadratic_periodic_spline.control_points))
    reverse_points = reverse_quadratic_periodic_spline.(t)
    reverse_boundary_nodes = [eachindex(t); 1]
    for i in 1:10000
        rng = StableRNG(i)
        p = (2rand(rng) - 1, 2rand(rng) - 1)
        δ = DT.distance_to_polygon(p, points, boundary_nodes)
        δ2 = DT.distance_to_polygon(p, reverse_points, reverse_boundary_nodes)
        cert = DT.point_position_relative_to_curve(quadratic_periodic_spline, p)
        cert2 = DT.point_position_relative_to_curve(reverse_quadratic_periodic_spline, p)
        if δ > 0
            @test DT.is_left(cert)
        elseif δ < 0
            @test DT.is_right(cert)
        else
            @test DT.is_on(cert)
        end
        if δ2 < 0
            @test DT.is_left(cert2)
        elseif δ2 > 0
            @test DT.is_right(cert2)
        else
            @test DT.is_on(cert2)
        end
    end

    ## Arclength
    for spl in (spline, periodic_spline, longer_spline, quadratic_spline, sextic_spline, quadratic_periodic_spline)
        @test DT.arc_length(spl) ≈ slow_arc_length(spl, 0, 1) rtol = 1e-4
        @inferred DT.arc_length(spl)
        @test DT.arc_length(spl, 0.1, 0.15) ≈ slow_arc_length(spl, 0.1, 0.15) rtol = 1e-4
        @test DT.arc_length(spl, 0.0, 0.0) ≈ 0.0
        @test DT.arc_length(spl, 1.0, 1.0) ≈ 0.0
        @test DT.arc_length(spl, 0.0, 1.0) ≈ slow_arc_length(spl, 0.0, 1.0) rtol = 1e-3
        @test DT.arc_length(spl, 0.3, 0.9) ≈ slow_arc_length(spl, 0.3, 0.9) rtol = 1e-3
        @test DT.arc_length(spl, 0.2, 0.21) ≈ slow_arc_length(spl, 0.2, 0.21) rtol = 1e-4
        for _ in 1:1000
            t₁, t₂ = rand(2)
            t₁, t₂ = minmax(t₁, t₂)
            @test DT.arc_length(spl, t₁, t₂) ≈ slow_arc_length(spl, t₁, t₂) rtol = 1e-2
        end
    end

    ## Curvature 
    for c in (spline, periodic_spline, longer_spline, quadratic_spline, sextic_spline, quadratic_periodic_spline)
        for t in LinRange(0, 1, 100)
            ∂x, ∂y = DT.differentiate(c, t)
            ∂²x, ∂²y = DT.twice_differentiate(c, t)
            vec1 = [∂x, ∂y, 0.0]
            vec2 = [∂²x, ∂²y, 0.0]
            cur1 = DT.curvature(c, t)
            cur2 = cross(vec1, vec2)[3] / norm(vec1)^3
            @test cur1 ≈ cur2 rtol = 1e-5 atol = 1e-5
        end
    end

    ## Total variation 
    for c in (spline, periodic_spline, longer_spline, quadratic_spline, sextic_spline, quadratic_periodic_spline)
        @test c.orientation_markers == DT.orientation_markers(c)
        TV = DT.total_variation(c)
        @inferred DT.total_variation(c)
        TVslow = slow_total_absolute_curvature(c, 0, 1)
        @test TV ≈ TVslow rtol = 1e-3 atol = 1e-3
        for _ in 1:1000
            t1, t2 = minmax(rand(2)...)
            TV = DT.total_variation(c, t1, t2)
            @inferred DT.total_variation(c, t1, t2)
            TVslow = slow_total_absolute_curvature(c, t1, t2)
            @test TV ≈ TVslow rtol = 1e-1 atol = 1e-1
        end
    end

    ## Equidistant split
    for c in (spline, periodic_spline, longer_spline, quadratic_spline, sextic_spline, quadratic_periodic_spline)
        for _ in 1:1000
            t1, t2 = minmax(rand(2)...)
            @inferred DT.get_equidistant_split(c, t1, t2)
            t = DT.get_equidistant_split(c, t1, t2)
            s1 = DT.arc_length(c, t1, t)
            s2 = DT.arc_length(c, t, t2)
            @test s1 ≈ s2 rtol = 1e-1 atol = 1e-1
            @test t1 ≤ t ≤ t2
            @test DT.arc_length(c, t1, t2) ≈ 2s1 rtol = 1e-1
        end
    end

    ## Equivariation split 
    for c in (spline, periodic_spline, longer_spline, quadratic_spline, sextic_spline, quadratic_periodic_spline)
        for _ in 1:1000
            t1, t2 = minmax(rand(2)...)
            @inferred DT.get_equivariation_split(c, t1, t2)
            t, T = DT.get_equivariation_split(c, t1, t2)
            T1 = DT.total_variation(c, t1, t)
            T2 = DT.total_variation(c, t, t2)
            @test T1 ≈ T2 rtol = 1e-1 atol = 1e-1
            @test t1 ≤ t ≤ t2
            @test T ≈ T1
        end
    end

    ## Inverse 
    for c in (spline, periodic_spline, longer_spline, quadratic_spline, sextic_spline, quadratic_periodic_spline)
        for t in LinRange(0, 1, 2500)
            t′ = DT.get_inverse(c, c(t))
            if c ∉ (periodic_spline, quadratic_periodic_spline)
                @test t ≈ t′
            else
                @test c(t) ⪧ c(t′)
            end
        end
    end

    ## get_circle_intersection
    c = sextic_spline
    t₁, t₂, r, = 0.0, 1.0, 1.0
    t′, q′ = DT.get_circle_intersection(c, t₁, t₂, r)
    @test t′ ≈ 0.5096019203840768 && q′ ⪧ (0.47649648275716505, 0.8518859813755716)
    @test DT.dist(c(t₁), q′) ≈ r rtol = 1e-3
    t₁, t₂, r, = 1.0, 0.0, 1.0
    t′, q′ = DT.get_circle_intersection(c, t₁, t₂, r)
    @test t′ ≈ 0.49019803960792163 && q′ ⪧ (0.4997446106449962, 0.8384596166783065)
    @test DT.dist(c(t₁), q′) ≈ r rtol = 1e-2
    t₁, t₂, r = 0.0, 0.1, 0.2
    t′, q′ = DT.get_circle_intersection(c, t₁, t₂, r)
    @test t′ ≈ 0.04410882176435287 && q′ ⪧ (0.9666126168475349, 0.1970440391610145)
    @test DT.dist(c(t₁), q′) ≈ r rtol = 1e-2

    ## ==
    @test spline == spline
    @test spline ≠ periodic_spline
    @test spline ≠ longer_spline
    @test spline ≠ quadratic_spline
    @test spline ≠ sextic_spline
    @test spline ≠ quadratic_periodic_spline
    @test periodic_spline == periodic_spline
    @test periodic_spline ≠ longer_spline
    ctrl1 = [(0.0, 1.3), (17.3, 5.0), (-1.0, 2.0), (50.0, 23.0), (17.3, -2.0), (27.3, 50.1)]
    ctrl2 = [(5.3, 1.3), (17.5, 23.0), (17.3, 200.0), (173.0, 1.3), (0.0, 0.0)]
    @test BSpline(ctrl1) == BSpline(ctrl1)
    @test BSpline(ctrl1) ≠ BSpline(ctrl1, degree=2)
    @test BSpline(ctrl1) ≠ BSpline(ctrl2)
end

@testset "CatmullRomSpline" begin
    @testset "CatmullRomSplineSegment" begin
        ## Evaluation 
        p₀ = (0.0, 0.0)
        p₁ = (1.0, 1.0)
        p₂ = (2.0, -1.0)
        p₃ = (3.0, 0.0)
        α = [0.0, 0.5, 1.0]
        τ = [0.0, 0.5, 1.0]
        spl = [DT.catmull_rom_spline_segment(p₀, p₁, p₂, p₃, α, τ) for α in α, τ in τ]
        for _spl in spl
            @test DT.is_curve_bounded(_spl)
            @test !DT.is_interpolating(_spl)
            @test !DT.is_piecewise_linear(_spl)
            @test _spl(0.0) == p₁
            @test _spl(1.0) == p₂
            @inferred _spl(rand())
            for _ in 1:10000
                t = rand()
                @test _spl(t) ⪧ _spl.a .* t^3 .+ _spl.b .* t^2 .+ _spl.c .* t .+ _spl.d
            end
        end
        t = LinRange(0, 1, 1500)
        fig = Figure(fontsize=44)
        for i in 1:3
            for j in 1:3
                ax = Axis(fig[i, j], xlabel=L"x", ylabel=L"y", title=L"α = %$(α[i]), τ = %$(τ[j])", titlealign=:left, width=400, height=400)
                lines!(ax, spl[i, j].(t), color=:blue, linewidth=6)
                scatter!(ax, [p₀, p₁, p₂, p₃], color=[:blue, :black, :red, :green], markersize=16)
            end
        end
        resize_to_layout!(fig)
        fig
        @test_reference "catmull_rom_segments.png" fig

        ## Differentiate 
        for _spl in spl
            for t in LinRange(0, 1, 15000)
                @test DT.differentiate(_spl, t) ⪧ 3.0 .* _spl.a .* t^2 .+ 2.0 .* _spl.b .* t .+ _spl.c
                @inferred DT.differentiate(_spl, t)
                @test DT.differentiate(_spl, t) ⪧ (_spl(t + 1e-6) .- _spl(t - 1e-6)) ./ 2e-6 atol = 1e-4 rtol = 1e-4
                @test DT.twice_differentiate(_spl, t) ⪧ 6.0 .* _spl.a .* t .+ 2.0 .* _spl.b
                @inferred DT.twice_differentiate(_spl, t)
                @test DT.thrice_differentiate(_spl, t) ⪧ 6.0 .* _spl.a
                @inferred DT.thrice_differentiate(_spl, t)
            end
        end
    end

    ## Extrapolations 
    p₁, p₂, p₃, p₄, x = (2.571, 4.812), (2.05, 17.81), (17.3, -25.3), (0.5, 0.3), 10.0
    perms_1234 =
        [
            4 3 2 1
            4 3 1 2
            4 2 3 1
            4 2 1 3
            4 1 3 2
            4 1 2 3
            3 4 2 1
            3 4 1 2
            3 2 4 1
            3 2 1 4
            3 1 4 2
            3 1 2 4
            2 4 3 1
            2 4 1 3
            2 3 4 1
            2 3 1 4
            2 1 4 3
            2 1 3 4
            1 4 3 2
            1 4 2 3
            1 3 4 2
            1 3 2 4
            1 2 4 3
            1 2 3 4
        ]
    th4 = [DT.thiele4((p₁, p₂, p₃, p₄)[[i, j, k, ℓ]]..., x)[1] for (i, j, k, ℓ) in eachrow(perms_1234)]
    th3 = [DT.thiele3((p₁, p₂, p₃, p₄)[[i, j, k]]..., x)[1] for (i, j, k, ℓ) in eachrow(perms_1234)]
    quad = [DT.quadratic_interp((p₁, p₂, p₃, p₄)[[i, j, k]]..., x)[1] for (i, j, k, ℓ) in eachrow(perms_1234)]
    th4th3quad =
        [
            -12.5908 -31.3945 44.1259
            -12.5908 -91.4788 3.2565
            -12.5908 -31.3945 44.1259
            -12.5908 1.95457 -1214.16
            -12.5908 -91.4788 3.2565
            -12.5908 1.95457 -1214.16
            -12.5908 -31.3945 44.1259
            -12.5908 -91.4788 3.2565
            -12.5908 -31.3945 44.1259
            -12.5908 -22.2832 -91.8257
            -12.5908 -91.4788 3.2565
            -12.5908 -22.2832 -91.8257
            -12.5908 -31.3945 44.1259
            -12.5908 1.95457 -1214.16
            -12.5908 -31.3945 44.1259
            -12.5908 -22.2832 -91.8257
            -12.5908 1.95457 -1214.16
            -12.5908 -22.2832 -91.8257
            -12.5908 -91.4788 3.2565
            -12.5908 1.95457 -1214.16
            -12.5908 -91.4788 3.2565
            -12.5908 -22.2832 -91.8257
            -12.5908 1.95457 -1214.16
            -12.5908 -22.2832 -91.8257
        ]
    @test [th4 th3 quad] ≈ th4th3quad rtol = 1e-5

    control_points = [(-9.0, 5.0), (-7.0, 6.0), (-6.0, 4.0), (-3.0, 5.0)]
    periodic_control_points = [(-2.0, -2.0), (0.0, 0.0), (15.0, 0.0), (13.7, 5.3), (-10.0, 0.0), (-2.0, -2.0)]
    spl = CatmullRomSpline(control_points)
    pspl = CatmullRomSpline(periodic_control_points)
    @test !DT.is_piecewise_linear(spl)
    @test DT.is_interpolating(spl)
    @test spl.left ⪧ (-11.0, 4.7894736842105265)
    @test spl.right ⪧ (0.0, 5.243243243243242)
    @test pspl.right ⪧ (0.0, 0.0)
    @test pspl.left ⪧ (-10.0, 0.0)

    ## get_segment
    control_points = [(-9.0, 5.0), (-7.0, 6.0), (-6.0, 4.0), (-3.0, 5.0)]
    periodic_control_points = [(-2.0, -2.0), (0.0, 0.0), (15.0, 0.0), (13.7, 5.3), (-10.0, 0.0), (-2.0, -2.0)]
    spl = CatmullRomSpline(control_points)
    pspl = CatmullRomSpline(periodic_control_points)
    for spl in (spl, pspl)
        for t in LinRange(0, 1, 15000)
            seg, i = DT.get_segment(spl, t)
            _seg, _i = slow_get_segment(spl.control_points, spl.knots, spl.alpha, spl.tension, t)
            @test seg ⪧ _seg && i == _i
        end
    end

    ## Construction 
    control_points = [(-9.0, 5.0), (-7.0, 6.0), (-6.0, 4.0), (-3.0, 5.0)]
    periodic_control_points = [(-2.0, -2.0), (0.0, 0.0), (15.0, 0.0), (13.7, 5.3), (-10.0, 0.0), (-2.0, -2.0)]
    for alpha in LinRange(0, 1, 10)
        for tension in LinRange(0, 1, 10)
            for lookup_steps in (nothing, 100, 500)
                lookup_steps = isnothing(lookup_steps) ? 5000 : lookup_steps
                if isnothing(lookup_steps)
                    spl = CatmullRomSpline(control_points)
                    pspl = CatmullRomSpline(periodic_control_points)
                    alpha = 1 / 2
                    tension = 0.0
                else
                    spl = CatmullRomSpline(control_points; lookup_steps, _alpha=alpha, _tension=tension)
                    pspl = CatmullRomSpline(periodic_control_points; lookup_steps, _alpha=alpha, _tension=tension)
                end
                knots = zeros(4)
                pknots = zeros(6)
                for i in 2:4
                    knots[i] = knots[i-1] + norm(control_points[i] .- control_points[i-1])^alpha
                end
                for i in 2:6
                    pknots[i] = pknots[i-1] + norm(periodic_control_points[i] .- periodic_control_points[i-1])^alpha
                end
                left = DT.extend_left_control_point(control_points)
                right = DT.extend_right_control_point(control_points)
                pleft = DT.extend_left_control_point(periodic_control_points)
                pright = DT.extend_right_control_point(periodic_control_points)
                knots ./= knots[end]
                pknots ./= pknots[end]
                lookup_table = Vector{NTuple{2,Float64}}(undef, lookup_steps)
                plookup_table = Vector{NTuple{2,Float64}}(undef, lookup_steps)
                for i in 1:lookup_steps
                    t = (i - 1) / (lookup_steps - 1)
                    lookup_table[i] = spl(t)
                    plookup_table[i] = pspl(t)
                end
                @test spl.control_points == control_points
                @test spl.knots == knots
                @test spl.left == left
                @test spl.right == right
                @test spl.lookup_table ⪧ lookup_table
                @test pspl.control_points == periodic_control_points
                @test spl.alpha == alpha
                @test spl.tension == tension
                @test pspl.knots == pknots
                @test pspl.left == pleft
                @test pspl.right == pright
                @test pspl.lookup_table ⪧ plookup_table
                @test pspl.alpha == alpha
                @test pspl.tension == tension
            end
        end
    end

    ## Evaluation 
    control_points = [(-9.0, 5.0), (-7.0, 6.0), (-6.0, 4.0), (-3.0, 5.0), (0.0, 0.0), (2.0, -5.0), (2.0, 10.0), (5.0, 12.0), (-9.0, 10.0)]
    periodic_control_points = [(-2.0, -2.0), (0.0, 0.0), (15.0, 0.0), (13.7, 5.3), (-10.0, 0.0), (-2.0, -2.0)]
    spl = CatmullRomSpline(control_points)
    pspl = CatmullRomSpline(periodic_control_points)
    @test DT.is_curve_bounded(spl)
    @test spl(0) ⪧ spl(1e-16) ⪧ (-9.0, 5.0)
    @test spl(1) ⪧ spl(1 - 1e-16) ⪧ (-9.0, 10.0)
    @test pspl(0) ⪧ pspl(1e-16) ⪧ (-2.0, -2.0)
    @test pspl(1) ⪧ pspl(1 - 1e-16) ⪧ (-2.0, -2.0)
    t_vals = LinRange(0.01, 0.99, 50)
    splt = spl.(t_vals)
    psplt = pspl.(t_vals)
    spltpsplt = [splt psplt]
    @test spltpsplt ⪧ [
        (-8.7329, 5.10981) (-1.74812, -1.89847)
        (-8.18293, 5.45889) (-1.37336, -1.57607)
        (-7.6486, 5.82385) (-1.094, -1.15227)
        (-7.1739, 6.0166) (-0.818099, -0.696327)
        (-6.81369, 5.84448) (-0.45374, -0.277504)
        (-6.59451, 5.29047) (0.0916427, 0.0350795)
        (-6.40898, 4.61491) (0.91763, 0.201985)
        (-6.13871, 4.0949) (2.00909, 0.24424)
        (-5.68115, 3.98246) (3.30874, 0.189587)
        (-5.06189, 4.2405) (4.75931, 0.0657693)
        (-4.35928, 4.64698) (6.30349, -0.0994708)
        (-3.65101, 4.9769) (7.88401, -0.278391)
        (-3.01473, 5.00525) (9.44359, -0.443247)
        (-2.46175, 4.63692) (10.9249, -0.566299)
        (-1.94357, 3.9911) (12.2708, -0.619801)
        (-1.45663, 3.15152) (13.4238, -0.576013)
        (-0.997403, 2.20188) (14.3268, -0.407191)
        (-0.562367, 1.22584) (14.9224, -0.0855925)
        (-0.148, 0.307086) (15.2263, 0.424077)
        (0.250636, -0.560964) (15.3684, 1.11489)
        (0.636968, -1.60017) (15.363, 1.9224)
        (1.00336, -2.70104) (15.2204, 2.78106)
        (1.34141, -3.72216) (14.9507, 3.62534)
        (1.64274, -4.52215) (14.5639, 4.38969)
        (1.89896, -4.95958) (14.0704, 5.00856)
        (2.08801, -4.88026) (13.4597, 5.42096)
        (2.1477, -4.18115) (12.5187, 5.64005)
        (2.1016, -2.9682) (11.2348, 5.69415)
        (1.98501, -1.36638) (9.66666, 5.60243)
        (1.83324, 0.499331) (7.87292, 5.38404)
        (1.68158, 2.50395) (5.91225, 5.05812)
        (1.56535, 4.52249) (3.84332, 4.64383)
        (1.51983, 6.42998) (1.7248, 4.16031)
        (1.58034, 8.10144) (-0.384649, 3.62672)
        (1.78218, 9.41189) (-2.42636, 3.06221)
        (2.17402, 10.2618) (-4.34167, 2.48592)
        (2.86201, 10.8564) (-6.07191, 1.91701)
        (3.69226, 11.2985) (-7.55842, 1.37463)
        (4.45364, 11.6246) (-8.74252, 0.877927)
        (4.93499, 11.8709) (-9.56556, 0.446053)
        (4.93347, 12.0636) (-9.96886, 0.0981565)
        (4.41621, 12.1301) (-9.91019, -0.182886)
        (3.46929, 12.0683) (-9.45274, -0.500459)
        (2.17522, 11.9027) (-8.68383, -0.841491)
        (0.616499, 11.658) (-7.68897, -1.18238)
        (-1.12436, 11.3588) (-6.55369, -1.49952)
        (-2.96485, 11.0298) (-5.36349, -1.76932)
        (-4.82247, 10.6957) (-4.20389, -1.96817)
        (-6.61471, 10.3812) (-3.16042, -2.07246)
        (-8.25906, 10.1109) (-2.31858, -2.05861)
    ] rtol = 1e-5
    for spl in (spl, pspl)
        for (i, t) in enumerate(spl.knots)
            @test spl(t) ⪧ spl.control_points[i]
        end
    end

    ## Differentiation 
    t = LinRange(0, 1, 15000)
    for spl in (spl, pspl)
        for t in t
            der1 = DT.differentiate(spl, t)
            @inferred DT.differentiate(spl, t)
            h = 1e-6
            if h < t < 1 - h
                der2 = (spl(t + h) .- spl(t - h)) ./ (2h)
            elseif t < h
                der2 = (spl(t + h) .- spl(t)) ./ h
            else
                der2 = (spl(t) .- spl(t - h)) ./ h
            end
            @test der1 ⪧ der2 rtol = 1e-4 atol = 1e-4
        end
    end

    ## Twice differentiation 
    t = LinRange(0, 1, 15000)
    for spl in (spl, pspl)
        for t in t
            der1 = DT.twice_differentiate(spl, t)
            @inferred DT.twice_differentiate(spl, t)
            h = 1e-6
            if 3h < t < 1 - 3h
                der2 = (spl(t + h) .- 2 .* spl(t) .+ spl(t - h)) ./ (h^2)
            elseif t < 3h
                # h = 1e-2
                der2 = (2.0 .* spl(t) .- 5.0 .* spl(t + h) .+ 4.0 .* spl(t + 2h) .- spl(t + 3h)) ./ h^2
            else
                #  h = 1e-2
                der2 = (2.0 .* spl(t) .- 5.0 .* spl(t - h) .+ 4.0 .* spl(t - 2h) .- spl(t - 3h)) ./ h^2
            end
            @test der1 ⪧ der2 rtol = 1e-1 atol = 1e-1
        end
    end

    ## Thrice differentiation
    t = LinRange(0, 1, 15000)
    for spl in (spl, pspl)
        for t in t
            h = 1e-6
            der1 = DT.thrice_differentiate(spl, t)
            @inferred DT.thrice_differentiate(spl, t)
            if h < t < 1 - h
                der2 = (DT.twice_differentiate(spl, t + h) .- DT.twice_differentiate(spl, t - h)) ./ (2h)
            elseif t < h
                der2 = (DT.twice_differentiate(spl, t + h) .- DT.twice_differentiate(spl, t)) ./ h
            else
                der2 = (DT.twice_differentiate(spl, t) .- DT.twice_differentiate(spl, t - h)) ./ h
            end
            @test der1 ⪧ der2 rtol = 1e-4 atol = 1e-4
        end
    end

    ## Closest point 
    if !USE_INEXACTPREDICATES
        for spl in (spl, pspl)
            for _ in 1:500
                p = randn(2) * 30 |> Tuple
                t′, q = DT.get_closest_point(spl, p)
                @test q ⪧ spl(t′)
                _t, _q = closest_point_on_curve(spl, p)
                if !(spl == pspl)
                    @test _t ≈ t′ rtol = 1e-1 atol = 1e-1
                elseif t ≠ 0 && t ≠ 1 && t′ ≠ 0 && t′ ≠ 1
                    @test _t ≈ t′ rtol = 1e-1 atol = 1e-1
                end
                @test q ⪧ _q rtol = 1e-1 atol = 1e-1
                @test DT.dist(p, _q) ≈ DT.dist(p, q) rtol = 1e-1 atol = 1e-1
            end
            for t in LinRange(0, 1, 150)
                p = spl(t)
                t′, q = DT.get_closest_point(spl, p)
                @test q ⪧ spl(t′)
                if !(spl == pspl)
                    @test t ≈ t′ rtol = 1e-1 atol = 1e-1
                elseif t ≠ 0 && t ≠ 1
                    @test t ≈ t′ rtol = 1e-1 atol = 1e-1
                end
                @test p ⪧ q rtol = 1e-1 atol = 1e-1
            end
            @test DT.get_closest_point(spl, spl.control_points[1]) == (0.0, spl.control_points[1])
            if spl ≠ pspl
                @test DT.get_closest_point(spl, spl.control_points[end]) == (1.0, spl.control_points[end])
            else
                @test DT.get_closest_point(spl, spl.control_points[end]) == (0.0, spl.control_points[end])
            end
        end
    end

    ## Sideof 
    t = LinRange(0, 1, 15000)
    rev_spl = CatmullRomSpline(reverse(control_points))
    d1 = (-5.0, 5.0) # left 
    d2 = (-5.0, 4.5) # left 
    d3 = (-5.0, 4.0) # right 
    d4 = (0.1, 0.0) # left 
    d5 = (5.0, 0.0) # right 
    d6 = (-5.0, 10.0) # left 
    d7 = (-3.0, 20.0) # right 
    d8 = (0.0, 8.0) # left 
    d9 = (5.0, 10.0) # right 
    d10 = (-5.0, 15.0) # right 
    d11 = (-5.0, -5.0) # right
    cert1 = DT.point_position_relative_to_curve(spl, d1)
    cert2 = DT.point_position_relative_to_curve(spl, d2)
    cert3 = DT.point_position_relative_to_curve(spl, d3)
    cert4 = DT.point_position_relative_to_curve(spl, d4)
    cert5 = DT.point_position_relative_to_curve(spl, d5)
    cert6 = DT.point_position_relative_to_curve(spl, d6)
    cert7 = DT.point_position_relative_to_curve(spl, d7)
    cert8 = DT.point_position_relative_to_curve(spl, d8)
    cert9 = DT.point_position_relative_to_curve(spl, d9)
    cert10 = DT.point_position_relative_to_curve(spl, d10)
    cert11 = DT.point_position_relative_to_curve(spl, d11)
    pcert1 = DT.point_position_relative_to_curve(rev_spl, d1)
    pcert2 = DT.point_position_relative_to_curve(rev_spl, d2)
    pcert3 = DT.point_position_relative_to_curve(rev_spl, d3)
    pcert4 = DT.point_position_relative_to_curve(rev_spl, d4)
    pcert5 = DT.point_position_relative_to_curve(rev_spl, d5)
    pcert6 = DT.point_position_relative_to_curve(rev_spl, d6)
    pcert7 = DT.point_position_relative_to_curve(rev_spl, d7)
    pcert8 = DT.point_position_relative_to_curve(rev_spl, d8)
    pcert9 = DT.point_position_relative_to_curve(rev_spl, d9)
    pcert10 = DT.point_position_relative_to_curve(rev_spl, d10)
    pcert11 = DT.point_position_relative_to_curve(rev_spl, d11)
    @test DT.is_left(cert1) && DT.is_left(cert2) && DT.is_right(cert3) && DT.is_left(cert4) && DT.is_right(cert5) && DT.is_left(cert6) && DT.is_right(cert7) && DT.is_left(cert8) && DT.is_right(cert9) && DT.is_right(cert10) && DT.is_right(cert11)
    @test DT.is_on(DT.point_position_relative_to_curve(spl, spl.control_points[1]))
    @test DT.is_on(DT.point_position_relative_to_curve(spl, spl.control_points[end]))
    @test DT.is_right(pcert1) && DT.is_right(pcert2) && DT.is_left(pcert3) && DT.is_right(pcert4) && DT.is_left(pcert5) && DT.is_right(pcert6) && DT.is_left(pcert7) && DT.is_right(pcert8) && DT.is_left(pcert9) && DT.is_left(pcert10) && DT.is_left(pcert11)
    @test DT.is_on(DT.point_position_relative_to_curve(rev_spl, rev_spl.control_points[1]))
    @test DT.is_on(DT.point_position_relative_to_curve(rev_spl, rev_spl.control_points[end]))

    t = LinRange(0, 1, 15000) |> collect
    pop!(t)
    rev_pspl = CatmullRomSpline(reverse(periodic_control_points))
    points = pspl.(t)
    boundary_nodes = [eachindex(t); 1]
    reverse_points = rev_pspl.(t)
    reverse_boundary_nodes = [eachindex(t); 1]
    for _ in 1:1000
        p = (30randn(), 30randn())
        δ = DT.distance_to_polygon(p, points, boundary_nodes)
        δ2 = DT.distance_to_polygon(p, reverse_points, reverse_boundary_nodes)
        cert = DT.point_position_relative_to_curve(pspl, p)
        cert2 = DT.point_position_relative_to_curve(rev_pspl, p)
        if δ > 0
            @test DT.is_left(cert)
        elseif δ < 0
            @test DT.is_right(cert)
        else
            @test DT.is_on(cert)
        end
        if δ2 < 0
            @test DT.is_left(cert2)
        elseif δ2 > 0
            @test DT.is_right(cert2)
        else
            @test DT.is_on(cert2)
        end
    end

    ## Arclength
    for spl in (spl, pspl)
        @test DT.arc_length(spl) ≈ slow_arc_length(spl, 0, 1) rtol = 1e-4
        @inferred DT.arc_length(spl)
        @test DT.arc_length(spl, 0.1, 0.15) ≈ slow_arc_length(spl, 0.1, 0.15) rtol = 1e-4
        @test DT.arc_length(spl, 0.0, 0.0) ≈ 0.0
        @test DT.arc_length(spl, 1.0, 1.0) ≈ 0.0
        @test DT.arc_length(spl, 0.0, 1.0) ≈ slow_arc_length(spl, 0.0, 1.0) rtol = 1e-3
        @test DT.arc_length(spl, 0.3, 0.9) ≈ slow_arc_length(spl, 0.3, 0.9) rtol = 1e-3
        @test DT.arc_length(spl, 0.2, 0.21) ≈ slow_arc_length(spl, 0.2, 0.21) rtol = 1e-4
        for _ in 1:1000
            t₁, t₂ = rand(2)
            t₁, t₂ = minmax(t₁, t₂)
            @test DT.arc_length(spl, t₁, t₂) ≈ slow_arc_length(spl, t₁, t₂) rtol = 1e-2
        end
    end

    ## Curvature 
    for c in (spl, pspl)
        for t in LinRange(0, 1, 1000)
            ∂x, ∂y = DT.differentiate(c, t)
            ∂²x, ∂²y = DT.twice_differentiate(c, t)
            vec1 = [∂x, ∂y, 0.0]
            vec2 = [∂²x, ∂²y, 0.0]
            cur1 = DT.curvature(c, t)
            cur2 = cross(vec1, vec2)[3] / norm(vec1)^3
            @test cur1 ≈ cur2 rtol = 1e-5 atol = 1e-5
            @inferred DT.curvature(c, t)
        end
    end

    ## Total variation 
    for c in (spl, pspl)
        @test c.orientation_markers == DT.orientation_markers(c)
        TV = DT.total_variation(c)
        @inferred DT.total_variation(c)
        TVslow = slow_total_absolute_curvature(c, 0, 1)
        @test TV ≈ TVslow rtol = 1e-1 atol = 1e-3
        for _ in 1:1000
            t1, t2 = minmax(rand(2)...)
            TV = DT.total_variation(c, t1, t2)
            @inferred DT.total_variation(c, t1, t2)
            TVslow = slow_total_absolute_curvature(c, t1, t2)
            @test TV ≈ TVslow rtol = 1e-1 atol = 1e-1
        end
    end

    ## Equidistant split
    for c in (spl, pspl)
        for _ in 1:1000
            t1, t2 = minmax(rand(2)...)
            @inferred DT.get_equidistant_split(c, t1, t2)
            t = DT.get_equidistant_split(c, t1, t2)
            s1 = DT.arc_length(c, t1, t)
            s2 = DT.arc_length(c, t, t2)
            @test s1 ≈ s2 rtol = 1e-1 atol = 1e-1
            @test t1 ≤ t ≤ t2
            @test DT.arc_length(c, t1, t2) ≈ 2s1 rtol = 1e-1
        end
    end

    ## Equivariation split 
    for c in (spl, pspl)
        for _ in 1:1000
            t1, t2 = minmax(rand(2)...)
            @inferred DT.get_equivariation_split(c, t1, t2)
            t, T = DT.get_equivariation_split(c, t1, t2)
            T1 = DT.total_variation(c, t1, t)
            T2 = DT.total_variation(c, t, t2)
            @test T1 ≈ T2 rtol = 1e-1 atol = 1e-1
            @test t1 ≤ t ≤ t2
            @test T ≈ T1
        end
    end

    # lines 
    c = CatmullRomSpline([(10.0, -3.0), (10.0, 12.0), (10.0, 32.0), (10.0, 64.0)])
    @test c.left == (10.0, -18.0)
    @test c.right == (10.0, 96.0)
    c = CatmullRomSpline([(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0)])
    @test c.left == (-1.0, 0.0)
    @test c.right == (4.0, 0.0)

    ## Inverse 
    for c in (spl, pspl)
        for t in LinRange(0, 1, 25000)
            t′ = DT.get_inverse(c, c(t))
            if c ≠ pspl
                @test t ≈ t′ rtol = 1e-3
            else
                @test c(t) ⪧ c(t′) rtol = 1e-2
            end
        end
    end

    ## get_circle_intersection
    c = spl
    t₁, t₂, r, = 0.0, 1.0, 1.0
    t′, q′ = DT.get_circle_intersection(c, t₁, t₂, r)
    @test t′ ≈ 0.03190638127625525 && q′ ⪧ (-8.130630223009339, 5.495953786102033)
    @test DT.dist(c(t₁), q′) ≈ r rtol = 1e-3
    r = 13.0
    t′, q′ = DT.get_circle_intersection(c, t₁, t₂, r)
    @test t′ ≈ 0.43806261252250456 && q′ ⪧ (1.1435184529348763, -3.13023349324047)
    @test DT.dist(c(t₁), q′) ≈ r rtol = 1e-3
    t₁, t₂, r, = 1.0, 0.0, 1.0
    t′, q′ = DT.get_circle_intersection(c, t₁, t₂, r)
    @test t′ ≈ 0.9865973255686293 && q′ ⪧ (-7.9932950700471235, 10.152688406611302)
    @test DT.dist(c(t₁), q′) ≈ r rtol = 1e-1
    t₁, t₂, r = 0.0, 0.1, 0.2
    t′, q′ = DT.get_circle_intersection(c, t₁, t₂, r)
    @test t′ ≈ 0.007051410282056411 && q′ ⪧ (-8.812636716331587, 5.071072149949431)
    @test DT.dist(c(t₁), q′) ≈ r rtol = 1e-2

    # == 
    @test spl == spl
    ctrl1 = [(10.0, -3.0), (10.0, 12.0), (10.0, 32.0), (10.0, 64.0)]
    ctrl2 = [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0)]
    @test CatmullRomSpline(ctrl1) == CatmullRomSpline(ctrl1)
    @test CatmullRomSpline(ctrl1) ≠ CatmullRomSpline(ctrl2)
    @test CatmullRomSpline(ctrl2) == CatmullRomSpline(ctrl2)
    @test CatmullRomSpline(ctrl1, _alpha=0.3) ≠ CatmullRomSpline(ctrl1)
    @test CatmullRomSpline(ctrl1, _tension=0.2) ≠ CatmullRomSpline(ctrl1)
end

@testset "angle_between" begin
    c₁ = BezierCurve([(1.3, -0.01), (1.2, 0.0), (1.0, 0.1), (0.0, 0.0)])
    c₂ = BezierCurve([(0.0, 0.0), (0.5, 0.005), (1.2, -0.1), (1.3, -0.01)])
    T₁₁ = DT.differentiate(c₁, 0.0)
    T₁₂ = .-DT.differentiate(c₁, 1.0)
    T₂₁ = DT.differentiate(c₂, 0.0)
    T₂₂ = .-DT.differentiate(c₂, 1.0)
    θ₁ = (dot(T₁₁, T₂₂) / (norm(T₁₁) * norm(T₂₂))) |> acos
    θ₂ = (dot(T₁₂, T₂₁) / (norm(T₁₂) * norm(T₂₁))) |> acos
    φ₂ = DT.angle_between(c₁, c₂)
    φ₁ = DT.angle_between(c₂, c₁)
    @test θ₁ ≈ φ₁
    @test θ₂ ≈ φ₂

    p, q = (1.0, 1.0), (0.0, 0.0)
    c₁ = LineSegment(p, q)
    c₂ = CatmullRomSpline([(0.0, 0.0), (1.0, 0.0), (0.0, -1.0), (-1.0, 0.0)])
    c₃ = BSpline([(-1.0, 0.0), (-1.0, 2.0), (-0.2, 1.5), (0.0, 0.5)])
    c₄ = LineSegment((0.0, 0.5), (1.0, 1.0))
    T₁ = (-1.0, -1.0)
    T₂ = (3.3784142300054416, 0.0)
    θ = π - acos(dot(T₁, T₂) / (norm(T₁) * norm(T₂)))
    ϕ = DT.angle_between(c₁, c₂)
    @test θ ≈ ϕ
    T₁ = (3.1328834115332677, -1.2976827983907686)
    T₂ = (0.0, 6.0)
    θ = 2π - acos(dot(T₁, T₂) / (norm(T₁) * norm(T₂)))
    ϕ = DT.angle_between(c₂, c₃)
    @test θ ≈ ϕ
    T₁ = (0.6, -3.0)
    T₂ = (1.0, 0.5)
    θ = π - acos(dot(T₁, T₂) / (norm(T₁) * norm(T₂)))
    ϕ = DT.angle_between(c₃, c₄)
    @test θ ≈ ϕ
    T₁ = (-1.0, -0.5)
    T₂ = (-1.0, -1.0)
    θ = 2π - acos(dot(T₁, T₂) / (norm(T₁) * norm(T₂)))
    ϕ = DT.angle_between(c₄, c₁)
    @test θ ≈ ϕ

    c = CircularArc((1.0, 1.0), (1.0, 1.0), (0.0, 0.0))
    @test DT.angle_between(c, c) ≈ π

    points = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0)]
    boundary1 = [1, 2, 3]
    boundary2 = [3, 1]
    PL1 = DT.PiecewiseLinear(points, boundary1)
    PL2 = DT.PiecewiseLinear(points, boundary2)
    θ = DT.angle_between(PL1, PL2) .|> rad2deg
    @test θ ≈ 45.0

    L1 = LineSegment((10.0, 0.0), (10.0, 10.0))
    L2 = LineSegment((10.0, 10.0), (0.0, 0.0))
    θ1 = DT.angle_between(L1, L2) .|> rad2deg
    @test θ1 ≈ 45.0

    A, B, C, D, E = (6.0, 0.0), (6.0, 8.0), (18.0, 14.0), (0.0, 0.0), (0.0, 6.0)
    L1 = LineSegment(D, A)
    L2 = LineSegment(A, B)
    L3 = LineSegment(B, C)
    L4 = LineSegment(C, E)
    L5 = LineSegment(E, D)
    θ1 = DT.angle_between(L1, L2) |> rad2deg
    @test θ1 ≈ 90.0
    θ1 = DT.angle_between(L2, L3) |> rad2deg
    @test θ1 ≈ 243.434948822922
    θ1 = DT.angle_between(L3, L4) |> rad2deg
    @test θ1 ≈ 2.6025622024998
    θ1 = DT.angle_between(L4, L5) |> rad2deg
    @test θ1 ≈ 113.962488745782
    θ1 = DT.angle_between(L5, L1) |> rad2deg
    @test θ1 ≈ 90.0
    @inferred DT.angle_between(L5, L1)

    A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V = (14.0, -2.0), (12.0, 2.0),
    (0.0, 6.0), (12.0, -2.0), (0.0, -2.0), (2.0, -6.0), (0.0, -8.0),
    (4.0, -8.0), (2.0, -4.0), (6.0, -4.0), (6.0, -14.0), (14.0, -14.0),
    (8.0, -12.0), (14.0, -12.0), (8.0, -8.0), (16.0, -8.0), (8.0, -4.0),
    (16.0, -4.0), (16.0, 6.0), (2.0, 6.0), (8.0, 4.0), (14.0, 4.0)
    LAB = LineSegment(A, B)
    LBC = LineSegment(B, C)
    LCD = LineSegment(C, D)
    LDE = LineSegment(D, E)
    LEF = LineSegment(E, F)
    LFG = LineSegment(F, G)
    LGH = LineSegment(G, H)
    LHI = LineSegment(H, I)
    LIJ = LineSegment(I, J)
    LJK = LineSegment(J, K)
    LKL = LineSegment(K, L)
    LLM = LineSegment(L, M)
    LMN = LineSegment(M, N)
    LNO = LineSegment(N, O)
    LOP = LineSegment(O, P)
    LPQ = LineSegment(P, Q)
    LQR = LineSegment(Q, R)
    LRS = LineSegment(R, S)
    LST = LineSegment(S, T)
    LTU = LineSegment(T, U)
    LUV = LineSegment(U, V)
    LVA = LineSegment(V, A)
    θABC = DT.angle_between(LAB, LBC) |> rad2deg
    θBCD = DT.angle_between(LBC, LCD) |> rad2deg
    θCDE = DT.angle_between(LCD, LDE) |> rad2deg
    θDEF = DT.angle_between(LDE, LEF) |> rad2deg
    θEFG = DT.angle_between(LEF, LFG) |> rad2deg
    θFGH = DT.angle_between(LFG, LGH) |> rad2deg
    θGHI = DT.angle_between(LGH, LHI) |> rad2deg
    θHIJ = DT.angle_between(LHI, LIJ) |> rad2deg
    θIJK = DT.angle_between(LIJ, LJK) |> rad2deg
    θJKL = DT.angle_between(LJK, LKL) |> rad2deg
    θKLM = DT.angle_between(LKL, LLM) |> rad2deg
    θLMN = DT.angle_between(LLM, LMN) |> rad2deg
    θMNO = DT.angle_between(LMN, LNO) |> rad2deg
    θNOP = DT.angle_between(LNO, LOP) |> rad2deg
    θOPQ = DT.angle_between(LOP, LPQ) |> rad2deg
    θPQR = DT.angle_between(LPQ, LQR) |> rad2deg
    θQRS = DT.angle_between(LQR, LRS) |> rad2deg
    θRST = DT.angle_between(LRS, LST) |> rad2deg
    θSTU = DT.angle_between(LST, LTU) |> rad2deg
    θTUV = DT.angle_between(LTU, LUV) |> rad2deg
    θUVA = DT.angle_between(LUV, LVA) |> rad2deg
    θVAB = DT.angle_between(LVA, LAB) |> rad2deg
    ABC = 135.0
    BCD = 15.2551187030578
    CDE = 326.3099324740202
    DEF = 63.434948822922
    EFG = 251.565051177078
    FGH = 45.0
    HIJ = 296.565051177078
    GHI = 63.434948822922
    IJK = 270.0
    JKL = 90.0
    KLM = 18.434948822922
    LMN = 341.565051177078
    MNO = 33.6900675259798
    NOP = 326.3099324740202
    OPQ = 26.565051177078
    PQR = 333.434948822922
    QRS = 90.0
    RST = 90.0
    STU = 18.434948822922
    TUV = 161.565051177078
    UVA = 270.0
    VAB = 333.434948822922
    @test θABC ≈ ABC
    @test θBCD ≈ BCD
    @test θCDE ≈ CDE
    @test θDEF ≈ DEF
    @test θEFG ≈ EFG
    @test θFGH ≈ FGH
    @test θGHI ≈ GHI
    @test θHIJ ≈ HIJ
    @test θIJK ≈ IJK
    @test θJKL ≈ JKL
    @test θKLM ≈ KLM
    @test θLMN ≈ LMN
    @test θMNO ≈ MNO
    @test θNOP ≈ NOP
    @test θOPQ ≈ OPQ
    @test θPQR ≈ PQR
    @test θQRS ≈ QRS
    @test θRST ≈ RST
    @test θSTU ≈ STU
    @test θTUV ≈ TUV
    @test θUVA ≈ UVA
    @test θVAB ≈ VAB

    t = LinRange(0, 1, 2500)
    c₁ = CircularArc((0.0, 0.0), (10.0, 0.0), (5.0, 0.0))
    c₂ = BSpline([(10.0, 0.0), (5.0, 1.0), (0.0, 10.0), (0.0, 8.0), (-5.0, 2.0)])
    c₃ = DT.PiecewiseLinear([(-5.0, 2.0), (-5.0, -4.0), (-7.0, -3.0), (-10.0, -3.0), (0.0, -10.0)], [1, 2, 3, 4, 5])
    c₄ = CatmullRomSpline([(0.0, -10.0), (5.0, -10.0), (10.0, -5.0), (12.0, -3.0)])
    c₅ = BezierCurve([(12.0, -3.0), (11.0, -4.0), (5.0, -7.0), (-1.0, -5.0)])
    c₆ = LineSegment((-1.0, -5.0), (0.0, 0.0))
    θ12 = DT.angle_between(c₁, c₂) |> rad2deg # 101.3099324740202
    θ23 = DT.angle_between(c₂, c₃) |> rad2deg # 140.19442890773482
    θ34 = DT.angle_between(c₃, c₄) |> rad2deg # 145.00797980144134
    θ45 = DT.angle_between(c₄, c₅) |> rad2deg # 4.3375617566185366e-15
    θ56 = DT.angle_between(c₅, c₆) |> rad2deg # 262.8749836510982
    θ61 = DT.angle_between(c₆, c₁) |> rad2deg # 348.6900675259798
    @test θ12 ≈ 101.3099324740202
    @test θ23 ≈ 140.19442890773482
    @test θ34 ≈ 145.00797980144134
    @test θ45 ≈ 4.3375617566185366e-15
    @test θ56 ≈ 262.8749836510982
    @test θ61 ≈ 348.6900675259798
    c₁ = LineSegment((0.0, 0.0), (-1.0, -5.0))
    c₂ = BezierCurve([(-1.0, -5.0), (5.0, -7.0), (11.0, -4.0), (12.0, -3.0)])
    c₃ = CatmullRomSpline([(12.0, -3.0), (12.0, -4.0), (5.0, -10.0), (0.0, -10.0)])
    c₄ = DT.PiecewiseLinear([(0.0, -10.0), (-10.0, -3.0), (-7.0, -3.0), (-5.0, -4.0), (-5.0, 2.0)], [1, 2, 3, 4, 5])
    c₅ = BSpline([(-5.0, 2.0), (0.0, 8.0), (0.0, 10.0), (5.0, 1.0), (10.0, 0.0)])
    c₆ = CircularArc((10.0, 0.0), (0.0, 0.0), (5.0, 0.0), positive=false)
    θ12 = DT.angle_between(c₁, c₂) |> rad2deg
    θ23 = DT.angle_between(c₂, c₃) |> rad2deg
    θ34 = DT.angle_between(c₃, c₄) |> rad2deg
    θ45 = DT.angle_between(c₄, c₅) |> rad2deg
    θ56 = DT.angle_between(c₅, c₆) |> rad2deg
    θ61 = DT.angle_between(c₆, c₁) |> rad2deg
    @test θ12 ≈ 360 - 262.8749836510982
    @test θ23 ≈ 315.0
    @test θ34 ≈ 360 - 145.00797980144134
    @test θ45 ≈ 360 - 140.19442890773482
    @test θ56 ≈ 360 - 101.3099324740202
    @test θ61 ≈ 360 - 348.6900675259798
end