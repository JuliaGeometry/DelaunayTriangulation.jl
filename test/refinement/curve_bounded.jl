using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using ..DelaunayTriangulation: EllipticalArc
using CairoMakie
using DataStructures
using Random
using LinearAlgebra
using ReferenceTests
using Test
using StructEquality
using DelimitedFiles
using StableRNGs

const SACM = DT.SmallAngleComplexMember
const SAC = DT.SmallAngleComplex

# Examples:
#   I: A piecewise linear contiguous section.
#   II: A piecewise linear boundary curve with multiple sections.
#   III: A piecewise linear boundary curve split into multiple curves.
#   IV: A contiguous section made up of a circular arc.
#   V: A contiguous section made up of a Bezier curve.
#   VI: A sectioned curve made up of a circular arc, BSpline, and a piecewise linear curve.
#   VII: A sectioned curve made up of a circular arc and a BSpline.
#   VIII: A sectioned curve made up of a piecewise linear curve, an elliptical arc, a piecewise linear curve, and a Catmull-Rom spline.
#   IX: A multiply-connected domain made up of a (1) piecewise linear curve and (2) an interior circle.
#   X: A multiply-connected domain made up of a (1) piecewise linear curve and an elliptical arc, and (2) a BSpline and a piecewise linear curve.
#   XI: A multiply-connected domain made up of a (1) piecewise linear curve and an elliptical arc, (2) a BSpline, (3) a piecewise linear curve, (4) a Catmull-Rom spline and a Bezier curve, (5) a piecewise linear curve, and (6) a circle
#   XII: A multiply-connected domain that is initially intersecting after the coarse-discretisation
names = (
    "piecewise_contiguous",
    "piecewise_multiple_sections",
    "piecewise_multiple_curves",
    "contiguous_circular_arc",
    "contiguous_bezier",
    "circular_bspline_piecewise_linear_curve",
    "circular_bspline",
    "piecewise_linear_elliptical_piecewise_linear_catmull",
    "multiply_connected_piecewise_linear_interior_circle",
    "multiply_connected_piecewise_linear_elliptical_bspline_piecewise_linear",
    "multiply_connected_piecewise_linear_elliptical_bspline_piecewise_linear_catmull_bezier_piecewise_linear_circle",
    "multiply_connected_intersecting",
)

curve_I = [1, 2, 3, 4, 1]
points_I = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
fpoints_I = flatten_boundary_nodes(points_I, curve_I)
points_I_extra = copy(points_I)
for _ in 1:50
    push!(points_I_extra, (rand(), rand()))
end
fpoints_I_extra = flatten_boundary_nodes(points_I_extra, curve_I)
points_I_extra_segments = copy(points_I_extra)
push!(points_I_extra_segments, (0.1, 0.1), (0.9, 0.9))
n = length(points_I_extra_segments)
segments_I = Set([(n, n - 1)])
fpoints_I_extra_segments = flatten_boundary_nodes(points_I_extra_segments, curve_I, segments_I)

curve_II = [[1, 2, 3, 4, 5], [5, 6, 7, 8, 9], [9, 10, 11, 1]]
points_II = [
    (0.0, 0.0), (0.25, 0.0), (0.5, 0.0), (0.75, 0.0), (1.0, 0.0),
    (1.0, 0.25), (1.0, 0.5), (1.0, 0.75), (1.0, 1.0),
    (0.75, 0.75), (0.25, 0.25),
]
fpoints_II = flatten_boundary_nodes(points_II, curve_II)
points_II_extra = copy(points_II)
push!(points_II_extra, (0.5, 0.2), (0.7, 0.3), (0.1, 0.1), (0.99, 0.99))
fpoints_II_extra = flatten_boundary_nodes(points_II_extra, curve_II)
points_II_extra_segments = copy(points_II_extra)
n = length(points_II_extra_segments)
push!(points_II_extra_segments, (0.55, 0.1), (0.6, 0.2), (0.6, 0.4), (0.67, 0.3), (0.3, 0.05), (0.4, 0.06), (0.9, 0.1))
segments_II = Set([(n + 1, n + 2), (n + 1, n + 3), (n + 1, n + 4), (n + 1, n + 5), (n + 1, n + 6), (n + 1, n + 7)])
fpoints_II_extra_segments = flatten_boundary_nodes(points_II_extra_segments, curve_II, segments_II)

curve_III = [[[1, 2, 3, 4, 5], [5, 6, 7, 8, 9], [9, 10, 11, 1]], [[15, 14, 13, 12], [12, 15]]]
points_III = [
    (0.0, 0.0), (0.25, 0.0), (0.5, 0.0), (0.75, 0.0), (1.0, 0.0),
    (1.0, 0.25), (1.0, 0.5), (1.0, 0.75), (1.0, 1.0),
    (0.0, 1.0), (0.0, 0.5),
    (0.25, 0.25), (0.75, 0.25), (0.75, 0.75), (0.25, 0.75),
]
fpoints_III = flatten_boundary_nodes(points_III, curve_III)
points_III_extra = copy(points_III)
push!(points_III_extra, (0.0001, 0.0001), (0.999, 0.999), (0.6, 0.2), (0.5, 0.9))
fpoints_III_extra = flatten_boundary_nodes(points_III_extra, curve_III)
points_III_extra_segments = copy(points_III_extra)
n = length(points_III_extra_segments)
push!(points_III_extra_segments, (0.2, 0.01), (0.2, 0.99))
segments_III = Set([(n + 2, n + 1)])
fpoints_III_extra_segments = flatten_boundary_nodes(points_III_extra_segments, curve_III, segments_III)

curve_IV = [CircularArc((1.0, 0.0), (1.0, 0.0), (0.0, 0.0))]
points_IV = NTuple{2,Float64}[]
fpoints_IV = flatten_boundary_nodes(points_IV, curve_IV)
points_IV_extra = copy(points_IV)
push!(points_IV_extra, (0.99, 0.14), (0.0, 0.99), (-0.99, 0.0), (0.99, 0.0), (0.0, -0.99), (0.0, 0.0), (0.5, -0.5))
fpoints_IV_extra = flatten_boundary_nodes(points_IV_extra, curve_IV)
points_IV_extra_segments = copy(points_IV_extra)
push!(points_IV_extra_segments, (-0.57, -0.57), (0.0, 0.9291912), (-0.58, 0.8), (-0.42, 0.9), (-0.18, -0.98), (0.04, -0.998))
n = length(points_IV_extra_segments)
segments_IV = Set([(n - 1, n), (n - 2, n - 3), (n - 4, n - 5)])
fpoints_IV_extra_segments = flatten_boundary_nodes(points_IV_extra_segments, curve_IV, segments_IV)

curve_V = [BezierCurve([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)])]
points_V = [(0.0, 0.0), (0.2, 0.25)]
fpoints_V = flatten_boundary_nodes(points_V, curve_V)
points_V_extra = copy(points_V)
push!(points_V_extra, (0.01, 0.01), (0.25, 0.25), (0.5, 0.5), (0.25, 0.55), (0.5, 0.25), (0.6, 0.6))

curve_VI = [
    [CircularArc((1.0, 0.0), (0.0, 1.0), (0.0, 0.0))],
    [BSpline([(0.0, 1.0), (-1.0, 2.0), (-2.0, 0.0), (-2.0, -1.0), (0.0, -2.0)])],
    [5, 6, 10],
]
points_VI = [(0.1, 0.1), (0.15, 0.15), (0.23, 0.23), (0.009, 0.11), (0.0, -2.0), (0.2, -1.7), (0.000591, 0.00019), (0.111, -0.005), (-0.0001, -0.00991), (1.0, 0.0)]
fpoints_VI = flatten_boundary_nodes(points_VI, curve_VI)
points_VI_extra = copy(points_VI)
push!(points_VI_extra, (0.0, 0.0), (0.999, 0.0), (-1.0, -1.0), (0.5, 0.0), (-1.0, 0.0), (-1.0, -1.2))

curve_VII = [
    [CircularArc((2.0, 0.0), (-2.0, 0.0), (0.0, 0.0))],
    [BSpline([(-2.0, 0.0), (-2.0, -1.0), (0.0, -1.0), (1.0, -1.0), (2.0, -1.0), (2.0, 0.0)])],
]
points_VII = [(2.0, 0.0), (0.0, 0.5)]
fpoints_VII = flatten_boundary_nodes(points_VII, curve_VII)
points_VII_extra = copy(points_VII)
push!(points_VII_extra, (0.0, 0.0), (1.0, 0.0), (-1.0, 1.0), (0.0, 1.95), (1.0, 0.8))

curve_VIII = [
    [1, 2, 3, 4, 5],
    [DT.EllipticalArc((0.0, 0.0), (2.0, -2.0), (1.0, -1.0), sqrt(2), sqrt(2), 45.0)],
    [6, 7, 8, 9, 10],
    [CatmullRomSpline([(10.0, -3.0), (20.0, 0.0), (18.0, 0.0), (10.0, 0.0)], lookup_steps=5000)],
]
points_VIII = [
    (10.0, 0.0), (8.0, 0.0), (4.0, 0.0), (2.0, 2.0), (0.0, 0.0), (2.0, -2.0),
    (2.5, -2.0), (3.5, -2.0), (4.5, -3.0), (10.0, -3.0), (10.0, -0.2), (14.0, -0.05),
]
fpoints_VIII = flatten_boundary_nodes(points_VIII, curve_VIII)
points_VIII_extra = copy(points_VIII)
push!(points_VIII_extra, (5.0, -0.01), (10.0, -1.0), (15.0, -1.95), (0.0, -1.0), (2.0, -1.5), (1.0, 0.5), (10.0, -2.0))

curve_IX =
    [
        [
            [1, 2, 3, 4, 5, 6, 7, 1],
        ],
        [
            [CircularArc((0.6, 0.5), (0.6, 0.5), (0.5, 0.5), positive=false)],
        ],
    ]
points_IX = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.5, 1.5), (0.0, 1.0), (0.0, 0.5), (0.0, 0.2)]
fpoints_IX = flatten_boundary_nodes(points_IX, curve_IX)
points_IX_extra = copy(points_IX)
push!(points_IX_extra, (0.5, 0.01), (0.01, 0.01), (0.01, 0.5), (0.01, 1.0), (0.5, 1.49), (0.99, 1.0), (0.99, 0.5), (0.75, 0.5), (0.5, 0.75))

curve_X = [
    [
        [1, 2, 3], [DT.EllipticalArc((2.0, 0.0), (-2.0, 0.0), (0.0, 0.0), 2, 1 / 2, 0.0)],
    ],
    [
        [BSpline(reverse([(1.0, 0.2), (0.0, 0.4), (0.0, 0.3), (-1.0, 0.2)]))], reverse([4, 5, 6, 7, 8]),
    ],
]
points_X = [
    (-2.0, 0.0), (0.0, 0.0), (2.0, 0.0), (-1.0, 0.2), (-1.0, 0.1), (0.0, 0.1), (1.0, 0.1), (1.0, 0.2),
]
fpoints_X = flatten_boundary_nodes(points_X, curve_X)
points_X_extra = copy(points_X)
push!(points_X_extra, (0.0, 0.01), (-0.5, 0.27), (1.0, 0.275), (1.5, 0.2), (0.0, 0.4), (-1.9, 0.01), (1.9, 0.01), (-1.5, 0.2), (-1.0, 0.4), (0.0, 0.45), (0.0, 0.49), (0.0, 0.099))

curve_XI = [
    [
        [1, 2, 3], [DT.EllipticalArc((2.0, 0.0), (-2.0, 0.0), (0.0, 0.0), 2, 1 / 2, 0.0)],
    ],
    [
        [BSpline([(0.0, 0.4), (1.0, 0.2), (0.0, 0.1), (-1.0, 0.2), (0.0, 0.4)])],
    ],
    [
        [4, 5, 6, 7, 4],
    ],
    [
        [BezierCurve(reverse([(-1.0, -3.0), (-1.0, -2.5), (0.0, -2.5), (0.0, -2.0)]))], [CatmullRomSpline(reverse([(0.0, -2.0), (1.0, -3.0), (0.0, -4.0), (-1.0, -3.0)]))],
    ],
    [
        [12, 11, 10, 12],
    ],
    [
        [CircularArc((1.1, -3.0), (1.1, -3.0), (0.0, -3.0), positive=false)],
    ],
]
points_XI = [(-2.0, 0.0), (0.0, 0.0), (2.0, 0.0), (-2.0, -5.0), (2.0, -5.0), (2.0, -1 / 10), (-2.0, -1 / 10), (-1.0, -3.0), (0.0, -4.0), (0.0, -2.3), (-0.5, -3.5), (0.9, -3.0)]
fpoints_XI = flatten_boundary_nodes(points_XI, curve_XI)
points_XI_extra = copy(points_XI)
push!(points_XI_extra, (-1.0, -4.0), (-1.0, -2.0), (0.0, -4.1), (1.0, -4.0), (1.0, -2.0), (1.0, 0.4), (0.0, 0.49), (1.99, -2.0), (-1.99, -2.0), (-1.99, -3.71))

ctrl = [
    (0.0, 0.0), (2.0, 0.0), (1.6, -0.1),
    (0.3, -0.2), (-0.31, -0.35), (-0.2, 1.0),
    (0.0, 0.8), (0.2, 0.6), (0.4, 0.4),
    (2.0, 0.4), (0.0, 0.0),
]
reverse!(ctrl)
points_XII = [
    (-0.1, 0.8), (-0.15, -0.15), (0.3, -0.1),
    (0.0, -0.1), (-0.1, 0.0),
    (0.4, 0.2), (0.2, 0.4), (0.0, 0.6),
]
curve_XII = [[[BSpline(ctrl, lookup_steps=25000)]], [[1, 8, 7, 6, 5, 4, 3, 2, 1]]]
fpoints_XII = flatten_boundary_nodes(points_XII, curve_XII)
points_XII_extra = copy(points_XII)
push!(points_XII_extra, (0.5, -0.15), (1.0, -0.1), (1.0, 0.25), (0.0, 0.75), (-0.01, 0.0))

@testset "multiple_sections/curves" begin
    @test !DT.has_multiple_sections(curve_I) && !DT.has_multiple_curves(curve_I)
    @test DT.has_multiple_sections(curve_II) && !DT.has_multiple_curves(curve_II)
    @test DT.has_multiple_sections(curve_III) && DT.has_multiple_curves(curve_III)
    @test !DT.has_multiple_sections(curve_IV) && !DT.has_multiple_curves(curve_IV)
    @test !DT.has_multiple_sections(curve_V) && !DT.has_multiple_curves(curve_V)
    @test DT.has_multiple_sections(curve_VI) && !DT.has_multiple_curves(curve_VI)
    @test DT.has_multiple_sections(curve_VII) && !DT.has_multiple_curves(curve_VII)
    @test DT.has_multiple_sections(curve_VIII) && !DT.has_multiple_curves(curve_VIII)
    @test DT.has_multiple_sections(curve_IX) && DT.has_multiple_curves(curve_IX)
    @test DT.has_multiple_sections(curve_X) && DT.has_multiple_curves(curve_X) && DT.has_multiple_sections(curve_X[1]) && !DT.has_multiple_curves(curve_X[1]) && DT.has_multiple_sections(curve_X[2])
    @test DT.has_multiple_sections(curve_XI) && DT.has_multiple_curves(curve_XI)
    @test DT.has_multiple_sections(curve_XII) && DT.has_multiple_curves(curve_XII)
end

@testset "is_curve_bounded" begin
    @test !DT.is_curve_bounded(())
    @test !DT.is_curve_bounded(Int[])
    @test !DT.is_curve_bounded(1)
    @test DT.is_curve_bounded((CircularArc((1.0, 0.0), (0.0, 1.0), (0.0, 0.0)),))
    curve1 = [CircularArc((1.0, 0.0), (0.0, 1.0), (0.0, 0.0))]
    curve2 = [1, 2, 3, 4, 5, 6, 7, 1]
    curve3 = [[1, 2, 3, 4, 5, 6, 7, 8], [8, 9, 10]]
    curve4 = [[[1, 2, 3, 4, 5, 6], [1, 2, 4, 5, 1]], [[1, 2, 3, 4, 5, 6]]]
    curve5 = [[CircularArc((1.0, 0.0), (0.0, 1.0), (0.0, 0.0))]]
    curve6 = [[[CircularArc((1.0, 0.0), (0.0, 1.0), (0.0, 0.0))]]]
    @test DT.is_curve_bounded(curve1)
    @test !DT.is_curve_bounded(curve2)
    @test !DT.is_curve_bounded(curve3)
    @test !DT.is_curve_bounded(curve4)
    @test DT.is_curve_bounded(curve5)
    @test DT.is_curve_bounded(curve6)

    @test !DT.is_curve_bounded(curve_I)
    @test !DT.is_curve_bounded(curve_II)
    @test !DT.is_curve_bounded(curve_III)
    @test DT.is_curve_bounded(curve_IV)
    @test DT.is_curve_bounded(curve_V)
    @test DT.is_curve_bounded(curve_VI)
    @test DT.is_curve_bounded(curve_VII)
    @test DT.is_curve_bounded(curve_VIII)
    @test DT.is_curve_bounded(curve_IX)
    @test DT.is_curve_bounded(curve_X)
    @test DT.is_curve_bounded(curve_XI)
end

@testset "to_boundary_curves" begin
    curves_I_tuple = DT.to_boundary_curves(points_I, curve_I)
    @test curves_I_tuple == (DT.PiecewiseLinear(points_I, curve_I),)

    curves_II_tuple = DT.to_boundary_curves(points_II, curve_II)
    @test curves_II_tuple == (DT.PiecewiseLinear(points_II, curve_II[1]), DT.PiecewiseLinear(points_II, curve_II[2]), DT.PiecewiseLinear(points_II, curve_II[3]))

    curves_III_tuple = DT.to_boundary_curves(points_III, curve_III)
    @test curves_III_tuple == (DT.PiecewiseLinear(points_III, curve_III[1][1]), DT.PiecewiseLinear(points_III, curve_III[1][2]), DT.PiecewiseLinear(points_III, curve_III[1][3]), DT.PiecewiseLinear(points_III, curve_III[2][1]), DT.PiecewiseLinear(points_III, curve_III[2][2]))

    curves_IV_tuple = DT.to_boundary_curves(points_IV, curve_IV)
    @test curves_IV_tuple == (curve_IV[1],)

    curves_V_tuple = DT.to_boundary_curves(points_V, curve_V)
    @test curves_V_tuple == (curve_V[1],)

    curves_VI_tuple = DT.to_boundary_curves(points_VI, curve_VI)
    @test curves_VI_tuple == (curve_VI[1][1], curve_VI[2][1], DT.PiecewiseLinear(points_VI, curve_VI[3]))

    curves_VII_tuple = DT.to_boundary_curves(points_VII, curve_VII)
    @test curves_VII_tuple == (curve_VII[1][1], curve_VII[2][1])

    curves_VIII_tuple = DT.to_boundary_curves(points_VIII, curve_VIII)
    @test curves_VIII_tuple == (DT.PiecewiseLinear(points_VIII, curve_VIII[1]), curve_VIII[2][1], DT.PiecewiseLinear(points_VIII, curve_VIII[3]), curve_VIII[4][1])

    curves_IX_tuple = DT.to_boundary_curves(points_IX, curve_IX)
    @test curves_IX_tuple == (DT.PiecewiseLinear(points_IX, curve_IX[1][1]), curve_IX[2][1][1])

    curves_X_tuple = DT.to_boundary_curves(points_X, curve_X)
    @test curves_X_tuple == (DT.PiecewiseLinear(points_X, curve_X[1][1]), curve_X[1][2][1], curve_X[2][1][1], DT.PiecewiseLinear(points_X, curve_X[2][2]))

    curves_XI_tuple = DT.to_boundary_curves(points_XI, curve_XI)
    @test curves_XI_tuple == (
        DT.PiecewiseLinear(points_XI, curve_XI[1][1]), curve_XI[1][2][1], curve_XI[2][1][1], DT.PiecewiseLinear(points_XI, curve_XI[3][1]), curve_XI[4][1][1], curve_XI[4][2][1],
        DT.PiecewiseLinear(points_XI, curve_XI[5][1]), curve_XI[6][1][1],
    )
end

@testset "get_skeleton" begin
    curve_I_skeleton = DT.get_skeleton(curve_I, Int)
    @inferred DT.get_skeleton(curve_I, Int)
    @test curve_I_skeleton ⊢ Int[]

    curve_II_skeleton = DT.get_skeleton(curve_II, Int32)
    @inferred DT.get_skeleton(curve_II, Int32)
    @test !(curve_II_skeleton ⊢ [Int[], Int[], Int[]])
    @test curve_II_skeleton ⊢ [Int32[], Int32[], Int32[]]

    curve_III_skeleton = DT.get_skeleton(curve_III, Int)
    @inferred DT.get_skeleton(curve_III, Int)
    @test curve_III_skeleton ⊢ [[Int[], Int[], Int[]], [Int[], Int[]]]

    curve_IV_skeleton = DT.get_skeleton(curve_IV, Int)
    @test curve_IV_skeleton ⊢ Int[]

    curve_V_skeleton = DT.get_skeleton(curve_V, Int)
    @test curve_V_skeleton ⊢ Int[]

    curve_VI_skeleton = DT.get_skeleton(curve_VI, Int)
    @test curve_VI_skeleton ⊢ [Int[], Int[], Int[]]

    curve_VII_skeleton = DT.get_skeleton(curve_VII, Int)
    @test curve_VII_skeleton ⊢ [Int[], Int[]]

    curve_VIII_skeleton = DT.get_skeleton(curve_VIII, Int32)
    @test curve_VIII_skeleton ⊢ [Int32[], Int32[], Int32[], Int32[]]

    curve_IX_skeleton = DT.get_skeleton(curve_IX, Int)
    @test curve_IX_skeleton ⊢ [[Int[]], [Int[]]]

    curve_X_skeleton = DT.get_skeleton(curve_X, Int)
    @test curve_X_skeleton ⊢ [[Int[], Int[]], [Int[], Int[]]]

    curve_XI_skeleton = DT.get_skeleton(curve_XI, Int32)
    @test curve_XI_skeleton ⊢ [[Int32[], Int32[]], [Int32[]], [Int32[]], [Int32[], Int32[]], [Int32[]], [Int32[]]]
end

@testset "convert_boundary_curves!" begin
    points = deepcopy(points_I)
    boundary_nodes = deepcopy(curve_I)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int)
    @test boundary_curves == (DT.PiecewiseLinear(points, boundary_nodes),) && new_boundary_nodes == boundary_nodes

    points = deepcopy(points_II)
    boundary_nodes = deepcopy(curve_II)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int32)
    @test boundary_curves == (DT.PiecewiseLinear(points, boundary_nodes[1]), DT.PiecewiseLinear(points, boundary_nodes[2]), DT.PiecewiseLinear(points, boundary_nodes[3])) && new_boundary_nodes == boundary_nodes

    points = deepcopy(points_III)
    boundary_nodes = deepcopy(curve_III)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int)
    @test boundary_curves == (DT.PiecewiseLinear(points, boundary_nodes[1][1]), DT.PiecewiseLinear(points, boundary_nodes[1][2]), DT.PiecewiseLinear(points, boundary_nodes[1][3]), DT.PiecewiseLinear(points, boundary_nodes[2][1]), DT.PiecewiseLinear(points, boundary_nodes[2][2])) && new_boundary_nodes == boundary_nodes

    points = deepcopy(points_IV)
    boundary_nodes = deepcopy(curve_IV)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int)
    @test boundary_curves == DT.to_boundary_curves(points, boundary_nodes) && new_boundary_nodes == [1, 1] && points == [(1.0, 0.0)]

    points = deepcopy(points_V)
    boundary_nodes = deepcopy(curve_V)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int)
    @test boundary_curves == DT.to_boundary_curves(points, boundary_nodes) && new_boundary_nodes == [1, 1] && points == [(0.0, 0.0), (0.2, 0.25)]

    points = deepcopy(points_VI)
    boundary_nodes = deepcopy(curve_VI)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int)
    @test boundary_curves == DT.to_boundary_curves(points, boundary_nodes) &&
          new_boundary_nodes == [[10, 11], [11, 5], [5, 6, 10]] &&
          points == [(0.1, 0.1), (0.15, 0.15), (0.23, 0.23), (0.009, 0.11), (0.0, -2.0), (0.2, -1.7), (0.000591, 0.00019), (0.111, -0.005), (-0.0001, -0.00991), (1.0, 0.0), (0.0, 1.0)]

    points = deepcopy(points_VII)
    boundary_nodes = deepcopy(curve_VII)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int)
    @test boundary_curves == DT.to_boundary_curves(points, boundary_nodes) &&
          new_boundary_nodes == [[1, 3], [3, 1]] &&
          points == [(2.0, 0.0), (0.0, 0.5), (-2.0, 0.0)]

    points = deepcopy(points_VIII)
    boundary_nodes = deepcopy(curve_VIII)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int)
    @test boundary_curves == DT.to_boundary_curves(points, boundary_nodes) &&
          new_boundary_nodes == [[1, 2, 3, 4, 5], [5, 6], [6, 7, 8, 9, 10], [10, 1]] &&
          points == [
              (10.0, 0.0), (8.0, 0.0), (4.0, 0.0), (2.0, 2.0), (0.0, 0.0), (2.0, -2.0),
              (2.5, -2.0), (3.5, -2.0), (4.5, -3.0), (10.0, -3.0), (10.0, -0.2), (14.0, -0.05),
          ]

    points = deepcopy(points_IX)
    boundary_nodes = deepcopy(curve_IX)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int)
    @test boundary_curves == DT.to_boundary_curves(points, boundary_nodes) &&
          new_boundary_nodes == [[[1, 2, 3, 4, 5, 6, 7, 1]], [[8, 8]]] &&
          points == [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.5, 1.5), (0.0, 1.0), (0.0, 0.5), (0.0, 0.2), (0.6, 0.5)]

    points = deepcopy(points_X)
    boundary_nodes = deepcopy(curve_X)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int)
    @test boundary_curves == DT.to_boundary_curves(points, boundary_nodes) &&
          new_boundary_nodes == [[[1, 2, 3], [3, 1]], [[4, 8], [8, 7, 6, 5, 4]]] &&
          points == [(-2.0, 0.0), (0.0, 0.0), (2.0, 0.0), (-1.0, 0.2), (-1.0, 0.1), (0.0, 0.1), (1.0, 0.1), (1.0, 0.2)]

    points = deepcopy(points_XI)
    boundary_nodes = deepcopy(curve_XI)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int)
    @test boundary_curves == DT.to_boundary_curves(points, boundary_nodes) &&
          new_boundary_nodes == [[[1, 2, 3], [3, 1]], [[13, 13]], [[4, 5, 6, 7, 4]], [[14, 8], [8, 14]], [[12, 11, 10, 12]], [[15, 15]]] &&
          points == [(-2.0, 0.0), (0.0, 0.0), (2.0, 0.0), (-2.0, -5.0), (2.0, -5.0), (2.0, -1 / 10), (-2.0, -1 / 10), (-1.0, -3.0), (0.0, -4.0), (0.0, -2.3), (-0.5, -3.5), (0.9, -3.0), (0.0, 0.4), (0.0, -2.0), (1.1, -3.0)]
end

@testset "BoundaryEnricher" begin
    all_points = deepcopy.((points_I, points_II, points_III, points_IV, points_V, points_VI, points_VII, points_VIII, points_IX, points_X, points_XI))
    all_boundary_nodes = deepcopy.((curve_I, curve_II, curve_III, curve_IV, curve_V, curve_VI, curve_VII, curve_VIII, curve_IX, curve_X, curve_XI))
    for (points, boundary_nodes) in zip(all_points, all_boundary_nodes)
        enricher = DT.BoundaryEnricher(points, boundary_nodes; IntegerType=Int)
        @test get_points(enricher) == enricher.points == points
        @test get_boundary_nodes(enricher) == enricher.boundary_nodes
        @test DT.get_boundary_curves(enricher) == enricher.boundary_curves
        _boundary_curves, _boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int)
        @test DT.get_boundary_curves(enricher) == _boundary_curves
        hierarchy = DT.get_polygon_hierarchy(enricher)
        _hierarchy = DT.construct_polygon_hierarchy(points, _boundary_nodes, _boundary_curves)
        DT.expand_bounds!(_hierarchy, DT.ε(Float64))
        @test compare_trees(hierarchy, _hierarchy)
        for i in eachindex(_boundary_curves)
            @test DT.get_boundary_curve(enricher, i) == _boundary_curves[i]
        end
        @test DT.each_boundary_edge(enricher) == keys(enricher.parent_map)
        @test maximum_total_variation(get_points(enricher), get_boundary_nodes(enricher), DT.get_boundary_curves(enricher)) < π / 2
        @test DT.get_queue(enricher) ⊢ DT.Queue{Int}()
        @test DT.get_boundary_edge_map(enricher) == enricher.boundary_edge_map == DT.construct_boundary_edge_map(DT.get_boundary_nodes(enricher))
        @test !DT.has_segments(enricher)
        @test DT.get_segments(enricher) == Set{NTuple{2,Int}}()
        @test !DT.is_segment(enricher, 1, 2)
    end
    points, boundary_nodes, segments = deepcopy(points_I_extra_segments), deepcopy(curve_I), deepcopy(segments_I)
    enricher = DT.BoundaryEnricher(points, boundary_nodes, segments)
    @test DT.get_parent(enricher, 1, 2) == DT.get_parent(enricher, 2, 3) == DT.get_parent(enricher, 3, 4) == DT.get_parent(enricher, 4, 1) == 1
    curve_index_map = DT.get_curve_index_map(enricher)
    @test curve_index_map == Dict(1 => 1)
    @test DT.get_boundary_edge_map(enricher, 2, 3) == (get_boundary_nodes(enricher), 2)
    @test DT.has_segments(enricher)
    @test DT.is_segment(enricher, 55, 56) && DT.is_segment(enricher, 56, 55)
    @test !DT.is_segment(enricher, 51, 50)
    @test DT.get_segments(enricher) === segments
    points, boundary_nodes, segments = deepcopy(points_II_extra_segments), deepcopy(curve_II), deepcopy(segments_II)
    enricher = DT.BoundaryEnricher(points, boundary_nodes, segments)
    @test DT.has_segments(enricher)
    @test DT.get_parent_map(enricher) == Dict(
        (1, 2) => 1,
        (2, 3) => 1,
        (3, 4) => 1,
        (4, 5) => 1,
        (5, 6) => 2,
        (6, 7) => 2,
        (7, 8) => 2,
        (8, 9) => 2,
        (9, 10) => 3,
        (10, 11) => 3,
        (11, 1) => 3,
    )
    @test length(DT.get_parent_map(enricher)) == 11
    curve_index_map = DT.get_curve_index_map(enricher)
    @test curve_index_map == Dict(1 => 1, 2 => 1, 3 => 1)
    @test DT.map_curve_index(enricher, 3) == 1
    rects, els = get_dt_rectangles(enricher.spatial_tree.tree)
    all_edges = Set(DT.get_edge(box) for box in els)
    @test all(e -> e ∈ all_edges || reverse(e) ∈ all_edges, segments)
    points, boundary_nodes = deepcopy(points_XI), deepcopy(curve_XI)
    enricher = DT.BoundaryEnricher(points, boundary_nodes)
    @test DT.get_parent(enricher, 1, 2) == 1
    @test DT.get_parent(enricher, 17, 16) == 2
    @test DT.get_parent(enricher, 21, 13) == 3
    @test DT.get_parent(enricher, 4, 5) == 4
    @test DT.get_parent(enricher, 22, 24) == 5
    @test DT.get_parent(enricher, 26, 25) == 6
    @test DT.get_parent(enricher, 11, 10) == 7
    @test DT.get_parent(enricher, 29, 32) == 8
    DT.update_parent_map!(enricher, 29, 32, 60)
    @test !haskey(DT.get_parent_map(enricher), (29, 32))
    @test DT.get_parent(enricher, 29, 60) == 8
    @test DT.get_parent(enricher, 60, 32) == 8
    curve_index_map = DT.get_curve_index_map(enricher)
    @test curve_index_map == Dict(
        1 => 1,
        2 => 1,
        3 => 2,
        4 => 3,
        5 => 4,
        6 => 4,
        7 => 5,
        8 => 6,
    )
    @test DT.map_curve_index(enricher, 1) == 1
    @test DT.map_curve_index(enricher, 7) == 5
    @test DT.get_boundary_edge_map(enricher, 14, 23) == ((4, 1), 1)
    @test DT.get_boundary_edge_map(enricher, 11, 10) == ((5, 1), 2)
    @test DT.reorient_edge(enricher, 23, 14) == (14, 23)
end

@testset "coarse_discretisation" begin
    # I 
    points = deepcopy(points_I)
    boundary_nodes = deepcopy(curve_I)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int)
    DT.coarse_discretisation!(points, new_boundary_nodes, boundary_curves)
    enricher = DT.BoundaryEnricher(deepcopy(points_I), deepcopy(curve_I))
    @test boundary_curves == DT.get_boundary_curves(enricher) && new_boundary_nodes == get_boundary_nodes(enricher) && get_points(enricher) == points_I
    @test maximum_total_variation(get_points(enricher), get_boundary_nodes(enricher), DT.get_boundary_curves(enricher)) < π / 2

    # II 
    points = deepcopy(points_II)
    boundary_nodes = deepcopy(curve_II)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int)
    DT.coarse_discretisation!(points, new_boundary_nodes, boundary_curves)
    enricher = DT.BoundaryEnricher(deepcopy(points_II), deepcopy(curve_II))
    @test boundary_curves == DT.get_boundary_curves(enricher) && new_boundary_nodes == get_boundary_nodes(enricher) && get_points(enricher) == points_II
    @test maximum_total_variation(get_points(enricher), get_boundary_nodes(enricher), DT.get_boundary_curves(enricher)) < π / 2

    # III
    points = deepcopy(points_III)
    boundary_nodes = deepcopy(curve_III)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int)
    DT.coarse_discretisation!(points, new_boundary_nodes, boundary_curves)
    enricher = DT.BoundaryEnricher(deepcopy(points_III), deepcopy(curve_III))
    @test boundary_curves == DT.get_boundary_curves(enricher) && new_boundary_nodes == get_boundary_nodes(enricher) && get_points(enricher) == points_III
    @test maximum_total_variation(get_points(enricher), get_boundary_nodes(enricher), DT.get_boundary_curves(enricher)) < π / 2

    # IV
    points = deepcopy(points_IV)
    boundary_nodes = deepcopy(curve_IV)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int)
    DT.coarse_discretisation!(points, new_boundary_nodes, boundary_curves)
    @test length(points) == 8
    @test maximum_total_variation(points, new_boundary_nodes, boundary_curves) < π / 2
    all_t = DT.get_inverse.(Ref(curve_IV[1]), get_point(points, new_boundary_nodes...)) |> collect
    all_t[end] = 1.0
    @test issorted(all_t)
    boundary_nodes = new_boundary_nodes
    @test DT.num_boundary_edges(boundary_nodes) == 8
    for i in 1:num_boundary_edges(boundary_nodes)
        j = i + 1
        k = j + 1 > length(boundary_nodes) ? 2 : j + 1
        u, v, w = get_boundary_nodes(boundary_nodes, i), get_boundary_nodes(boundary_nodes, j), get_boundary_nodes(boundary_nodes, k)
        p, q, r = get_point(points, u, v, w)
        t₁, t₂, t₃ = DT.get_inverse(curve_IV[1], p), DT.get_inverse(curve_IV[1], q), DT.get_inverse(curve_IV[1], r)
        if t₃ == 0.0
            t₃ = 1.0 # periodic
        end
        if t₂ == 0.0
            t₂ = 1.0
            t₃ += 1 # periodic
        end
        Tθ₁ = DT.total_variation(curve_IV[1], t₁, t₂)
        Tθ₂ = DT.total_variation(curve_IV[1], t₂, t₃)
        @test Tθ₁ ≈ Tθ₂
    end

    # V 
    points = deepcopy(points_V)
    boundary_nodes = deepcopy(curve_V)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int)
    DT.coarse_discretisation!(points, new_boundary_nodes, boundary_curves; n=5) # 5 gets mapped to 8
    @test length(points) == 9
    all_t = DT.get_inverse.(Ref(curve_V[1]), get_point(points, new_boundary_nodes...)) |> collect
    all_t[end] = 1.0
    @test issorted(all_t)
    boundary_nodes = new_boundary_nodes
    @test DT.num_boundary_edges(boundary_nodes) == 8
    for i in 1:(num_boundary_edges(boundary_nodes)-1)
        j = i + 1
        k = j + 1
        u, v, w = get_boundary_nodes(boundary_nodes, i), get_boundary_nodes(boundary_nodes, j), get_boundary_nodes(boundary_nodes, k)
        p, q, r = get_point(points, u, v, w)
        t₁, t₂, t₃ = DT.get_inverse(curve_V[1], p), DT.get_inverse(curve_V[1], q), DT.get_inverse(curve_V[1], r)
        if t₃ == 0.0
            t₃ = 1.0 # periodic
        end
        if t₂ == 0.0
            t₂ = 0.0
            t₃ += 1 # periodic
        end
        Tθ₁ = slow_total_absolute_curvature(curve_V[1], t₁, t₂)
        Tθ₂ = slow_total_absolute_curvature(curve_V[1], t₂, t₃)
        @test Tθ₁ ≈ Tθ₂ rtol = 1.0e-2
    end

    # VI 
    points = deepcopy(points_VI)
    boundary_nodes = deepcopy(curve_VI)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int)
    DT.coarse_discretisation!(points, new_boundary_nodes, boundary_curves; n=64)
    @test length(points) == 137
    all_t = DT.get_inverse.(Ref(curve_VI[1][1]), get_point(points, get_boundary_nodes(new_boundary_nodes, 1)...)) |> collect
    @test issorted(all_t)
    all_t = DT.get_inverse.(Ref(curve_VI[2][1]), get_point(points, get_boundary_nodes(new_boundary_nodes, 2)...)) |> collect
    @test issorted(all_t)
    @test get_boundary_nodes(new_boundary_nodes, 3) == [5, 6, 10]
    boundary_nodes = new_boundary_nodes
    @test DT.num_boundary_edges(boundary_nodes[1]) == 64
    @test DT.num_boundary_edges(boundary_nodes[2]) == 64
    for idx in 1:2
        section_nodes = get_boundary_nodes(boundary_nodes, idx)
        for i in 1:(num_boundary_edges(section_nodes)-1)
            j = i + 1
            k = j + 1
            u, v, w = get_boundary_nodes(section_nodes, i), get_boundary_nodes(section_nodes, j), get_boundary_nodes(section_nodes, k)
            p, q, r = get_point(points, u, v, w)
            t₁, t₂, t₃ = DT.get_inverse(curve_VI[idx][1], p), DT.get_inverse(curve_VI[idx][1], q), DT.get_inverse(curve_VI[idx][1], r)
            Tθ₁ = DT.total_variation(curve_VI[idx][1], t₁, t₂)
            Tθ₂ = DT.total_variation(curve_VI[idx][1], t₂, t₃)
            @test Tθ₁ ≈ Tθ₂ rtol = 1.0e-2
        end
    end

    # VII 
    points = deepcopy(points_VII)
    boundary_nodes = deepcopy(curve_VII)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int)
    DT.coarse_discretisation!(points, new_boundary_nodes, boundary_curves; n=512)
    @test length(points) == 1025
    all_t = DT.get_inverse.(Ref(curve_VII[1][1]), get_point(points, get_boundary_nodes(new_boundary_nodes, 1)...)) |> collect
    @test issorted(all_t)
    all_t = DT.get_inverse.(Ref(curve_VII[2][1]), get_point(points, get_boundary_nodes(new_boundary_nodes, 2)...)) |> collect
    @test issorted(all_t)
    boundary_nodes = new_boundary_nodes
    @test DT.num_boundary_edges(boundary_nodes[1]) == 512
    @test DT.num_boundary_edges(boundary_nodes[2]) == 512
    for idx in 1:2
        section_nodes = get_boundary_nodes(boundary_nodes, idx)
        for i in 1:(num_boundary_edges(section_nodes)-1)
            j = i + 1
            k = j + 1
            u, v, w = get_boundary_nodes(section_nodes, i), get_boundary_nodes(section_nodes, j), get_boundary_nodes(section_nodes, k)
            p, q, r = get_point(points, u, v, w)
            t₁, t₂, t₃ = DT.get_inverse(curve_VII[idx][1], p), DT.get_inverse(curve_VII[idx][1], q), DT.get_inverse(curve_VII[idx][1], r)
            Tθ₁ = DT.total_variation(curve_VII[idx][1], t₁, t₂)
            Tθ₂ = DT.total_variation(curve_VII[idx][1], t₂, t₃)
            @test Tθ₁ ≈ Tθ₂ rtol = 1.0e-1
        end
    end

    # VIII 
    points = deepcopy(points_VIII)
    boundary_nodes = deepcopy(curve_VIII)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int)
    DT.coarse_discretisation!(points, new_boundary_nodes, boundary_curves) # default is 0
    @test length(points) == 18
    @test maximum_total_variation(points, new_boundary_nodes, boundary_curves) < π / 2
    # ignore 1 and 3 since those don't get enriched 
    all_t = DT.get_inverse.(Ref(curve_VIII[2][1]), get_point(points, get_boundary_nodes(new_boundary_nodes, 2)...)) |> collect
    @test issorted(all_t)
    all_t = DT.get_inverse.(Ref(curve_VIII[4][1]), get_point(points, get_boundary_nodes(new_boundary_nodes, 4)...)) |> collect
    @test issorted(all_t)
    boundary_nodes = new_boundary_nodes
    @test DT.num_boundary_edges(boundary_nodes[2]) == 4
    @test DT.num_boundary_edges(boundary_nodes[4]) == 4
    @test boundary_nodes[1] == [1, 2, 3, 4, 5]
    @test boundary_nodes[3] == [6, 7, 8, 9, 10]
    for idx in [2, 4]
        section_nodes = get_boundary_nodes(boundary_nodes, idx)
        for i in 1:(num_boundary_edges(section_nodes)-1)
            j = i + 1
            k = j + 1
            u, v, w = get_boundary_nodes(section_nodes, i), get_boundary_nodes(section_nodes, j), get_boundary_nodes(section_nodes, k)
            p, q, r = get_point(points, u, v, w)
            t₁, t₂, t₃ = DT.get_inverse(curve_VIII[idx][1], p), DT.get_inverse(curve_VIII[idx][1], q), DT.get_inverse(curve_VIII[idx][1], r)
            Tθ₁ = DT.total_variation(curve_VIII[idx][1], t₁, t₂)
            Tθ₂ = DT.total_variation(curve_VIII[idx][1], t₂, t₃)
            @test Tθ₁ ≈ Tθ₂ rtol = 1.0e-2
        end
    end

    # IX 
    points = deepcopy(points_IX)
    boundary_nodes = deepcopy(curve_IX)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int)
    DT.coarse_discretisation!(points, new_boundary_nodes, boundary_curves, n=1) # check that 1 becomes 4
    @test length(points) == 11
    # 1 gets ignored because it's a piecewise linear curve 
    all_t = DT.get_inverse.(Ref(curve_IX[2][1][1]), get_point(points, new_boundary_nodes[2][1]...)) |> collect
    all_t[end] = 1.0 # periodic
    @test issorted(all_t)
    boundary_nodes = new_boundary_nodes
    @test DT.num_boundary_edges(boundary_nodes[2][1]) == 4
    @test boundary_nodes[1][1] == [1, 2, 3, 4, 5, 6, 7, 1]
    for idx in (2,)
        curve_nodes = get_boundary_nodes(boundary_nodes, idx)
        section_nodes = get_boundary_nodes(curve_nodes, 1)
        for i in 1:(num_boundary_edges(section_nodes)-1)
            j = i + 1
            k = j + 1
            u, v, w = get_boundary_nodes(section_nodes, i), get_boundary_nodes(section_nodes, j), get_boundary_nodes(section_nodes, k)
            p, q, r = get_point(points, u, v, w)
            t₁, t₂, t₃ = DT.get_inverse(curve_IX[idx][1][1], p), DT.get_inverse(curve_IX[idx][1][1], q), DT.get_inverse(curve_IX[idx][1][1], r)
            if t₃ == 0.0
                t₃ = 1.0 # periodic
            end
            if t₂ == 0.0
                t₂ = 1.0
                t₃ += 1 # periodic
            end
            Tθ₁ = DT.total_variation(curve_IX[idx][1][1], t₁, t₂)
            Tθ₂ = DT.total_variation(curve_IX[idx][1][1], t₂, t₃)
            @test Tθ₁ ≈ Tθ₂
        end
    end

    # X
    points = deepcopy(points_X)
    boundary_nodes = deepcopy(curve_X)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int)
    DT.coarse_discretisation!(points, new_boundary_nodes, boundary_curves; n=32)
    @test length(points) == 70
    # 1 and 4 are piecewise linear curves 
    all_t = DT.get_inverse.(Ref(curve_X[1][2][1]), get_point(points, new_boundary_nodes[1][2]...)) |> collect
    @test issorted(all_t)
    all_t = DT.get_inverse.(Ref(curve_X[2][1][1]), get_point(points, new_boundary_nodes[2][1]...)) |> collect
    @test issorted(all_t)
    boundary_nodes = new_boundary_nodes
    @test DT.num_boundary_edges(boundary_nodes[1][2]) == 32
    @test DT.num_boundary_edges(boundary_nodes[2][1]) == 32
    @test boundary_nodes[1][1] == [1, 2, 3]
    @test boundary_nodes[2][2] == [8, 7, 6, 5, 4]
    for (idx, idx2) in zip((1, 2), (2, 1))
        curve_nodes = get_boundary_nodes(boundary_nodes, idx)
        section_nodes = get_boundary_nodes(curve_nodes, idx2)
        for i in 1:(num_boundary_edges(section_nodes)-1)
            j = i + 1
            k = j + 1
            u, v, w = get_boundary_nodes(section_nodes, i), get_boundary_nodes(section_nodes, j), get_boundary_nodes(section_nodes, k)
            p, q, r = get_point(points, u, v, w)
            t₁, t₂, t₃ = DT.get_inverse(curve_X[idx][idx2][1], p), DT.get_inverse(curve_X[idx][idx2][1], q), DT.get_inverse(curve_X[idx][idx2][1], r)
            Tθ₁ = DT.total_variation(curve_X[idx][idx2][1], t₁, t₂)
            Tθ₂ = DT.total_variation(curve_X[idx][idx2][1], t₂, t₃)
            @test Tθ₁ ≈ Tθ₂ rtol = 1.0e-1
        end
    end

    # XI
    points = deepcopy(points_XI)
    boundary_nodes = deepcopy(curve_XI)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, Int32)
    DT.coarse_discretisation!(points, new_boundary_nodes, boundary_curves; n=128)
    @test length(points) == 650
    # (1, 1), (3, 1), and (5, 1) are piecewise linear 
    # (2, 1) and (6, 1) are periodic 
    all_t = DT.get_inverse.(Ref(curve_XI[1][2][1]), get_point(points, new_boundary_nodes[1][2]...)) |> collect
    @test issorted(all_t)
    all_t = DT.get_inverse.(Ref(curve_XI[2][1][1]), get_point(points, new_boundary_nodes[2][1]...)) |> collect
    all_t[end] = 1.0
    @test issorted(all_t)
    all_t = DT.get_inverse.(Ref(curve_XI[4][1][1]), get_point(points, new_boundary_nodes[4][1]...)) |> collect
    @test issorted(all_t)
    all_t = DT.get_inverse.(Ref(curve_XI[6][1][1]), get_point(points, new_boundary_nodes[6][1]...)) |> collect
    all_t[end] = 1.0
    @test issorted(all_t)
    boundary_nodes = new_boundary_nodes
    @test DT.num_boundary_edges(boundary_nodes[1][2]) == 128
    @test DT.num_boundary_edges(boundary_nodes[2][1]) == 128
    @test DT.num_boundary_edges(boundary_nodes[4][1]) == 128
    @test DT.num_boundary_edges(boundary_nodes[6][1]) == 128
    @test boundary_nodes[1][1] == [1, 2, 3]
    @test boundary_nodes[3][1] == [4, 5, 6, 7, 4]
    @test boundary_nodes[5][1] == [12, 11, 10, 12]
    for (idx, idx2) in zip((1, 2, 4, 4, 6), (2, 1, 1, 2, 1))
        curve_nodes = get_boundary_nodes(boundary_nodes, idx)
        section_nodes = get_boundary_nodes(curve_nodes, idx2)
        for i in 1:(num_boundary_edges(section_nodes)-1)
            j = i + 1
            k = j + 1
            u, v, w = get_boundary_nodes(section_nodes, i), get_boundary_nodes(section_nodes, j), get_boundary_nodes(section_nodes, k)
            p, q, r = get_point(points, u, v, w)
            t₁, t₂, t₃ = DT.get_inverse(curve_XI[idx][idx2][1], p), DT.get_inverse(curve_XI[idx][idx2][1], q), DT.get_inverse(curve_XI[idx][idx2][1], r)
            if t₃ == 0.0
                t₃ = 1.0 # periodic
            end
            if t₂ == 0.0
                t₂ = 1.0
                t₃ += 1 # periodic
            end
            Tθ₁ = DT.total_variation(curve_XI[idx][idx2][1], t₁, t₂)
            Tθ₂ = DT.total_variation(curve_XI[idx][idx2][1], t₂, t₃)
            @test Tθ₁ ≈ Tθ₂ rtol = 1.0e-2
        end
    end
end

@testset "Calling on BoundaryEnricher and performance checking" begin
    enricher_I = DT.BoundaryEnricher(deepcopy(points_I), deepcopy(curve_I))
    enricher_II = DT.BoundaryEnricher(deepcopy(points_II), deepcopy(curve_II))
    enricher_III = DT.BoundaryEnricher(deepcopy(points_III), deepcopy(curve_III))
    enricher_IV = DT.BoundaryEnricher(deepcopy(points_IV), deepcopy(curve_IV))
    enricher_V = DT.BoundaryEnricher(deepcopy(points_V), deepcopy(curve_V))
    enricher_VI = DT.BoundaryEnricher(deepcopy(points_VI), deepcopy(curve_VI))
    enricher_VII = DT.BoundaryEnricher(deepcopy(points_VII), deepcopy(curve_VII))
    enricher_VIII = DT.BoundaryEnricher(deepcopy(points_VIII), deepcopy(curve_VIII))
    enricher_IX = DT.BoundaryEnricher(deepcopy(points_IX), deepcopy(curve_IX))
    enricher_X = DT.BoundaryEnricher(deepcopy(points_X), deepcopy(curve_X))
    enricher_XI = DT.BoundaryEnricher(deepcopy(points_XI), deepcopy(curve_XI))

    @test DT.is_piecewise_linear(enricher_I, 1)
    @test DT.is_piecewise_linear(enricher_II, 1)
    @test DT.is_piecewise_linear(enricher_III, 1)
    @test !DT.is_piecewise_linear(enricher_IV, 1)
    @test DT.is_piecewise_linear(enricher_XI, 1) && !DT.is_piecewise_linear(enricher_XI, 2) &&
          !DT.is_piecewise_linear(enricher_XI, 3) && DT.is_piecewise_linear(enricher_XI, 4) && !DT.is_piecewise_linear(enricher_XI, 5) &&
          !DT.is_piecewise_linear(enricher_XI, 6) && DT.is_piecewise_linear(enricher_XI, 7) && !DT.is_piecewise_linear(enricher_XI, 8)
    @inferred DT.is_piecewise_linear(enricher_I, 1)
    @inferred DT.is_piecewise_linear(enricher_II, 3)
    @inferred DT.is_piecewise_linear(enricher_III, 3)
    @inferred DT.is_piecewise_linear(enricher_IV, 4)
    @inferred DT.is_piecewise_linear(enricher_V, 5)
    @inferred DT.is_piecewise_linear(enricher_VI, 6)
    @inferred DT.is_piecewise_linear(enricher_VII, 10)
    @inferred DT.is_piecewise_linear(enricher_VIII, 1)
    @inferred DT.is_piecewise_linear(enricher_IX, 1)
    @inferred DT.is_piecewise_linear(enricher_X, 1)
    @inferred DT.is_piecewise_linear(enricher_XI, 1)

    curve = enricher_IV.boundary_curves[1]
    p = curve(0.5)
    @test DT.get_inverse(enricher_IV, 1, p) ≈ 0.5
    @inferred DT.get_inverse(enricher_IV, 1, p)
    curve = enricher_XI.boundary_curves[2]
    p = curve(0.381)
    @test DT.get_inverse(enricher_XI, 2, p) ≈ 0.381
    @inferred DT.get_inverse(enricher_XI, 2, p)
    curve = enricher_XI.boundary_curves[5]
    p = curve(0.998)
    @test DT.get_inverse(enricher_XI, 5, p) ≈ 0.998
    @inferred DT.get_inverse(enricher_XI, 5, p)
    curve = enricher_X.boundary_curves[3]
    p = curve(0.111)
    @test DT.get_inverse(enricher_X, 3, p) ≈ 0.111
    @inferred DT.get_inverse(enricher_X, 3, p)

    curve = enricher_IV.boundary_curves[1]
    t₁, t₂ = 0.2, 0.8
    @test DT.get_equivariation_split(enricher_IV, 1, t₁, t₂) == DT.get_equivariation_split(curve, t₁, t₂)
    @inferred DT.get_equivariation_split(enricher_IV, 1, t₁, t₂)
    curve = enricher_XI.boundary_curves[2]
    t₁, t₂ = 0.2, 0.8
    @test DT.get_equivariation_split(enricher_XI, 2, t₁, t₂) == DT.get_equivariation_split(curve, t₁, t₂)
    @inferred DT.get_equivariation_split(enricher_XI, 2, t₁, t₂)
    curve = enricher_XI.boundary_curves[3]
    t₁, t₂ = 0.2, 0.8
    @test DT.get_equivariation_split(enricher_XI, 3, t₁, t₂) == DT.get_equivariation_split(curve, t₁, t₂)
    @inferred DT.get_equivariation_split(enricher_XI, 3, t₁, t₂)
    curve = enricher_X.boundary_curves[3]
    t₁, t₂ = 0.2, 0.55
    @test DT.get_equivariation_split(enricher_X, 3, t₁, t₂) == DT.get_equivariation_split(curve, t₁, t₂)
    @inferred DT.get_equivariation_split(enricher_X, 3, t₁, t₂)

    curve = enricher_IV.boundary_curves[1]
    t = rand()
    @test DT.eval_boundary_curve(enricher_IV, 1, t) == curve(t)
    @inferred DT.eval_boundary_curve(enricher_IV, 1, t)
    curve = enricher_XI.boundary_curves[2]
    t = rand()
    @test DT.eval_boundary_curve(enricher_XI, 2, t) == curve(t)
    @inferred DT.eval_boundary_curve(enricher_XI, 2, t)
    curve = enricher_XI.boundary_curves[8]
    t = rand()
    @test DT.eval_boundary_curve(enricher_XI, 8, t) == curve(t)
    @inferred DT.eval_boundary_curve(enricher_XI, 8, t)
    curve = enricher_X.boundary_curves[3]
    t = rand()
    @test DT.eval_boundary_curve(enricher_X, 3, t) == curve(t)
    @inferred DT.eval_boundary_curve(enricher_X, 3, t)

    curve = enricher_IV.boundary_curves[1]
    p = (10randn(), 10randn())
    @test DT.point_position_relative_to_curve(rt(), enricher_IV, 1, p) == DT.point_position_relative_to_curve(rt(), curve, p)
    @test DT.point_position_relative_to_curve(enricher_IV, 1, p) == DT.point_position_relative_to_curve(AdaptiveKernel(), curve, p)
    @inferred DT.point_position_relative_to_curve(rt(), enricher_IV, 1, p)
    curve = enricher_XI.boundary_curves[2]
    p = (10randn(), 10randn())
    @test DT.point_position_relative_to_curve(rt(), enricher_XI, 2, p) == DT.point_position_relative_to_curve(rt(), curve, p)
    @inferred DT.point_position_relative_to_curve(rt(), enricher_XI, 2, p)
    curve = enricher_XI.boundary_curves[8]
    p = (10randn(), 10randn())
    @test DT.point_position_relative_to_curve(rt(), enricher_XI, 8, p) == DT.point_position_relative_to_curve(rt(), curve, p)
    @inferred DT.point_position_relative_to_curve(rt(), enricher_XI, 8, p)
    curve = enricher_X.boundary_curves[3]
    p = (10randn(), 10randn())
    @test DT.point_position_relative_to_curve(rt(), enricher_X, 3, p) == DT.point_position_relative_to_curve(rt(), curve, p)
    @inferred DT.point_position_relative_to_curve(rt(), enricher_X, 3, p)
    @test DT.point_position_relative_to_curve(DT.get_boundary_curves(enricher_X), 3, p) == DT.point_position_relative_to_curve(AdaptiveKernel(), DT.get_boundary_curves(enricher_X), 3, p)

    curve = enricher_IV.boundary_curves[1]
    t₁, t₂, r = rand(3)
    @test DT.get_circle_intersection(enricher_IV, 1, t₁, t₂, r) == DT.get_circle_intersection(curve, t₁, t₂, r)
    @inferred DT.get_circle_intersection(enricher_IV, 1, t₁, t₂, r)
    curve = enricher_XI.boundary_curves[2]
    t₁, t₂, r = rand(3)
    @test DT.get_circle_intersection(enricher_XI, 2, t₁, t₂, r) == DT.get_circle_intersection(curve, t₁, t₂, r)
    @inferred DT.get_circle_intersection(enricher_XI, 2, t₁, t₂, r)
    curve = enricher_XI.boundary_curves[8]
    t₁, t₂, r = rand(3)
    @test DT.get_circle_intersection(enricher_XI, 8, t₁, t₂, r) == DT.get_circle_intersection(curve, t₁, t₂, r)
    @inferred DT.get_circle_intersection(enricher_XI, 8, t₁, t₂, r)
    curve = enricher_X.boundary_curves[3]
    t₁, t₂, r = rand(3)
    @test DT.get_circle_intersection(enricher_X, 3, t₁, t₂, r) == DT.get_circle_intersection(curve, t₁, t₂, r)
    @inferred DT.get_circle_intersection(enricher_X, 3, t₁, t₂, r)

    curve = enricher_IV.boundary_curves[1]
    t₁, t₂ = rand(2)
    @test DT.get_equidistant_split(enricher_IV, 1, t₁, t₂) == DT.get_equidistant_split(curve, t₁, t₂)
    @inferred DT.get_equidistant_split(enricher_IV, 1, t₁, t₂)
    curve = enricher_XI.boundary_curves[2]
    t₁, t₂ = rand(2)
    @test DT.get_equidistant_split(enricher_XI, 2, t₁, t₂) == DT.get_equidistant_split(curve, t₁, t₂)
    @inferred DT.get_equidistant_split(enricher_XI, 2, t₁, t₂)
    curve = enricher_XI.boundary_curves[8]
    t₁, t₂ = rand(2)
    @test DT.get_equidistant_split(enricher_XI, 8, t₁, t₂) == DT.get_equidistant_split(curve, t₁, t₂)
    @inferred DT.get_equidistant_split(enricher_XI, 8, t₁, t₂)
    curve = enricher_X.boundary_curves[3]
    t₁, t₂ = rand(2)
    @test DT.get_equidistant_split(enricher_X, 3, t₁, t₂) == DT.get_equidistant_split(curve, t₁, t₂)
    @inferred DT.get_equidistant_split(enricher_X, 3, t₁, t₂)
end

@testset "Visibility testing" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        # Simple domain 
        a, b, c, d, e, f, g, h, ii, jj, kk, ℓ, m, n, o, p, r, s, t, u, v, w, z, a1, b1, c1, d1, e1, f1, g1 =
            (0.0, 0.0), (5.0, 0.0), (10.0, 0.0), (10.0, 5.0),
            (10.0, 10.0), (5.0, 10.0),
            (0.0, 10.0), (0.0, 5.0), (2.5, 0.0), (7.5, 0.0),
            (10.0, 2.5), (0.0, 2.5), (0.0, 7.5), (2.5, 10.0),
            (7.5, 10.0), (10.0, 7.5), (3.5, 1.1),
            (2.0, 1.5), (3.0, 0.5), (3.5, 1.5), (4.0, 0.5),
            (5.5, 0.5), (6.5, 2.0), (3.3, 3.0),
            (4.0, 0.8), (2.8, 0.2), (4.4, 0.2), (3.0, 0.6),
            (3.4, -0.2), (3.4, 0.0)
        geo1 = [[g, m, h, ℓ, a, ii, b, jj, c, kk, d, p, e, o, f, n, g]]
        geo2 = [[a1, z, w, v, u, t, s, a1]]
        boundary_nodes, points = convert_boundary_points_to_indices([geo1, geo2]; existing_points=[r, c1, b1, d1, e1, f1, g1])
        enricher = DT.BoundaryEnricher(points, boundary_nodes)
        i, j, k = findfirst(==(ii), points), findfirst(==(b), points), findfirst(==(r), points)
        @test DT.is_invisible(DT.test_visibility(PT(), enricher, i, j, k))
        k = findfirst(==(c1), points)
        @test DT.is_visible(DT.test_visibility(PT(), enricher, i, j, k))
        k = findfirst(==(d1), points)
        @test DT.is_visible(DT.test_visibility(PT(), enricher, i, j, k))
        k = findfirst(==(e1), points)
        @test DT.is_invisible(DT.test_visibility(PT(), enricher, i, j, k))
        k = findfirst(==(b1), points)
        @test DT.is_invisible(DT.test_visibility(PT(), enricher, i, j, k))
        k = findfirst(==(f1), points)
        @test DT.is_invisible(DT.test_visibility(PT(), enricher, i, j, k))
        k = findfirst(==(g1), points)
        @test DT.is_visible(DT.test_visibility(PT(), enricher, i, j, k))
        @inferred DT.test_visibility(PT(), enricher, i, j, k)

        # Curve I 
        points, boundary_nodes = deepcopy(points_I), deepcopy(curve_I)
        enricher = DT.BoundaryEnricher(points, boundary_nodes)
        push!(points, (0.5, 0.2), (0.2, -0.1))
        @test DT.is_visible(DT.test_visibility(PT(), enricher, 1, 2, 5))
        @test DT.is_invisible(DT.test_visibility(PT(), enricher, 1, 2, 6))

        # Curve II 
        points, boundary_nodes = deepcopy(points_II), deepcopy(curve_II)
        enricher = DT.BoundaryEnricher(points, boundary_nodes)
        n = DT.num_points(get_points(enricher))
        push!(points, (0.05, -0.05), (0.6, 0.1), (1.0, 0.1), (0.025, 0.05))
        @test DT.is_invisible(DT.test_visibility(PT(), enricher, 1, 2, n + 1))
        @test DT.is_visible(DT.test_visibility(PT(), enricher, 3, 4, n + 2))
        @test DT.is_visible(DT.test_visibility(PT(), enricher, 5, 6, n + 3))
        @test DT.is_invisible(DT.test_visibility(PT(), enricher, 1, 2, n + 4))

        # Curve III
        points, boundary_nodes = deepcopy(points_III), deepcopy(curve_III)
        enricher = DT.BoundaryEnricher(points, boundary_nodes)
        @test DT.is_visible(DT.test_visibility(PT(), enricher, 1, 2, 11))
        @test DT.is_visible(DT.test_visibility(PT(), enricher, 9, 10, 14))
        @test DT.is_visible(DT.test_visibility(PT(), enricher, 9, 10, 15))

        # Curve IX
        points, boundary_nodes = deepcopy(points_IX), deepcopy(curve_IX)
        enricher = DT.BoundaryEnricher(points, boundary_nodes)
        @test DT.is_invisible(DT.test_visibility(PT(), enricher, 1, 2, 11))
        @test DT.is_visible(DT.test_visibility(PT(), enricher, 1, 2, 13))
        @test DT.is_visible(DT.test_visibility(PT(), enricher, 1, 2, 14))
        @test DT.is_visible(DT.test_visibility(PT(), enricher, 1, 2, 10))
        @test DT.is_visible(DT.test_visibility(PT(), enricher, 2, 3, 10))
        @test DT.is_visible(DT.test_visibility(PT(), enricher, 2, 3, 12))
        @test DT.is_visible(DT.test_visibility(PT(), enricher, 2, 3, 15))
        @test DT.is_visible(DT.test_visibility(PT(), enricher, 2, 3, 14))
        @test DT.is_invisible(DT.test_visibility(PT(), enricher, 2, 3, 9))

        # Curve X
        points, boundary_nodes = deepcopy(points_X), deepcopy(curve_X)
        enricher = DT.BoundaryEnricher(points, boundary_nodes)
        @test DT.is_visible(DT.test_visibility(PT(), enricher, 2, 3, 9))
        @test DT.is_invisible(DT.test_visibility(PT(), enricher, 2, 3, 14))
        @test DT.is_invisible(DT.test_visibility(PT(), enricher, 2, 3, 13))

        # An interior segment 
        points, boundary_nodes, segments = deepcopy(points_I_extra_segments), deepcopy(curve_I), deepcopy(segments_I)
        enricher = DT.BoundaryEnricher(points, boundary_nodes, segments)
        push!(enricher.points, (0.5, 0.8))
        vis = DT.test_visibility(PT(), enricher, 1, 2, length(enricher.points))
        @test DT.is_visible(vis)
        push!(enricher.points, (0.5, 0.0))
        DT.split_edge!(enricher, 1, 2, length(enricher.points))
        vis = DT.test_visibility(PT(), enricher, length(enricher.points), 2, length(enricher.points) - 1)
        @test DT.is_invisible(vis)
        enricher.points[end-1] = (0.5, 0.49)
        vis = DT.test_visibility(PT(), enricher, length(enricher.points), 2, length(enricher.points) - 1)
        @test DT.is_visible(vis)
    end
end

@testset "SmallAngleComplexes" begin
    member = SACM(1, 7)
    nmember = DT.replace_next_edge(member, 12)
    @test DT.get_parent_curve(nmember) == 1 && DT.get_next_edge(nmember) == 12

    points, boundary_nodes = deepcopy(points_I), deepcopy(curve_I)
    enricher = DT.BoundaryEnricher(points, boundary_nodes)
    points, boundary_nodes, boundary_curves = DT.get_points(enricher), DT.get_boundary_nodes(enricher), DT.get_boundary_curves(enricher)
    complexes = DT.get_small_angle_complexes(points, boundary_nodes, boundary_curves)
    @test isempty(complexes)
    @inferred DT.get_small_angle_complexes(points, boundary_nodes, boundary_curves)
    @test DT.get_small_angle_complexes(enricher) == complexes

    points, boundary_nodes = deepcopy(points_II), deepcopy(curve_II)
    enricher = DT.BoundaryEnricher(points, boundary_nodes)
    points, boundary_nodes, boundary_curves = DT.get_points(enricher), DT.get_boundary_nodes(enricher), DT.get_boundary_curves(enricher)
    complexes = DT.get_small_angle_complexes(points, boundary_nodes, boundary_curves)
    _complexes = Dict(
        9 => [SAC(9, [SACM(2, 8), SACM(3, 10)])],
        1 => [SAC(1, [SACM(3, 11), SACM(1, 2)])],
    )
    @test complexes == _complexes
    @test DT.get_small_angle_complexes(enricher) == complexes
    @test DT.is_small_angle_complex_apex(enricher, 9)
    @test DT.is_small_angle_complex_apex(enricher, 1)
    @test !DT.is_small_angle_complex_apex(enricher, 5)
    flag98, apex98, cid98, mid98 = DT.is_small_angle_complex_member(enricher, 9, 8)
    flag910, apex910, cid910, mid910 = DT.is_small_angle_complex_member(enricher, 9, 10)
    flag111, apex111, cid111, mid111 = DT.is_small_angle_complex_member(enricher, 1, 11)
    flag12, apex12, cid12, mid12 = DT.is_small_angle_complex_member(enricher, 1, 2)
    @test all((flag98, flag910, flag111, flag12))
    @test apex98 == apex910 == 9 && apex111 == apex12 == 1
    @test cid98 == cid910 == cid111 == cid12 == 1
    @test mid98 == 1 && mid910 == 2 && mid111 == 1 && mid12 == 2
    flag510, apex510, cid510, mid510 = DT.is_small_angle_complex_member(enricher, 5, 10)
    @test !flag510 && apex510 == cid510 == mid510 == 0
    flag59, apex59, cid59, mid59 = DT.is_small_angle_complex_member(enricher, 5, 9)
    @test !flag59 && apex59 == cid59 == mid59 == 0
    DT.replace_next_edge!(enricher, 9, 1, 2, 17)
    _complexes = Dict(
        9 => [SAC(9, [SACM(2, 8), SACM(3, 17)])],
        1 => [SAC(1, [SACM(3, 11), SACM(1, 2)])],
    )
    @test DT.get_small_angle_complexes(enricher) == _complexes

    A, B, C, D, E, F, G, H, I, J, K = (0.0, 0.0), (0.2, 1.4), (0.6, 1.2),
    (1.2, 0.2), (1.2, -0.2), (-1.4, -0.2),
    (-1.0, -0.6), (0.6, 1.0), (0.8, 0.6),
    (0.6, 0.4), (0.6, 0.2)
    points = [A, B, C, D, E, F, G, H, I, J, K]
    boundary_nodes = [[[1, 3, 2, 1]], [[1, 9, 8, 1]], [[1, 11, 10, 1]], [[1, 5, 4, 1]], [[1, 6, 7, 1]]]
    enricher = DT.BoundaryEnricher(points, boundary_nodes)
    points, boundary_nodes, boundary_curves = get_points(enricher), get_boundary_nodes(enricher), DT.get_boundary_curves(enricher)
    complexes = DT.get_small_angle_complexes(points, boundary_nodes, DT.get_boundary_curves(enricher))
    _complexes = Dict(
        1 => [
            SAC(1, [SACM(1, 2), SACM(1, 3), SACM(2, 8), SACM(2, 9), SACM(3, 10), SACM(3, 11), SACM(4, 4), SACM(4, 5)]),
            SAC(1, [SACM(5, 7), SACM(5, 6)]),
        ],
    )
    @test complexes == _complexes
    @test DT.get_small_angle_complexes(enricher) == complexes
    @test DT.is_small_angle_complex_apex(enricher, 1)
    @test !DT.is_small_angle_complex_apex(enricher, 2)
    len = DT.get_minimum_edge_length(complexes[1][1], points)
    p = get_point(points, 1)
    qs = get_point(points, 2, 3, 8, 9, 10, 11, 4, 5)
    @test len ≈ minimum([norm(p .- q) for q in qs])
    len = DT.get_minimum_edge_length(complexes[1][2], points)
    p = get_point(points, 1)
    qs = get_point(points, 7, 6)
    @test len ≈ minimum([norm(p .- q) for q in qs])
    @test DT.get_small_angle_complexes(enricher, 1) == complexes[1]
    flag12, apex12, cid12, mid12 = DT.is_small_angle_complex_member(enricher, 1, 2)
    flag13, apex13, cid13, mid13 = DT.is_small_angle_complex_member(enricher, 1, 3)
    flag18, apex18, cid18, mid18 = DT.is_small_angle_complex_member(enricher, 1, 8)
    flag19, apex19, cid19, mid19 = DT.is_small_angle_complex_member(enricher, 1, 9)
    flag110, apex110, cid110, mid110 = DT.is_small_angle_complex_member(enricher, 1, 10)
    flag111, apex111, cid111, mid111 = DT.is_small_angle_complex_member(enricher, 1, 11)
    flag14, apex14, cid14, mid14 = DT.is_small_angle_complex_member(enricher, 1, 4)
    flag15, apex15, cid15, mid15 = DT.is_small_angle_complex_member(enricher, 1, 5)
    flag17, apex17, cid17, mid17 = DT.is_small_angle_complex_member(enricher, 1, 7)
    flag16, apex16, cid16, mid16 = DT.is_small_angle_complex_member(enricher, 1, 6)
    @test all((flag12, flag13, flag18, flag19, flag110, flag111, flag14, flag15, flag17, flag16))
    @test apex12 == apex13 == apex18 == apex19 == apex110 == apex111 == apex14 == apex15 == apex17 == apex16 == 1
    @test cid12 == cid13 == cid18 == cid19 == cid110 == cid111 == cid14 == cid15 == 1
    @test cid17 == cid16 == 2
    @test mid12 == 1 && mid13 == 2 && mid18 == 3 && mid19 == 4 && mid110 == 5 && mid111 == 6 && mid14 == 7 && mid15 == 8
    @test mid17 == 1 && mid16 == 2
    flag23, apex23, cid23, mid23 = DT.is_small_angle_complex_member(enricher, 2, 3)
    @test !flag23 && apex23 == cid23 == mid23 == 0
    flag21, apex21, cid21, mid21 = DT.is_small_angle_complex_member(enricher, 2, 1)
    @inferred DT.is_small_angle_complex_member(enricher, 2, 1)
    flag31, apex31, cid31, mid31 = DT.is_small_angle_complex_member(enricher, 3, 1)
    flag81, apex81, cid81, mid81 = DT.is_small_angle_complex_member(enricher, 8, 1)
    flag91, apex91, cid91, mid91 = DT.is_small_angle_complex_member(enricher, 9, 1)
    flag101, apex101, cid101, mid101 = DT.is_small_angle_complex_member(enricher, 10, 1)
    flag111, apex111, cid111, mid111 = DT.is_small_angle_complex_member(enricher, 11, 1)
    flag41, apex41, cid41, mid41 = DT.is_small_angle_complex_member(enricher, 4, 1)
    flag51, apex51, cid51, mid51 = DT.is_small_angle_complex_member(enricher, 5, 1)
    flag71, apex71, cid71, mid71 = DT.is_small_angle_complex_member(enricher, 7, 1)
    flag61, apex61, cid61, mid61 = DT.is_small_angle_complex_member(enricher, 6, 1)
    @test all((flag21, flag31, flag81, flag91, flag101, flag111, flag41, flag51, flag71, flag61))
    @test apex21 == apex31 == apex81 == apex91 == apex101 == apex111 == apex41 == apex51 == apex71 == apex61 == 1
    @test cid21 == cid31 == cid81 == cid91 == cid101 == cid111 == cid41 == cid51 == 1
    @test cid71 == cid61 == 2
    @test mid21 == 1 && mid31 == 2 && mid81 == 3 && mid91 == 4 && mid101 == 5 && mid111 == 6 && mid41 == 7 && mid51 == 8
    @test mid71 == 1 && mid61 == 2
    DT.replace_next_edge!(enricher, 1, 2, 2, 20)
    _complexes = Dict(
        1 => [
            SAC(1, [SACM(1, 2), SACM(1, 3), SACM(2, 8), SACM(2, 9), SACM(3, 10), SACM(3, 11), SACM(4, 4), SACM(4, 5)]),
            SAC(1, [SACM(5, 7), SACM(5, 20)]),
        ],
    )
    @test DT.get_small_angle_complexes(enricher) == _complexes

    points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    boundary_nodes = [[1, 2], [2, 3], [3, 4], [4, 1]]
    enricher = DT.BoundaryEnricher(points, boundary_nodes)
    points, boundary_nodes, boundary_curves = get_points(enricher), get_boundary_nodes(enricher), DT.get_boundary_curves(enricher)
    complexes = DT.get_small_angle_complexes(points, boundary_nodes, DT.get_boundary_curves(enricher))
    @test isempty(complexes)
    @test DT.get_small_angle_complexes(enricher) == complexes

    points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    boundary_nodes = [[[1, 2], [2, 3], [3, 4], [4, 1]]]
    enricher = DT.BoundaryEnricher(points, boundary_nodes)
    points, boundary_nodes, boundary_curves = get_points(enricher), get_boundary_nodes(enricher), DT.get_boundary_curves(enricher)
    complexes = DT.get_small_angle_complexes(points, boundary_nodes, DT.get_boundary_curves(enricher))
    @test isempty(complexes)
    @test DT.get_small_angle_complexes(enricher) == complexes

    curves = [curve_III, curve_IV, curve_V, curve_VI, curve_VII, curve_VIII, curve_IX, curve_X, curve_XI]
    points = [points_III, points_IV, points_V, points_VI, points_VII, points_VIII, points_IX, points_X, points_XI]

    for i in eachindex(curves, points)
        points_i, boundary_nodes_i = deepcopy(points[i]), deepcopy(curves[i])
        enricher = DT.BoundaryEnricher(points_i, boundary_nodes_i)
        _points, boundary_nodes, boundary_curves = DT.get_points(enricher), DT.get_boundary_nodes(enricher), DT.get_boundary_curves(enricher)
        complexes = DT.get_small_angle_complexes(_points, boundary_nodes, boundary_curves)
        @test isempty(complexes)
        @test DT.get_small_angle_complexes(enricher) == complexes
    end
end

@testset "construct_segment_map" begin
    segments = segments_II
    points = points_II_extra_segments
    segment_map = DT.construct_segment_map(segments, points, Int)
    @test segment_map == Dict(
        16 => [18, 17, 19, 22, 21, 20],
        17 => [16],
        18 => [16],
        19 => [16],
        20 => [16],
        21 => [16],
        22 => [16],
    )
end

@testset "SmallAngleComplex with segments (from curve_II)" begin
    points, boundary_nodes, segments = deepcopy(points_II_extra_segments), deepcopy(curve_II), deepcopy(segments_II)
    enricher = DT.BoundaryEnricher(points, boundary_nodes, segments)
    complexes = DT.get_small_angle_complexes(enricher)
    _complexes = Dict(
        9 => [SAC(9, [SACM(2, 8), SACM(3, 10)])],
        1 => [SAC(1, [SACM(3, 11), SACM(1, 2)])],
        16 => [SAC(16, [SACM(0, 18), SACM(0, 17), SACM(0, 19), SACM(0, 22)]), SAC(16, [SACM(0, 21), SACM(0, 20)])],
    )
    @test complexes == _complexes

    points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.5, 0.5), (0.5, 0.9), (0.5, 0.1)]
    boundary_nodes = [1, 2, 3, 4, 1]
    segments = Set([(5, 6), (5, 7)])
    enricher = DT.BoundaryEnricher(points, boundary_nodes, segments)
    complexes = DT.get_small_angle_complexes(enricher)
    @test isempty(complexes)
end

@testset "SmallAngleComplex ==" begin
    member1 = SACM(1, 2)
    member2 = SACM(1, 3)
    @test member1 ≠ member2
    member2 = SACM(1, 2)
    @test member1 == member2
    member2 = SACM(2, 2)
    @test member1 ≠ member2

    member1 = SACM(1, 2)
    member2 = SACM(2, 4)
    member3 = SACM(3, 5)
    apex = 1
    complex1 = SAC(apex, [member1, member2, member3])
    complex2 = SAC(apex, [member1, member2, member3])
    @test complex1 == complex2
    complex2 = SAC(apex, [member1])
    @test complex1 ≠ complex2
    complex2 = SAC(2, [member1, member2, member3])
    @test complex1 ≠ complex2
end

@testset "split_boundary_edge!" begin
    enricher_III = DT.BoundaryEnricher(deepcopy(points_III), deepcopy(curve_III))
    push!(enricher_III.points, (0.67, 0.67))
    ii, jj, rr = 7, 8, DT.num_points(enricher_III.points)
    DT.split_boundary_edge!(enricher_III, ii, jj, rr)
    @test enricher_III.boundary_nodes == [[[1, 2, 3, 4, 5], [5, 6, 7, 16, 8, 9], [9, 10, 11, 1]], [[15, 14, 13, 12], [12, 15]]]
    for (e, (pos, i)) in enricher_III.boundary_edge_map
        j = i + 1
        u = enricher_III.boundary_nodes[pos[1]][pos[2]][i]
        v = enricher_III.boundary_nodes[pos[1]][pos[2]][j]
        @test e == (u, v)
    end
    ctr = 1
    for i in eachindex(enricher_III.boundary_nodes)
        for j in eachindex(enricher_III.boundary_nodes[i])
            for k in 1:(length(enricher_III.boundary_nodes[i])-1)
                u = enricher_III.boundary_nodes[i][j][k]
                v = enricher_III.boundary_nodes[i][j][k+1]
                p = enricher_III.parent_map[(u, v)]
                @test p == ctr
            end
            ctr += 1
        end
    end
    rects, els = get_dt_rectangles(enricher_III.spatial_tree.tree)
    all_edges = Set(DT.get_edge(el) for el in els)
    @test (ii, jj) ∉ all_edges && (jj, ii) ∉ all_edges &&
          (ii, rr) ∈ all_edges && (jj, rr) ∈ all_edges
end

@testset "split_subcurve! (standard)" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        enricher_I = DT.BoundaryEnricher(deepcopy(points_I), deepcopy(curve_I))
        t, Δθ, ct = DT.compute_split_position(enricher_I, 3, 4, PT())
        @test isnan(t) && Δθ == 0.0 && ct ⪧ DT.midpoint(get_point(enricher_I.points, 3, 4)...)
        DT.split_subcurve!(enricher_I, 3, 4, PT())
        @test enricher_I.points[end] ⪧ DT.midpoint(get_point(enricher_I.points, 3, 4)...)

        enricher_II = DT.BoundaryEnricher(deepcopy(points_II), deepcopy(curve_II))
        t, Δθ, ct = DT.compute_split_position(enricher_II, 7, 8, PT())
        @test isnan(t) && Δθ ≈ 0.0 && ct ⪧ DT.midpoint(get_point(enricher_II.points, 7, 8)...)
        DT.split_subcurve!(enricher_II, 7, 8, PT())
        @test enricher_II.points[end] ⪧ DT.midpoint(get_point(enricher_II.points, 7, 8)...)

        enricher_III = DT.BoundaryEnricher(deepcopy(points_III), deepcopy(curve_III))
        t, Δθ, ct = DT.compute_split_position(enricher_III, 14, 13, PT())
        @test isnan(t) && Δθ ≈ 0.0 && ct ⪧ DT.midpoint(get_point(enricher_III.points, 14, 13)...)
        DT.split_subcurve!(enricher_III, 14, 13, PT())
        @test enricher_III.points[end] ⪧ DT.midpoint(get_point(enricher_III.points, 14, 13)...)

        enricher_IV = DT.BoundaryEnricher(deepcopy(points_IV), deepcopy(curve_IV))
        t, Δθ, ct = DT.compute_split_position(enricher_IV, 2, 7, PT())
        origθ = DT.total_variation(curve_IV[1], DT.get_inverse(curve_IV[1], get_point(enricher_IV.points, 2)), DT.get_inverse(curve_IV[1], get_point(enricher_IV.points, 7)))
        @test Δθ ≈ origθ / 2
        @test 0.5 ≤ t ≤ 0.625 && t ≈ 0.5625 && ct ⪧ curve_IV[1](t)
        DT.split_subcurve!(enricher_IV, 2, 7, PT())
        @test enricher_IV.points[end] ⪧ ct

        enricher_V = DT.BoundaryEnricher(deepcopy(points_V), deepcopy(curve_V))
        t, Δθ, ct = DT.compute_split_position(enricher_V, 3, 5, PT())
        origθ = DT.total_variation(curve_V[1], DT.get_inverse(curve_V[1], get_point(enricher_V.points, 3)), DT.get_inverse(curve_V[1], get_point(enricher_V.points, 5)))
        @test Δθ ≈ origθ / 2 rtol = 1.0e-3
        @test 0.4999 ≤ t ≤ 0.708 && t ≈ 0.5994483483424031 && ct ⪧ curve_V[1](t)
        DT.split_subcurve!(enricher_V, 3, 5, PT())
        @test enricher_V.points[end] ⪧ ct

        enricher_VI = DT.BoundaryEnricher(deepcopy(points_VI), deepcopy(curve_VI))
        t, Δθ, ct = DT.compute_split_position(enricher_VI, 12, 14, PT())
        origθ = DT.total_variation(curve_VI[1][1], DT.get_inverse(curve_VI[1][1], get_point(enricher_VI.points, 12)), DT.get_inverse(curve_VI[1][1], get_point(enricher_VI.points, 14)))
        @test Δθ ≈ origθ / 2 rtol = 1.0e-3
        @test 0.5 ≤ t ≤ 0.75 && t ≈ 0.625 && ct ⪧ curve_VI[1][1](t)
        DT.split_subcurve!(enricher_VI, 12, 14, PT())
        @test enricher_VI.points[end] ⪧ ct

        enricher_VII = DT.BoundaryEnricher(deepcopy(points_VII), deepcopy(curve_VII))
        t, Δθ, ct = DT.compute_split_position(enricher_VII, 8, 7, PT())
        origθ = DT.total_variation(curve_VII[2][1], DT.get_inverse(curve_VII[2][1], get_point(enricher_VII.points, 8)), DT.get_inverse(curve_VII[2][1], get_point(enricher_VII.points, 7)))
        @test Δθ ≈ origθ / 2 rtol = 1.0e-3
        @test 0.1 ≤ t ≤ 0.501 && t ≈ 0.1599962634542922 && ct ⪧ curve_VII[2][1](t)
        DT.split_subcurve!(enricher_VII, 8, 7, PT())
        @test enricher_VII.points[end] ⪧ ct

        enricher_VIII = DT.BoundaryEnricher(deepcopy(points_VIII), deepcopy(curve_VIII))
        t, Δθ, ct = DT.compute_split_position(enricher_VIII, 7, 8, PT())
        @test Δθ ≈ 0.0 && isnan(t) && ct ⪧ DT.midpoint(get_point(enricher_VIII.points, 7, 8)...)
        DT.split_subcurve!(enricher_VIII, 7, 8, PT())
        @test enricher_VIII.points[end] ⪧ DT.midpoint(get_point(enricher_VIII.points, 7, 8)...)
    end
end

@testset "split_subcurve! (SmallAngleComplex)" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        points, boundary_nodes = deepcopy(points_II), deepcopy(curve_II)
        enricher = DT.BoundaryEnricher(points, boundary_nodes)
        emin = DT.get_minimum_edge_length(enricher.small_angle_complexes[9][1], enricher.points)
        ℓ = DT.balanced_power_of_two_ternary_split(emin)
        DT.split_subcurve!(enricher, 8, 9, PT())
        @test DT.dist(enricher.points[9], enricher.points[13]) ≈ DT.dist(enricher.points[9], enricher.points[12]) ≈ ℓ
        _complexes = Dict(
            9 => [SAC(9, [SACM(2, 12), SACM(3, 13)])],
            1 => [SAC(1, [SACM(3, 11), SACM(1, 2)])],
        )
        @test DT.get_small_angle_complexes(enricher) == _complexes

        A, B, C, D, E, F, G, H, I, J, K = (0.0, 0.0), (0.2, 1.4), (0.6, 1.2),
        (1.2, 0.2), (1.2, -0.2), (-1.4, -0.2),
        (-1.0, -0.6), (0.6, 1.0), (0.8, 0.6),
        (0.6, 0.4), (0.6, 0.2)
        points = [A, B, C, D, E, F, G, H, I, J, K]
        boundary_nodes = [[[1, 3, 2, 1]], [[1, 9, 8, 1]], [[1, 11, 10, 1]], [[1, 5, 4, 1]], [[1, 6, 7, 1]]]
        enricher = DT.BoundaryEnricher(points, boundary_nodes)
        DT.split_subcurve!(enricher, 1, 3, PT())
        points = get_points(enricher)
        @test all(≈(DT.get_minimum_edge_length(enricher.small_angle_complexes[1][1], points)), DT.dist.(Ref(points[1]), get_point(points, 12, 13, 14, 15, 16, 17, 18, 19)))
        _complexes = Dict(
            1 => [
                SAC(1, [SACM(1, 12), SACM(1, 13), SACM(2, 14), SACM(2, 15), SACM(3, 16), SACM(3, 17), SACM(4, 18), SACM(4, 19)]),
                SAC(1, [SACM(5, 7), SACM(5, 6)]),
            ],
        )
        @test DT.get_small_angle_complexes(enricher) == _complexes
        DT.split_subcurve!(enricher, 1, 6, PT())
        @test all(≈(DT.get_minimum_edge_length(enricher.small_angle_complexes[1][2], points)), DT.dist.(Ref(points[1]), get_point(points, 20, 21)))
        _complexes = Dict(
            1 => [
                SAC(1, [SACM(1, 12), SACM(1, 13), SACM(2, 14), SACM(2, 15), SACM(3, 16), SACM(3, 17), SACM(4, 18), SACM(4, 19)]),
                SAC(1, [SACM(5, 20), SACM(5, 21)]),
            ],
        )
        @test DT.get_small_angle_complexes(enricher) == _complexes
    end
end

@testset "has_acute_neighbouring_angles and splitting small angles" begin
    A, B, C, D, E, F, G, H, I, J, K = (0.0, 0.0), (7.0, 0.0), (-0.2, 1.0), (1.0, 0.0),
    (2.0, 0.0), (3.0, 0.0), (4.0, 0.0), (4.5, 0.5), (3.0, 0.8),
    (1.6, 0.8), (0.4, 0.2)
    points = [E, F, G, B, H, I, C, A, D, E]
    boundary_nodes, points = convert_boundary_points_to_indices(points)
    enricher = DT.BoundaryEnricher(points, boundary_nodes)
    num_adjoin, adjoin_vert = DT.has_acute_neighbouring_angles(rt(), enricher, 1, 2)
    @inferred DT.has_acute_neighbouring_angles(rt(), enricher, 1, 2)
    @test num_adjoin == 0 && adjoin_vert == 0
    num_adjoin, adjoin_vert = DT.has_acute_neighbouring_angles(rt(), enricher, 2, 3)
    @test num_adjoin == adjoin_vert == 0
    num_adjoin, adjoin_vert = DT.has_acute_neighbouring_angles(rt(), enricher, 3, 4)
    @test num_adjoin == 1 && adjoin_vert == 4
    num_adjoin, adjoin_vert = DT.has_acute_neighbouring_angles(rt(), enricher, 4, 5)
    @test num_adjoin == 1 && adjoin_vert == 4
    num_adjoin, adjoin_vert = DT.has_acute_neighbouring_angles(rt(), enricher, 5, 6)
    @test num_adjoin == adjoin_vert == 0
    num_adjoin, adjoin_vert = DT.has_acute_neighbouring_angles(rt(), enricher, 6, 7)
    @test num_adjoin == 1 && adjoin_vert == 7
    num_adjoin, adjoin_vert = DT.has_acute_neighbouring_angles(rt(), enricher, 7, 8)
    @test num_adjoin == 1 && adjoin_vert == 7
    num_adjoin, adjoin_vert = DT.has_acute_neighbouring_angles(rt(), enricher, 8, 9)
    @test num_adjoin == adjoin_vert == 0
    t12, Δθ12, ct12 = DT.compute_split_position(enricher, 1, 2, rt())
    @test (t12, Δθ12, ct12)[2:3] == DT._compute_split_position_standard(enricher, 1, 2)[2:3]
    t23, Δθ23, ct23 = DT.compute_split_position(enricher, 2, 3, rt())
    @test (t23, Δθ23, ct23)[2:3] == DT._compute_split_position_standard(enricher, 2, 3)[2:3]
    t34, Δθ34, ct34 = DT.compute_split_position(enricher, 3, 4, rt())
    @test (t34, Δθ34, ct34) == DT._compute_split_position_acute(enricher, 3, 4, 1, 4)
    @test ispow2(t34 * DT.dist(DT.get_point(enricher.points, 3, 4)...))
    t45, Δθ45, ct45 = DT.compute_split_position(enricher, 4, 5, rt())
    @test (t45, Δθ45, ct45) == DT._compute_split_position_acute(enricher, 4, 5, 1, 4)
    @test ispow2(t45 * DT.dist(DT.get_point(enricher.points, 4, 5)...))
    t56, Δθ56, ct56 = DT.compute_split_position(enricher, 5, 6, rt())
    @test (t56, Δθ56, ct56)[2:3] == DT._compute_split_position_standard(enricher, 5, 6)[2:3]
    t67, Δθ67, ct67 = DT.compute_split_position(enricher, 6, 7, rt())
    @test (t67, Δθ67, ct67) == DT._compute_split_position_acute(enricher, 6, 7, 1, 7)
    @test ispow2((1 - t67) * DT.dist(DT.get_point(enricher.points, 6, 7)...))
    t78, Δθ78, ct78 = DT.compute_split_position(enricher, 7, 8, rt())
    @test (t78, Δθ78, ct78) == DT._compute_split_position_acute(enricher, 7, 8, 1, 7)
    @test ispow2(t78 * DT.dist(DT.get_point(enricher.points, 7, 8)...))
    t89, Δθ89, ct89 = DT.compute_split_position(enricher, 8, 9, rt())
    @test (t89, Δθ89, ct89)[2:3] == DT._compute_split_position_standard(enricher, 8, 9)[2:3]

    A, B, C, D, E = (0.0, 0.0), (3.0, 0.0), (0.5, 1.0), (2.5, 1.0), (1.5, 1.5)
    points = [A, B, C, D, E]
    boundary_nodes = [5, 3, 1, 2, 4, 5]
    enricher = DT.BoundaryEnricher(points, boundary_nodes)
    num_adjoin, adjoin_vert = DT.has_acute_neighbouring_angles(rt(), enricher, 1, 2)
    @test num_adjoin == 2 && adjoin_vert == 0
    num_adjoin, adjoin_vert = DT.has_acute_neighbouring_angles(rt(), enricher, 3, 1)
    @test num_adjoin == 1 && adjoin_vert == 1
    t12, Δθ12, ct12 = DT.compute_split_position(enricher, 1, 2, rt())
    @test (t12, Δθ12, ct12) == DT._compute_split_position_acute(enricher, 1, 2, 2, 0)
    @test t12 == DT.compute_concentric_shell_quarternary_split_position(get_point(enricher.points, 1, 2)...)
    t23, Δθ23, ct23 = DT.compute_split_position(enricher, 3, 1, rt())
    @test (t23, Δθ23, ct23) == DT._compute_split_position_acute(enricher, 3, 1, 1, 1)
    @test (1 - t23) * DT.dist(DT.get_point(enricher.points, 3, 1)...) ≈ 1 / 2
end

@testset "enrich_boundary!" begin
    @testset "enrich_boundary! (no points or extra segments)" begin
        point_sets = [points_I, points_II, points_III, points_IV, points_V, points_VI, points_VII, points_VIII, points_IX, points_X, points_XI, points_XII]
        curve_sets = [curve_I, curve_II, curve_III, curve_IV, curve_V, curve_VI, curve_VII, curve_VIII, curve_IX, curve_X, curve_XI, curve_XII]
        for (points, curve) in zip(point_sets, curve_sets)
            enricher = DT.BoundaryEnricher(deepcopy(points), deepcopy(curve))
            DT.enrich_boundary!(enricher; predicates=rt())
            @test all_diametral_circles_are_empty(enricher) == 1
            @test allunique(get_points(enricher))
            @test all_points_are_inside(enricher, points, curve)
        end

        A, B, C, D, E, F, G, H, I, J, K = (0.0, 0.0), (7.0, 0.0), (2.0, 1.0), (1.0, 0.0),
        (2.0, 0.0), (3.0, 0.0), (4.0, 0.0), (4.5, 0.5), (3.0, 0.8),
        (1.6, 0.8), (0.4, 0.2)
        points = [E, F, G, B, H, I, C, J, K, A, D, E]
        boundary_nodes, points = convert_boundary_points_to_indices(points)
        enricher = DT.BoundaryEnricher(points, boundary_nodes)
        DT.enrich_boundary!(enricher; predicates=rt())
        @test all_diametral_circles_are_empty(enricher) == 1
        @test all_points_are_inside(enricher, points, boundary_nodes)
    end

    @testset "enrich_boundary! (extra points, no segments)" begin
        point_sets = [points_I_extra, points_II_extra, points_III_extra, points_IV_extra, points_V_extra, points_VI_extra, points_VII_extra, points_VIII_extra, points_IX_extra, points_X_extra, points_XI_extra, points_XII_extra]
        curve_sets = [curve_I, curve_II, curve_III, curve_IV, curve_V, curve_VI, curve_VII, curve_VIII, curve_IX, curve_X, curve_XI, curve_XII]
        for (points, curve) in zip(point_sets, curve_sets)
            enricher = DT.BoundaryEnricher(deepcopy(points), deepcopy(curve))
            DT.enrich_boundary!(enricher; predicates=rt())
            @test all_diametral_circles_are_empty(enricher) == 1
            @test allunique(get_points(enricher))
            @test all_points_are_inside(enricher, points, curve)
        end
    end

    @testset "enrich_boundary! (extra points, extra segments)" begin
        for PT in subtypes(DT.AbstractPredicateKernel)
            point_sets = [points_I_extra_segments, points_II_extra_segments, points_III_extra_segments, points_IV_extra_segments]
            curve_sets = [curve_I, curve_II, curve_III, curve_IV]
            segment_sets = [segments_I, segments_II, segments_III, segments_IV]
            for i in eachindex(point_sets)
                enricher = DT.BoundaryEnricher(deepcopy(point_sets[i]), deepcopy(curve_sets[i]), deepcopy(segment_sets[i]))
                DT.enrich_boundary!(enricher; predicates=PT())
                @test all_diametral_circles_are_empty(enricher) == 1
                @test allunique(get_points(enricher))
                @test all_points_are_inside(enricher, point_sets[i], curve_sets[i])
            end
        end
    end
end

@testset "triangulate_curve_bounded" begin
    @testset "triangulate_curve_bounded (no points or extra segments)" begin
        for PT in (DT.ExactKernel, DT.AdaptiveKernel)
            point_sets = deepcopy.([points_I, points_II, points_III, points_IV, points_V, points_VI, points_VII, points_VIII, points_IX, points_X, points_XI, points_XII])
            curve_sets = deepcopy.([curve_I, curve_II, curve_III, curve_IV, curve_V, curve_VI, curve_VII, curve_VIII, curve_IX, curve_X, curve_XI, curve_XII])
            for i in eachindex(point_sets, curve_sets)
                points, curve = deepcopy(point_sets[i]), deepcopy(curve_sets[i])
                tri = triangulate(points; boundary_nodes=curve, enrich=i ≤ 3, predicates=PT())
                @test validate_triangulation(tri)
                @test is_conformal(tri; predicates=PT())
                @test DT.get_boundary_enricher(tri) == DT.enrich_boundary!(DT.BoundaryEnricher(deepcopy(point_sets[i]), deepcopy(curve_sets[i])), predicates=PT())
                i > 3 && @test DT.is_curve_bounded(tri)
                i ≤ 3 && @test !DT.is_curve_bounded(tri)
            end
        end
    end

    @testset "triangulate_curve_bounded (extra points, no segments)" begin
        for PT in (DT.ExactKernel, DT.AdaptiveKernel)
            point_sets = deepcopy.([points_I_extra, points_II_extra, points_III_extra, points_IV_extra, points_V_extra, points_VI_extra, points_VII_extra, points_VIII_extra, points_IX_extra, points_X_extra, points_XI_extra, points_XII_extra])
            curve_sets = deepcopy.([curve_I, curve_II, curve_III, curve_IV, curve_V, curve_VI, curve_VII, curve_VIII, curve_IX, curve_X, curve_XI, curve_XII])
            for i in eachindex(point_sets, curve_sets)
                points, curve = deepcopy(point_sets[i]), deepcopy(curve_sets[i])
                tri = triangulate(points; boundary_nodes=curve, enrich=i ≤ 3, predicates=PT())
                @test validate_triangulation(tri)
                @test is_conformal(tri; predicates=PT())
                i ∉ (2, 11) && @test DT.get_boundary_enricher(tri) == DT.enrich_boundary!(DT.BoundaryEnricher(deepcopy(point_sets[i]), deepcopy(curve_sets[i])), predicates=PT()) # i ≠ 2 since we deliberately included some boundary points in the extra points, which triangulate then sees and mutates
                i > 3 && @test DT.is_curve_bounded(tri)
                i ≤ 3 && @test !DT.is_curve_bounded(tri)
            end
        end
    end

    @testset "triangulate_curve_bounded (extra points, extra segments)" begin
        for PT in (DT.ExactKernel, DT.AdaptiveKernel)
            point_sets = deepcopy.([points_I_extra_segments, points_II_extra_segments, points_III_extra_segments, points_IV_extra_segments])
            curve_sets = deepcopy.([curve_I, curve_II, curve_III, curve_IV])
            segment_sets = deepcopy.([segments_I, segments_II, segments_III, segments_IV])
            for i in eachindex(point_sets, curve_sets, segment_sets)
                points, curve, segments = deepcopy(point_sets[i]), deepcopy(curve_sets[i]), deepcopy(segment_sets[i])
                tri = triangulate(points; boundary_nodes=curve, segments=segments, enrich=i ≤ 3, predicates=PT())
                @test validate_triangulation(tri)
                @test is_conformal(tri; predicates=PT())
                i ≠ 2 && @test DT.get_boundary_enricher(tri) == DT.enrich_boundary!(DT.BoundaryEnricher(deepcopy(point_sets[i]), deepcopy(curve_sets[i]), deepcopy(segment_sets[i])), predicates=PT())
                i > 3 && @test DT.is_curve_bounded(tri)
                i ≤ 3 && @test !DT.is_curve_bounded(tri)
            end
        end
    end
end

@testset "refine_curve_bounded with circumcenters" begin
    point_sets_no_extra = deepcopy.([points_I, points_II, points_III, points_IV, points_V, points_VI, points_VII, points_VIII, points_IX, points_X, points_XI, points_XII])
    point_sets_extra_points = deepcopy.([points_I_extra, points_II_extra, points_III_extra, points_IV_extra, points_V_extra, points_VI_extra, points_VII_extra, points_VIII_extra, points_IX_extra, points_X_extra, points_XI_extra, points_XII_extra])
    point_sets_extra_segments = deepcopy.([points_I_extra_segments, points_II_extra_segments, points_III_extra_segments, points_IV_extra_segments])
    segment_sets = deepcopy.([segments_I, segments_II, segments_III, segments_IV])
    point_sets = (point_sets_no_extra, point_sets_extra_points, point_sets_extra_segments)
    curve_sets = deepcopy.([curve_I, curve_II, curve_III, curve_IV, curve_V, curve_VI, curve_VII, curve_VIII, curve_IX, curve_X, curve_XI, curve_XII])
    point_names = ("default", "extra_points", "extra_segments")
    _rng_num(
        idx1,
        idx2, idx3, idx4, idx5, curve_idx, point_idx,
    ) = 2^idx1 * 3^idx2 * 5^idx3 * 7^idx4 * 11^idx5 * 13^curve_idx * 17^point_idx

    @testset "all_examples" begin
        for PT in (DT.ExactKernel, DT.AdaptiveKernel)
            max_area_opts = [
                (1.0e-2, 1.0e-3),
                (1.0e-2, 1.0e-3),
                (1.0e-2, 1.0e-3),
                (1.0e-2, 1.0e-3),
                (1.0e-2, 1.0e-3),
                (1.0e-1, 1.0e-2),
                (1.0e-1, 1.0e-2),
                (1.0e-1, 1.0e-2),
                (1.0e-2, 1.0e-3),
                (1.0e-1, 1.0e-2),
                (1.0e-1, 1.0e-2),
                (1.0e-1, 1.0e-2),
            ]
            for curve_idx in 4:lastindex(curve_sets)
                for point_idx in 1:3
                    point_idx == 3 && curve_idx ≥ 5 && continue # no extra segments for curves ≥ 5
                    for (idx1, use_lens) in enumerate((false, true))
                        for (idx2, min_angle) in enumerate((20.0, 27.5, 30.0))
                            for (idx3, min_area) in enumerate((1.0e-12,))
                                for (idx4, max_area) in enumerate(max_area_opts[curve_idx])
                                    for (idx5, seditious_angle) in enumerate((10.0, 20.0))
                                        @info "Testing curve-bounded refinement with circumcenters. use_lens: $use_lens; min_angle: $min_angle; min_area: $min_area; max_area: $max_area; seditious_angle: $seditious_angle; curve: $curve_idx; point set: $point_idx"
                                        rng = StableRNG(abs(_rng_num(idx1, idx2, idx3, idx4, idx5, curve_idx, point_idx)))
                                        points, curve = deepcopy(point_sets[point_idx][curve_idx]), deepcopy(curve_sets[curve_idx])
                                        if point_idx ≤ 2
                                            tri = triangulate(points; boundary_nodes=curve, enrich=curve_idx ≤ 3, rng, predicates=PT())
                                        else
                                            segments = deepcopy(segment_sets[curve_idx])
                                            tri = triangulate(points; boundary_nodes=curve, segments=segments, enrich=curve_idx ≤ 3, rng, predicates=PT())
                                        end
                                        custom_constraint = (_tri, T) -> curve_idx ≠ 5 ? false : begin
                                            i, j, k = triangle_vertices(T)
                                            p, q, r = get_point(_tri, i, j, k)
                                            c = (p .+ q .+ r) ./ 3
                                            x, y = getxy(c)
                                            return (x + y^2 < 1 / 4) && DT.triangle_area(p, q, r) > 1.0e-4 / 2
                                        end
                                        refine!(tri; min_angle, min_area, max_area, custom_constraint, seditious_angle, use_circumcenter=true, use_lens, rng, predicates=PT())
                                        args = DT.RefinementArguments(tri; min_angle, min_area, max_area, seditious_angle, custom_constraint, use_circumcenter=true, use_lens, predicates=PT())
                                        @test validate_refinement(tri, args, warn=false)
                                        if _rng_num(idx1, idx2, idx3, idx4, idx5, curve_idx, point_idx) == _rng_num(1, 3, 1, 2, 2, curve_idx, point_idx)
                                            fig, ax, sc = triplot(tri)
                                            @test_reference "refine_curve_bounded_example_$(curve_idx)_$(names[curve_idx])_$(point_names[point_idx])_$(abs(_rng_num(1, 3, 1, 2, 2, curve_idx, point_idx))).png" fig by = psnr_equality(7)
                                        elseif _rng_num(idx1, idx2, idx3, idx4, idx5, curve_idx, point_idx) == _rng_num(2, 3, 1, 2, 2, curve_idx, point_idx)
                                            fig, ax, sc = triplot(tri)
                                            @test_reference "refine_curve_bounded_example_$(curve_idx)_$(names[curve_idx])_$(point_names[point_idx])_$(abs(_rng_num(2, 3, 1, 2, 2, curve_idx, point_idx))).png" fig by = psnr_equality(7)
                                        end
                                        @test tri.boundary_enricher.boundary_edge_map == tri.boundary_edge_map
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "algorithm_terminates appropriately" begin
        for PT in (DT.ExactKernel, DT.AdaptiveKernel)
            curve_idx = 12
            point_idx = 2
            points, curve = deepcopy(point_sets[point_idx][curve_idx]), deepcopy(curve_sets[curve_idx])
            tri = triangulate(points; boundary_nodes=curve, predicates=PT())
            refine!(tri; max_area=1.0e-3, max_points=500, use_circumcenter=true, predicates=PT())
            @test DT.num_solid_vertices(tri) == 500
            @test validate_triangulation(tri)
            @test !validate_refinement(tri; max_area=1.0e-3, max_points=500, use_circumcenter=true, warn=false)
        end
    end
end

@testset "adding segment to a curve-bounded domain with no existing segments" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        curve = [
            [
                [1, 2, 3], [EllipticalArc((2.0, 0.0), (-2.0, 0.0), (0.0, 0.0), 2, 1 / 2, 0.0)],
            ],
            [
                [BSpline([(0.0, 0.4), (1.0, 0.2), (0.0, 0.1), (-1.0, 0.2), (0.0, 0.4)])],
            ],
            [
                [4, 5, 6, 7, 4],
            ],
            [
                [BezierCurve([(0.0, -2.0), (0.0, -2.5), (-1.0, -2.5), (-1.0, -3.0)])], [CatmullRomSpline([(-1.0, -3.0), (0.0, -4.0), (1.0, -3.0), (0.0, -2.0)])],
            ],
            [
                [12, 11, 10, 12],
            ],
            [
                [CircularArc((1.1, -3.0), (1.1, -3.0), (0.0, -3.0), positive=false)],
            ],
        ]
        points = [(-2.0, 0.0), (0.0, 0.0), (2.0, 0.0), (-2.0, -5.0), (2.0, -5.0), (2.0, -1 / 10), (-2.0, -1 / 10), (-1.0, -3.0), (0.0, -4.0), (0.0, -2.3), (-0.5, -3.5), (0.9, -3.0)]
        rng = StableRNG(123)
        tri = triangulate(points; boundary_nodes=curve, rng, predicates=PT())
        refine!(tri; max_area=1.0e-2, predicates=PT())
        @test validate_triangulation(tri, predicates=PT())
        r = DT.num_points(tri)
        add_point!(tri, -3 / 2, -4.0, concavity_protection=true, predicates=PT())
        add_point!(tri, -3 / 2, -1.0, concavity_protection=true, predicates=PT())
        @test validate_triangulation(tri, predicates=PT())
        add_segment!(tri, r + 1, r + 2, predicates=PT())
        @test validate_triangulation(tri, predicates=PT())
        @test tri.boundary_enricher.segments ∈ (Set(((r + 1, r + 2),)), Set(((r + 2, r + 1),)))
    end
end

@testset "Custom structs" begin
    # Used to be a broken example
    DT = DelaunayTriangulation
    struct Custom2Points
        points::Vector{NTuple{2,Float64}}
    end
    struct Custom2Segments
        segments::Set{NTuple{2,Int}}
    end
    struct Custom2Triangles
        triangles::Set{NTuple{3,Int}}
    end
    Base.eachindex(points::Custom2Points) = Base.eachindex(points.points)
    Base.iterate(points::Custom2Points, state...) = Base.iterate(points.points, state...)
    Base.length(points::Custom2Points) = length(points.points)
    Base.getindex(points::Custom2Points, i) = points.points[i]
    DT.number_type(::Type{Custom2Points}) = Float64
    Base.iterate(triangles::Custom2Triangles, state...) = Base.iterate(triangles.triangles, state...)
    Base.sizehint!(triangles::Custom2Triangles, n) = sizehint!(triangles.triangles, n)
    Base.eltype(::Type{Custom2Triangles}) = NTuple{3,Int}
    Base.length(triangles::Custom2Triangles) = length(triangles.triangles)
    Custom2Triangles() = Custom2Triangles(Set{NTuple{3,Int}}())
    Base.iterate(segments::Custom2Segments, state...) = Base.iterate(segments.segments, state...)
    Base.eltype(::Type{Custom2Segments}) = NTuple{2,Int}
    Base.rand(rng::AbstractRNG, segments::Custom2Segments) = rand(rng, segments.segments)
    Base.length(segments::Custom2Segments) = length(segments.segments)
    Custom2Segments() = Custom2Segments(Set{NTuple{2,Int}}())
    DT.contains_edge(e::NTuple{2,Int}, Es::Custom2Segments) = e ∈ Es.segments
    Base.empty!(triangles::Custom2Triangles) = empty!(triangles.triangles)
    Base.empty!(segments::Custom2Segments) = empty!(segments.segments)
    Base.push!(segments::Custom2Segments, e::NTuple{2,Int}) = push!(segments.segments, e)
    Base.delete!(segments::Custom2Segments, e::NTuple{2,Int}) = delete!(segments.segments, e)
    Base.pop!(points::Custom2Points) = pop!(points.points)
    DT.push_point!(points::Custom2Points, x, y) = push!(points.points, (x, y))
    Base.delete!(triangles::Custom2Triangles, T::NTuple{3,Int}) = delete!(triangles.triangles, T)
    DT.set_point!(points::Custom2Segments, i, x, y) = points.points[i] = (x, y)
    Base.push!(triangles::Custom2Triangles, T::NTuple{3,Int}) = push!(triangles.triangles, T)
    points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    append!(points, [(clamp(0.5 + randn(), -2, 2), clamp(0.5 + randn(), -2, 2)) for _ in 1:100])
    inner_circle = CircularArc((0.5, 0.25), (0.5, 0.25), (0.5, 0.5), positive=false)
    boundary_nodes = [[[1, 2, 3, 4, 1]], [[inner_circle]]]
    points = Custom2Points(points)
    tri = triangulate(points; boundary_nodes, TrianglesType=Custom2Triangles, EdgesType=Custom2Segments)
    @test DT.validate_triangulation(tri)
    refine!(tri; max_area = 0.2)
    @test DT.validate_triangulation(tri)
end