using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using DataStructures
using Random
using LinearAlgebra
using ReferenceTests
using Test
using StructEquality
using DelimitedFiles
using StableRNGs

for _ in 1:20 # Run many times to make sure the segfault is gone
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
    curve_I = [1, 2, 3, 4, 1]
    points_I = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]

    curve_II = [[1, 2, 3, 4, 5], [5, 6, 7, 8, 9], [9, 10, 11, 1]]
    points_II = [
        (0.0, 0.0), (0.25, 0.0), (0.5, 0.0), (0.75, 0.0), (1.0, 0.0),
        (1.0, 0.25), (1.0, 0.5), (1.0, 0.75), (1.0, 1.0),
        (0.75, 0.75), (0.25, 0.25)
    ]

    curve_III = [[[1, 2, 3, 4, 5], [5, 6, 7, 8, 9], [9, 10, 11, 1]], [[15, 14, 13, 12], [12, 15]]]
    points_III = [
        (0.0, 0.0), (0.25, 0.0), (0.5, 0.0), (0.75, 0.0), (1.0, 0.0),
        (1.0, 0.25), (1.0, 0.5), (1.0, 0.75), (1.0, 1.0),
        (0.0, 1.0), (0.0, 0.5),
        (0.25, 0.25), (0.75, 0.25), (0.75, 0.75), (0.25, 0.75), (0.5, 0.5)
    ]

    curve_IV = [CircularArc((1.0, 0.0), (1.0, 0.0), (0.0, 0.0))]
    points_IV = NTuple{2,Float64}[]

    curve_V = [BezierCurve([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)])]
    points_V = [(0.0, 0.0), (0.2, 0.25)]

    curve_VI = [
        [CircularArc((1.0, 0.0), (0.0, 1.0), (0.0, 0.0))],
        [BSpline([(0.0, 1.0), (-1.0, 2.0), (-2.0, 0.0), (-2.0, -1.0), (0.0, -2.0)])],
        [5, 6, 10]
    ]
    points_VI = [(0.1, 0.1), (0.15, 0.15), (0.23, 0.23), (0.009, 0.11), (0.0, -2.0), (0.2, -1.7), (0.000591, 0.00019), (0.111, -0.005), (-0.0001, -0.00991), (1.0, 0.0)]

    curve_VII = [
        [CircularArc((2.0, 0.0), (-2.0, 0.0), (0.0, 0.0))],
        [BSpline([(-2.0, 0.0), (-2.0, -1.0), (0.0, -1.0), (1.0, -1.0), (2.0, -1.0), (2.0, 0.0)])]
    ]
    points_VII = [(2.0, 0.0), (0.0, 0.5)]

    curve_VIII = [
        [1, 2, 3, 4, 5],
        [DT.EllipticalArc((0.0, 0.0), (2.0, -2.0), (1.0, -1.0), sqrt(2), sqrt(2), 45.0)],
        [6, 7, 8, 9, 10],
        [CatmullRomSpline([(10.0, -3.0), (20.0, 0.0), (18.0, 0.0), (10.0, 0.0)])]
    ]
    points_VIII = [(10.0, 0.0), (8.0, 0.0), (4.0, 0.0), (2.0, 2.0), (0.0, 0.0), (2.0, -2.0),
        (2.5, -2.0), (3.5, -2.0), (4.5, -3.0), (10.0, -3.0), (10.0, 12.0), (14.0, 0.0)]

    curve_IX =
        [
            [
                [1, 2, 3, 4, 5, 6, 7, 1]
            ],
            [
                [CircularArc((0.6, 0.5), (0.6, 0.5), (0.5, 0.5), positive=false)]
            ],
        ]
    points_IX = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.5, 1.5), (0.0, 1.0), (0.0, 0.5), (0.0, 0.2)]

    curve_X = [
        [
            [1, 2, 3], [DT.EllipticalArc((2.0, 0.0), (-2.0, 0.0), (0.0, 0.0), 2, 1 / 2, 0.0)]
        ],
        [
            [BSpline(reverse([(1.0, 0.2), (0.0, 0.4), (0.0, 0.3), (-1.0, 0.2)]))], reverse([4, 5, 6, 7, 8])
        ]
    ]
    points_X = [
        (-2.0, 0.0), (0.0, 0.0), (2.0, 0.0), (-1.0, 0.2), (-1.0, 0.1), (0.0, 0.1), (1.0, 0.1), (1.0, 0.2)
    ]

    curve_XI = [
        [
            [1, 2, 3], [DT.EllipticalArc((2.0, 0.0), (-2.0, 0.0), (0.0, 0.0), 2, 1 / 2, 0.0)]
        ],
        [
            [BSpline([(0.0, 0.4), (1.0, 0.2), (0.0, 0.1), (-1.0, 0.2), (0.0, 0.4)])]
        ],
        [
            [4, 5, 6, 7, 4]
        ],
        [
            [BezierCurve(reverse([(-1.0, -3.0), (-1.0, -2.5), (0.0, -2.5), (0.0, -2.0)]))], [CatmullRomSpline(reverse([(0.0, -2.0), (1.0, -3.0), (0.0, -4.0), (-1.0, -3.0)]))]
        ],
        [
            [12, 11, 10, 12]
        ],
        [
            [CircularArc((1.1, -3.0), (1.1, -3.0), (0.0, -3.0), positive=false)]
        ]
    ]
    points_XI = [(-2.0, 0.0), (0.0, 0.0), (2.0, 0.0), (-2.0, -5.0), (2.0, -5.0), (2.0, -1 / 10), (-2.0, -1 / 10), (-1.0, -3.0), (0.0, -4.0), (0.0, -2.3), (-0.5, -3.5), (0.9, -3.0)]

    @testset "Points" begin
        points = rand(2, 500)
        hierarchy = DT.construct_polygon_hierarchy(points)
        @test DT.get_polygon_orientations(hierarchy) ⊢ BitVector([true])
        @test DT.get_polygon_orientation(hierarchy, 1)
        @test DT.get_bounding_boxes(hierarchy) ⊢ [DT.bounding_box(points)]
        @test DT.get_bounding_box(hierarchy, 1) ⊢ DT.bounding_box(points)
        @test DT.get_exterior_curve_indices(hierarchy) == Set(1)
        @test DT.get_trees(hierarchy) === hierarchy.trees
        @test DT.get_tree(hierarchy, 1) === hierarchy.trees[1]
        tree = DT.get_tree(hierarchy, 1)
        @test isnothing(DT.get_parent(tree))
        @test DT.get_children(tree) ⊢ Set{DT.PolygonTree{Int}}()
        @test DT.get_index(tree) == 1
        @test DT.get_height(tree) == 0
        @test !DT.has_parent(tree)
        @test !DT.has_children(tree)
    end

    @testset "Points with a piecewise linear boundary" begin
        points = deepcopy(points_I)
        boundary_nodes = deepcopy(curve_I)
        hierarchy = DT.construct_polygon_hierarchy(points, boundary_nodes)
        @test DT.get_polygon_orientations(hierarchy) ⊢ BitVector([true])
        @test DT.get_polygon_orientation(hierarchy, 1)
        @test DT.get_bounding_boxes(hierarchy) ⊢ [DT.BoundingBox(DT.polygon_bounds(points, boundary_nodes)...)]
        @test DT.get_bounding_box(hierarchy, 1) ⊢ DT.BoundingBox(DT.polygon_bounds(points, boundary_nodes)...)
        @test DT.get_exterior_curve_indices(hierarchy) == Set(1)
        @test DT.get_trees(hierarchy) === hierarchy.trees
        @test DT.get_tree(hierarchy, 1) === hierarchy.trees[1]
        tree = DT.get_tree(hierarchy, 1)
        @test isnothing(DT.get_parent(tree)) && !DT.has_parent(tree)
        @test DT.get_children(tree) ⊢ Set{DT.PolygonTree{Int}}()
        @test DT.get_index(tree) == 1
        @test DT.get_height(tree) == 0
        @test !DT.has_children(tree)

        points = deepcopy(points_II)
        boundary_nodes = deepcopy(curve_II)
        hierarchy = DT.construct_polygon_hierarchy(points, boundary_nodes)
        @test DT.get_polygon_orientations(hierarchy) ⊢ BitVector([true])
        @test DT.get_polygon_orientation(hierarchy, 1)
        @test DT.get_bounding_boxes(hierarchy) ⊢ [DT.BoundingBox(DT.polygon_bounds(points, boundary_nodes)...)]
        @test DT.get_bounding_box(hierarchy, 1) ⊢ DT.BoundingBox(DT.polygon_bounds(points, boundary_nodes)...)
        @test DT.get_exterior_curve_indices(hierarchy) == Set(1)
        @test DT.get_trees(hierarchy) === hierarchy.trees
        @test DT.get_tree(hierarchy, 1) === hierarchy.trees[1]
        tree = DT.get_tree(hierarchy, 1)
        @test isnothing(DT.get_parent(tree)) && !DT.has_parent(tree)
        @test DT.get_children(tree) ⊢ Set{DT.PolygonTree{Int}}()
        @test DT.get_index(tree) == 1
        @test DT.get_height(tree) == 0
        @test !DT.has_children(tree)

        points = deepcopy(points_III)
        boundary_nodes = deepcopy(curve_III)
        hierarchy = DT.construct_polygon_hierarchy(points, boundary_nodes)
        @test DT.get_polygon_orientations(hierarchy) ⊢ BitVector([true, false])
        @test DT.get_polygon_orientation(hierarchy, 1) && !DT.get_polygon_orientation(hierarchy, 2)
        @test DT.get_bounding_boxes(hierarchy) ⊢ [DT.BoundingBox(DT.polygon_bounds(points, get_boundary_nodes(boundary_nodes, 1))...), DT.BoundingBox(DT.polygon_bounds(points, get_boundary_nodes(boundary_nodes, 2))...)]
        @test DT.get_bounding_box(hierarchy, 1) ⊢ DT.BoundingBox(DT.polygon_bounds(points, get_boundary_nodes(boundary_nodes, 1))...)
        @test DT.get_bounding_box(hierarchy, 2) ⊢ DT.BoundingBox(DT.polygon_bounds(points, get_boundary_nodes(boundary_nodes, 2))...)
        @test DT.get_exterior_curve_indices(hierarchy) == Set(1)
        @test DT.get_trees(hierarchy) === hierarchy.trees
        @test DT.get_tree(hierarchy, 1) === hierarchy.trees[1]
        tree = DT.get_tree(hierarchy, 1)
        @test isnothing(DT.get_parent(tree)) && !DT.has_parent(tree)
        @test DT.has_children(tree) && length(DT.get_children(tree)) == 1
        @test DT.get_index(tree) == 1
        @test DT.get_height(tree) == 0
        children = DT.get_children(tree)
        tree = first(children)
        @test DT.get_parent(tree) === hierarchy.trees[1]
        @test DT.get_children(tree) ⊢ Set{DT.PolygonTree{Int}}()
        @test DT.get_index(tree) == 2
        @test DT.get_height(tree) == 1
        @test !DT.has_children(tree)
        hierarchy2 = DT.construct_polygon_hierarchy(points, boundary_nodes)

        A, B, C, D = (-9.0, 7.0), (-2.0, 7.0), (-2.0, -1.0), (-9.0, -1.0)
        E, F, G, H, I, J = (2.0, 6.0), (10.0, 6.0), (10.0, -5.0), (0.0, -5.0), (0.0, 5.0), (2.0, 5.0)
        K, L, M, N = (-1.0, 7.0), (1.0, 7.0), (1.0, 6.0), (-1.0, 6.0)
        O, P, Q = (-8.0, 6.0), (-3.0, 4.0), (-3.0, 6.0)
        R, S, T = (-8.0, 5.0), (-3.0, 0.0), (-8.0, 0.0)
        U, V, W, Z = (1.0, 4.0), (8.0, 4.0), (8.0, 3.0), (1.0, 3.0)
        A1, B1, C1 = (1.0, 2.0), (8.0, -4.0), (1.0, -4.0)
        D1, E1, F1, G1 = (-12.0, 12.0), (14.0, 12.0), (14.0, -14.0), (-12.0, -14.0)
        H1, I1, J1, K1 = (-22.0, 12.0), (-14.0, 12.0), (-14.0, 6.0), (-22.0, 6.0)
        L1, M1, N1, O1 = (-22.0, 4.0), (-14.0, 4.0), (-14.0, -14.0), (-22.0, -14.0)
        P1, Q1, R1, S1 = (-20.0, 10.0), (-16.0, 10.0), (-16.0, 8.0), (-20.0, 8.0)
        T1, U1, V1, W1 = (-10.0, -12.0), (-2.0, -12.0), (-2.0, -4.0), (-10.0, -4.0)
        Z1, A2, B2, C2 = (-8.0, -6.0), (-4.0, -6.0), (-4.0, -10.0), (-8.0, -10.0)
        D2, E2, F2, G2 = (-9.0, -5.0), (-9.0, -11.0), (-3.0, -11.0), (-3.0, -5.0)
        H2, I2, J2, K2, L2, M2, N2, O2 = (-9.0, -15.0), (-8.0, -15.0), (-8.0, -16.0), (-3.0, -16.0), (-3.0, -15.0), (-2.0, -15.0), (-2.0, -17.0), (-9.0, -17.0)
        _Q1 = [[A, B, C, D, A]]
        _Q2 = [[E, F, G, H, I, J, E]]
        _Q3 = [[K, L, M, N, K]]
        _Q4 = [[P1, Q1, R1, S1, P1]]
        _Q5 = [[L1, O1, N1, M1, L1]]
        _Q6 = [[H1, K1, J1, I1, H1]]
        _Q7 = [[D1, G1, F1, E1, D1]]
        _Q8 = [[A1, C1, B1, A1]]
        _Q9 = [[U, Z, W, V, U]]
        _Q10 = [[R, T, S, R]]
        _Q11 = [[O, P, Q, O]]
        _Q12 = [[W1, V1, U1, T1, W1]]
        _Q13 = [[Z1, A2, B2, C2, Z1]]
        _Q14 = [[D2, E2, F2, G2, D2]]
        _Q15 = [[H2, O2, N2, M2, L2, K2, J2, I2, H2]]
        boundary_nodes, points = convert_boundary_points_to_indices([_Q1, _Q2, _Q3, _Q4, _Q5, _Q6, _Q7, _Q8, _Q9, _Q10, _Q11, _Q12, _Q13, _Q14, _Q15])
        hierarchy = DT.construct_polygon_hierarchy(points, boundary_nodes)
        @test DT.get_exterior_curve_indices(hierarchy) == Set([6, 5, 7, 15])
        @test DT.get_bounding_boxes(hierarchy) ⊢ [
            DT.BoundingBox(-9.0, -2.0, -1.0, 7.0),
            DT.BoundingBox(0.0, 10.0, -5.0, 6.0),
            DT.BoundingBox(-1.0, 1.0, 6.0, 7.0),
            DT.BoundingBox(-20.0, -16.0, 8.0, 10.0),
            DT.BoundingBox(-22.0, -14.0, -14.0, 4.0),
            DT.BoundingBox(-22.0, -14.0, 6.0, 12.0),
            DT.BoundingBox(-12.0, 14.0, -14.0, 12.0),
            DT.BoundingBox(1.0, 8.0, -4.0, 2.0),
            DT.BoundingBox(1.0, 8.0, 3.0, 4.0),
            DT.BoundingBox(-8.0, -3.0, 0.0, 5.0),
            DT.BoundingBox(-8.0, -3.0, 4.0, 6.0),
            DT.BoundingBox(-10.0, -2.0, -12.0, -4.0),
            DT.BoundingBox(-8.0, -4.0, -10.0, -6.0),
            DT.BoundingBox(-9.0, -3.0, -11.0, -5.0),
            DT.BoundingBox(-9.0, -2.0, -17.0, -15.0)
        ]
        @test DT.get_polygon_orientations(hierarchy) == [
            false, false, false, false, true, true, true, true, true, true, true, false, false, true, true
        ]
        trees = DT.get_trees(hierarchy)
        tree6 = trees[6]
        tree5 = trees[5]
        tree7 = trees[7]
        tree15 = trees[15]
        @test DT.num_children(tree6) == 1 && DT.has_children(tree6)
        @test DT.get_height(tree6) == 0
        @test DT.get_index(tree6) == 6
        @test !DT.has_parent(tree6)
        tree4 = first(DT.get_children(tree6))
        @test DT.num_children(tree4) == 0 && !DT.has_children(tree4)
        @test DT.get_height(tree4) == 1
        @test DT.get_index(tree4) == 4
        @test DT.get_parent(tree4) == tree6
        @test DT.num_children(tree5) == 0 && !DT.has_children(tree5)
        @test DT.get_height(tree5) == 0
        @test DT.get_index(tree5) == 5
        @test !DT.has_parent(tree5)
        @test DT.num_children(tree7) == 4 && DT.has_children(tree7)
        @test DT.get_height(tree7) == 0
        @test DT.get_index(tree7) == 7
        @test !DT.has_parent(tree7)
        tree1, tree3, tree2, tree12 = get_child_from_tree.(Ref(tree7), (1, 3, 2, 12))
        @test DT.get_index(tree1) == 1 && DT.get_height(tree1) == 1 && DT.get_parent(tree1) == tree7
        @test DT.get_index(tree3) == 3 && DT.get_height(tree3) == 1 && DT.get_parent(tree3) == tree7
        @test DT.get_index(tree2) == 2 && DT.get_height(tree2) == 1 && DT.get_parent(tree2) == tree7
        @test DT.get_index(tree12) == 12 && DT.get_height(tree12) == 1 && DT.get_parent(tree12) == tree7
        @test DT.num_children(tree1) == 2 && DT.has_children(tree1)
        tree10, tree11 = get_child_from_tree.(Ref(tree1), (10, 11))
        @test DT.get_index(tree10) == 10 && DT.get_height(tree10) == 2 && DT.get_parent(tree10) == tree1 && !DT.has_children(tree10)
        @test DT.get_index(tree11) == 11 && DT.get_height(tree11) == 2 && DT.get_parent(tree11) == tree1 && !DT.has_children(tree11)
        @test !DT.has_children(tree3)
        @test DT.num_children(tree2) == 2 && DT.has_children(tree2)
        tree8, tree9 = get_child_from_tree.(Ref(tree2), (8, 9))
        @test DT.get_index(tree8) == 8 && DT.get_height(tree8) == 2 && DT.get_parent(tree8) == tree2 && !DT.has_children(tree8)
        @test DT.get_index(tree9) == 9 && DT.get_height(tree9) == 2 && DT.get_parent(tree9) == tree2 && !DT.has_children(tree9)
        @test DT.num_children(tree12) == 1
        tree14 = first(DT.get_children(tree12))
        @test DT.get_index(tree14) == 14 && DT.get_height(tree14) == 2 && DT.get_parent(tree14) == tree12 && DT.has_children(tree14) && DT.num_children(tree14) == 1
        tree13 = first(DT.get_children(tree14))
        @test DT.get_index(tree13) == 13 && DT.get_height(tree13) == 3 && DT.get_parent(tree13) == tree14 && !DT.has_children(tree13)
        @test DT.get_index(tree15) == 15 && DT.get_height(tree15) == 0 && !DT.has_children(tree15) && !DT.has_parent(tree15)
        DT.expand_bounds!(hierarchy)
        @test DT.get_bounding_boxes(hierarchy) ⊢ [
            DT.expand(DT.BoundingBox(-9.0, -2.0, -1.0, 7.0), 0.10),
            DT.expand(DT.BoundingBox(0.0, 10.0, -5.0, 6.0), 0.10),
            DT.expand(DT.BoundingBox(-1.0, 1.0, 6.0, 7.0), 0.10),
            DT.expand(DT.BoundingBox(-20.0, -16.0, 8.0, 10.0), 0.10),
            DT.expand(DT.BoundingBox(-22.0, -14.0, -14.0, 4.0), 0.10),
            DT.expand(DT.BoundingBox(-22.0, -14.0, 6.0, 12.0), 0.10),
            DT.expand(DT.BoundingBox(-12.0, 14.0, -14.0, 12.0), 0.10),
            DT.expand(DT.BoundingBox(1.0, 8.0, -4.0, 2.0), 0.10),
            DT.expand(DT.BoundingBox(1.0, 8.0, 3.0, 4.0), 0.10),
            DT.expand(DT.BoundingBox(-8.0, -3.0, 0.0, 5.0), 0.10),
            DT.expand(DT.BoundingBox(-8.0, -3.0, 4.0, 6.0), 0.10),
            DT.expand(DT.BoundingBox(-10.0, -2.0, -12.0, -4.0), 0.10),
            DT.expand(DT.BoundingBox(-8.0, -4.0, -10.0, -6.0), 0.10),
            DT.expand(DT.BoundingBox(-9.0, -3.0, -11.0, -5.0), 0.10),
            DT.expand(DT.BoundingBox(-9.0, -2.0, -17.0, -15.0), 0.10)
        ]
        @test all([
            DT.BoundingBox(-9.0, -2.0, -1.0, 7.0),
            DT.BoundingBox(0.0, 10.0, -5.0, 6.0),
            DT.BoundingBox(-1.0, 1.0, 6.0, 7.0),
            DT.BoundingBox(-20.0, -16.0, 8.0, 10.0),
            DT.BoundingBox(-22.0, -14.0, -14.0, 4.0),
            DT.BoundingBox(-22.0, -14.0, 6.0, 12.0),
            DT.BoundingBox(-12.0, 14.0, -14.0, 12.0),
            DT.BoundingBox(1.0, 8.0, -4.0, 2.0),
            DT.BoundingBox(1.0, 8.0, 3.0, 4.0),
            DT.BoundingBox(-8.0, -3.0, 0.0, 5.0),
            DT.BoundingBox(-8.0, -3.0, 4.0, 6.0),
            DT.BoundingBox(-10.0, -2.0, -12.0, -4.0),
            DT.BoundingBox(-8.0, -4.0, -10.0, -6.0),
            DT.BoundingBox(-9.0, -3.0, -11.0, -5.0),
            DT.BoundingBox(-9.0, -2.0, -17.0, -15.0)
        ] .∈ DT.get_bounding_boxes(hierarchy))
        @test traverse_tree(hierarchy.trees[7]) ≠ traverse_tree(hierarchy.trees[15])
        @test compare_trees(hierarchy, hierarchy)
        @test compare_trees(deepcopy(hierarchy), deepcopy(hierarchy))
        @test !compare_trees(hierarchy, hierarchy2)
        @test hierarchy.trees[7] ≠ hierarchy.trees[15]
        @test hierarchy == hierarchy
        @test deepcopy(hierarchy) == hierarchy
        @test deepcopy(hierarchy) == deepcopy(hierarchy)
        @test hierarchy ≠ hierarchy2
        @test deepcopy(hierarchy) ≠ hierarchy2
    end

    @testset "Curve-bounded" begin
        IntegerType = Int
        points = deepcopy(points_IV)
        boundary_nodes = deepcopy(curve_IV)
        boundary_curves, nnew_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, IntegerType)
        orig_points = deepcopy(points)
        orig_boundary_nodes = deepcopy(boundary_nodes)
        new_points, new_boundary_nodes = DT.polygonise(points, nnew_boundary_nodes, boundary_curves)
        @test points == orig_points # check that the original points are not modified
        @test boundary_nodes == orig_boundary_nodes # check that the original boundary nodes are not modified
        hierarchy = DT.construct_polygon_hierarchy(points, nnew_boundary_nodes, boundary_curves)
        DT.expand_bounds!(hierarchy, DT.ε)
        @test DT.get_bounding_boxes(hierarchy) ⊢ [DT.BoundingBox(-1 - 2DT.ε, 1 + 2DT.ε, -1 - 2DT.ε, 1 + 2DT.ε)] &&
              DT.get_polygon_orientations(hierarchy) ⊢ BitVector([1])
        trees = DT.get_trees(hierarchy)
        @test length(trees) == 1
        tree = trees[1]
        @test DT.get_index(tree) == 1 && DT.get_height(tree) == 0 && !DT.has_parent(tree) && !DT.has_children(tree)

        points = deepcopy(points_V)
        boundary_nodes = deepcopy(curve_V)
        boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, IntegerType)
        hierarchy = DT.construct_polygon_hierarchy(points, new_boundary_nodes, boundary_curves)
        @test DT.get_polygon_orientations(hierarchy) ⊢ BitVector([1])
        trees = DT.get_trees(hierarchy)
        @test length(trees) == 1
        tree = trees[1]
        @test DT.get_index(tree) == 1 && DT.get_height(tree) == 0 && !DT.has_parent(tree) && !DT.has_children(tree)

        points = deepcopy(points_VI)
        boundary_nodes = deepcopy(curve_VI)
        boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, IntegerType)
        hierarchy = DT.construct_polygon_hierarchy(points, new_boundary_nodes, boundary_curves)
        @test DT.get_polygon_orientations(hierarchy) ⊢ BitVector([1])
        trees = DT.get_trees(hierarchy)
        @test length(trees) == 1
        tree = trees[1]
        @test DT.get_index(tree) == 1 && DT.get_height(tree) == 0 && !DT.has_parent(tree) && !DT.has_children(tree)

        points = deepcopy(points_VII)
        boundary_nodes = deepcopy(curve_VII)
        boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, IntegerType)
        hierarchy = DT.construct_polygon_hierarchy(points, new_boundary_nodes, boundary_curves)
        @test DT.get_polygon_orientations(hierarchy) ⊢ BitVector([1])
        trees = DT.get_trees(hierarchy)
        @test length(trees) == 1
        tree = trees[1]
        @test DT.get_index(tree) == 1 && DT.get_height(tree) == 0 && !DT.has_parent(tree) && !DT.has_children(tree)

        points = deepcopy(points_VIII)
        boundary_nodes = deepcopy(curve_VIII)
        boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, IntegerType)
        hierarchy = DT.construct_polygon_hierarchy(points, new_boundary_nodes, boundary_curves)
        @test DT.get_polygon_orientations(hierarchy) ⊢ BitVector([1])
        trees = DT.get_trees(hierarchy)
        @test length(trees) == 1
        tree = trees[1]
        @test DT.get_index(tree) == 1 && DT.get_height(tree) == 0 && !DT.has_parent(tree) && !DT.has_children(tree)

        points = deepcopy(points_IX)
        boundary_nodes = deepcopy(curve_IX)
        boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, IntegerType)
        hierarchy = DT.construct_polygon_hierarchy(points, new_boundary_nodes, boundary_curves)
        @test DT.get_polygon_orientations(hierarchy) ⊢ BitVector([1, 0])
        trees = DT.get_trees(hierarchy)
        @test length(trees) == 1
        tree = trees[1]
        @test DT.get_index(tree) == 1 && DT.get_height(tree) == 0 && !DT.has_parent(tree) && DT.has_children(tree) && DT.num_children(tree) == 1
        tree = first(DT.get_children(tree))
        @test DT.get_index(tree) == 2 && DT.get_height(tree) == 1 && DT.get_parent(tree) == trees[1] && !DT.has_children(tree)

        points = deepcopy(points_X)
        boundary_nodes = deepcopy(curve_X)
        boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, IntegerType)
        hierarchy = DT.construct_polygon_hierarchy(points, new_boundary_nodes, boundary_curves)
        @test DT.get_polygon_orientations(hierarchy) ⊢ BitVector([1, 0])
        trees = DT.get_trees(hierarchy)
        @test length(trees) == 1
        tree = trees[1]
        @test DT.get_index(tree) == 1 && DT.get_height(tree) == 0 && !DT.has_parent(tree) && DT.has_children(tree) && DT.num_children(tree) == 1
        tree = first(DT.get_children(tree))
        @test DT.get_index(tree) == 2 && DT.get_height(tree) == 1 && DT.get_parent(tree) == trees[1] && !DT.has_children(tree)

        points = deepcopy(points_XI)
        boundary_nodes = deepcopy(curve_XI)
        boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(points, boundary_nodes, IntegerType)
        hierarchy = DT.construct_polygon_hierarchy(points, new_boundary_nodes, boundary_curves)
        @test DT.get_polygon_orientations(hierarchy) ⊢ BitVector([1, 0, 1, 1, 0, 0])
        trees = DT.get_trees(hierarchy)
        @test length(trees) == 2
        tree1 = trees[1]
        @test DT.get_index(tree1) == 1 && DT.get_height(tree1) == 0 && !DT.has_parent(tree1) && DT.has_children(tree1) && DT.num_children(tree1) == 1
        tree2 = first(DT.get_children(tree1))
        @test DT.get_index(tree2) == 2 && DT.get_height(tree2) == 1 && DT.get_parent(tree2) == tree1 && !DT.has_children(tree2)
        tree3 = trees[3]
        @test DT.get_index(tree3) == 3 && DT.get_height(tree3) == 0 && !DT.has_parent(tree3) && DT.has_children(tree3) && DT.num_children(tree3) == 1
        tree6 = first(DT.get_children(tree3))
        @test DT.get_index(tree6) == 6 && DT.get_height(tree6) == 1 && DT.get_parent(tree6) == tree3 && DT.has_children(tree6)
        tree4 = first(DT.get_children(tree6))
        @test DT.get_index(tree4) == 4 && DT.get_height(tree4) == 2 && DT.get_parent(tree4) == tree6 && DT.has_children(tree4)
        tree5 = first(DT.get_children(tree4))
        @test DT.get_index(tree5) == 5 && DT.get_height(tree5) == 3 && DT.get_parent(tree5) == tree4 && !DT.has_children(tree5)
    end

    @testset "Shared endpoint" begin
        points = [(-1.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 2.0), (-1.0, 2.0)]
        boundary_nodes = [[[1, 2, 3, 1]], [[3, 4, 5, 3]]]
        hierarchy = DT.construct_polygon_hierarchy(points, boundary_nodes)
        trees = DT.get_trees(hierarchy)
        @test length(trees) == 2
        tree1 = trees[1]
        @test DT.get_index(tree1) == 1 && DT.get_height(tree1) == 0 && !DT.has_parent(tree1) && !DT.has_children(tree1)
        tree2 = trees[2]
        @test DT.get_index(tree2) == 2 && DT.get_height(tree2) == 0 && !DT.has_parent(tree2) && !DT.has_children(tree2)
    end
end