using ..DelaunayTriangulation
using ..DelaunayTriangulation: EllipticalArc
const DT = DelaunayTriangulation

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
        [CircularArc((1.1, -3.0), (1.1, -3.0), (0.0, -3.0), positive = false)],
    ],
]
points = [(-2.0, 0.0), (0.0, 0.0), (2.0, 0.0), (-2.0, -5.0), (2.0, -5.0), (2.0, -1 / 10), (-2.0, -1 / 10), (-1.0, -3.0), (0.0, -4.0), (0.0, -2.3), (-0.5, -3.5), (0.9, -3.0)]
tri = triangulate(points; boundary_nodes = curve)
refine!(tri; max_area = 1.0e-2)

q = (-1.0, -4.0)
idx = find_polygon(tri, q)
@test idx == 3
q = (-1.0, 0.2)
idx = find_polygon(tri, q)
@test idx == 1
q = (-1 / 2, -3.0)
idx = find_polygon(tri, q)
@test idx == 4
q = (0.0, -4.05)
idx = find_polygon(tri, q)
@test idx == 6
q = (0.0, -5.0)
idx = find_polygon(tri, q)
@test idx == 0
