# Benchmarks for DelaunayTriangulation.jl, as a BenchmarkTools `SUITE`. See CONTRIBUTING.md for how to run them and how
# to compare two versions of the package. All inputs are generated from fixed seeds when this file is included, and any
# input that a benchmark mutates is copied in the setup of each sample, so that only the operation itself is timed. Only
# the public API is used, so that older versions of the package can be benchmarked with the same file.

using BenchmarkTools
using DelaunayTriangulation
using StableRNGs
const DT = DelaunayTriangulation

const SUITE = BenchmarkGroup()

## Inputs
function uniform_points(n; seed = n, a = 0.0, b = 1.0)
    rng = StableRNG(seed)
    return [(a + (b - a) * rand(rng), a + (b - a) * rand(rng)) for _ in 1:n]
end

# Points in 20 Gaussian clusters, with spreads from 1e-4 to 1e-1
function clustered_points(n; seed = n)
    rng = StableRNG(seed)
    centres = [(rand(rng), rand(rng)) for _ in 1:20]
    spreads = [10.0^(-rand(rng, 1:4)) for _ in 1:20]
    points = NTuple{2, Float64}[]
    for _ in 1:n
        c = rand(rng, 1:20)
        push!(points, (centres[c][1] + spreads[c] * randn(rng), centres[c][2] + spreads[c] * randn(rng)))
    end
    return points
end

# n points on a circle, so that all of them are on the convex hull
circle_points(n; centre = (0.0, 0.0), radius = 1.0) = [(centre[1] + radius * cos(θ), centre[2] + radius * sin(θ)) for θ in LinRange(0, 2π, n + 1)[1:(end - 1)]]

# A closed boundary through n points on a circle, counter-clockwise unless clockwise = true
function circle_boundary(n, centre, radius; clockwise = false)
    points = circle_points(n; centre, radius)
    clockwise && reverse!(points)
    push!(points, points[begin])
    return points
end

# A closed counter-clockwise boundary of the unit square with n points on each side
function square_boundary(n)
    t = LinRange(0, 1, n + 1)[1:(end - 1)]
    points = [[(s, 0.0) for s in t]; [(1.0, s) for s in t]; [(1 - s, 1.0) for s in t]; [(0.0, 1 - s) for s in t]]
    push!(points, points[begin])
    return points
end

# The unit square with n boundary points on each side and m random interior points
function square_domain(n, m; seed = 1)
    boundary_nodes, points = convert_boundary_points_to_indices(square_boundary(n))
    append!(points, uniform_points(m; seed, a = 0.01, b = 0.99))
    return points, boundary_nodes
end

# The unit square with four circular holes, with n points on each boundary curve and m random points, some of which are in the holes
function holes_domain(n, m; seed = 1)
    curves = [[square_boundary(n)]]
    for centre in ((0.25, 0.25), (0.75, 0.25), (0.25, 0.75), (0.75, 0.75))
        push!(curves, [circle_boundary(n, centre, 0.15; clockwise = true)])
    end
    boundary_nodes, points = convert_boundary_points_to_indices(curves)
    append!(points, uniform_points(m; seed, a = 0.01, b = 0.99))
    return points, boundary_nodes
end

# Three disjoint discs with n points on each boundary and m random points in each disc
function disjoint_domain(n, m; seed = 1)
    rng = StableRNG(seed)
    centres = ((0.0, 0.0), (3.0, 0.0), (0.0, 3.0))
    boundary_nodes, points = convert_boundary_points_to_indices([[circle_boundary(n, centre, 1.0)] for centre in centres])
    for centre in centres, _ in 1:m
        r, θ = 0.95sqrt(rand(rng)), 2π * rand(rng)
        push!(points, (centre[1] + r * cos(θ), centre[2] + r * sin(θ)))
    end
    return points, boundary_nodes
end

# m random points together with k horizontal lines of segments, each line through 50 points
function segments_input(m, k; seed = 1)
    points = uniform_points(m; seed)
    segments = Set{NTuple{2, Int}}()
    for y in LinRange(0.05, 0.95, k)
        n = length(points)
        append!(points, [(x, y) for x in LinRange(0.05, 0.95, 50)])
        foreach(i -> push!(segments, (n + i, n + i + 1)), 1:49)
    end
    return points, segments
end

# The domain from the curve-bounded tutorial, with an elliptical arc, a B-spline, a Bézier curve, a Catmull-Rom spline, a circular arc, and piecewise linear sections
function curve_domain()
    curve = [
        [[1, 2, 3], [EllipticalArc((2.0, 0.0), (-2.0, 0.0), (0.0, 0.0), 2, 1 / 2, 0.0)]],
        [[BSpline([(0.0, 0.4), (1.0, 0.2), (0.0, 0.1), (-1.0, 0.2), (0.0, 0.4)])]],
        [[4, 5, 6, 7, 4]],
        [[BezierCurve([(0.0, -2.0), (0.0, -2.5), (-1.0, -2.5), (-1.0, -3.0)])], [CatmullRomSpline([(-1.0, -3.0), (0.0, -4.0), (1.0, -3.0), (0.0, -2.0)])]],
        [[12, 11, 10, 12]],
        [[CircularArc((1.1, -3.0), (1.1, -3.0), (0.0, -3.0), positive = false)]],
    ]
    points = [(-2.0, 0.0), (0.0, 0.0), (2.0, 0.0), (-2.0, -5.0), (2.0, -5.0), (2.0, -1 / 10), (-2.0, -1 / 10), (-1.0, -3.0), (0.0, -4.0), (0.0, -2.3), (-0.5, -3.5), (0.9, -3.0)]
    return points, curve
end

# Some operations, applied to each of a collection of inputs
function find_triangles(tri, queries, rng)
    for q in queries
        find_triangle(tri, q; rng)
    end
    return tri
end
function find_polygons(tri, queries)
    for q in queries
        find_polygon(tri, q)
    end
    return tri
end
function get_nearest_neighbours(tri, queries, rng)
    for q in queries
        get_nearest_neighbour(tri, q; rng)
    end
    return tri
end
function add_points!(tri, points, rng)
    for p in points
        add_point!(tri, p; rng)
    end
    return tri
end
function delete_points!(tri, vertices, rng)
    for v in vertices
        delete_point!(tri, v; rng)
    end
    return tri
end
function add_segments!(tri, segments, rng)
    for (i, j) in segments
        add_segment!(tri, i, j; rng)
    end
    return tri
end

## Inputs that are used by several benchmarks
const POINTS_1K = uniform_points(1_000)
const POINTS_10K = uniform_points(10_000)
const POINTS_100K = uniform_points(100_000)
const TRI_10K = triangulate(POINTS_10K; rng = StableRNG(1))
const TRI_100K = triangulate(POINTS_100K; rng = StableRNG(1))
const HOLES_POINTS, HOLES_BOUNDARY_NODES = holes_domain(100, 10_000)
const TRI_HOLES = triangulate(HOLES_POINTS; boundary_nodes = HOLES_BOUNDARY_NODES, rng = StableRNG(1))
const QUERIES = uniform_points(1_000; seed = 2)

## Triangulating
# The segments can get rotated or split, and the indices in S can get shifted, so these are copied
SUITE["triangulate"] = BenchmarkGroup()
for n in (1_000, 10_000, 100_000)
    points = uniform_points(n)
    SUITE["triangulate"]["uniform, $n points"] = @benchmarkable triangulate($points; rng = StableRNG(1))
end
let points = clustered_points(10_000)
    SUITE["triangulate"]["clustered, 10000 points"] = @benchmarkable triangulate($points; rng = StableRNG(1))
end
let points = get_points(triangulate_rectangle(0, 1, 0, 1, 100, 100))
    SUITE["triangulate"]["grid, 100 × 100 points"] = @benchmarkable triangulate($points; rng = StableRNG(1))
end
# Weights of about the squared distance between neighbouring points, so that only some of the points are submerged
let weights = 1.0e-5randn(StableRNG(3), 10_000)
    SUITE["triangulate"]["weighted, 10000 points"] = @benchmarkable triangulate($POINTS_10K; weights = $weights, rng = StableRNG(1))
end
let (points, segments) = segments_input(10_000, 20)
    SUITE["triangulate"]["segments, 10000 points and 20 × 49 segments"] = @benchmarkable triangulate($points; segments = _segments, rng = StableRNG(1)) setup = (_segments = copy($segments))
end
let (points, boundary_nodes) = square_domain(100, 10_000)
    SUITE["triangulate"]["square boundary, 10400 points"] = @benchmarkable triangulate($points; boundary_nodes = $boundary_nodes, rng = StableRNG(1))
end
SUITE["triangulate"]["multiply connected, 4 holes, 12000 points"] = @benchmarkable triangulate($HOLES_POINTS; boundary_nodes = $HOLES_BOUNDARY_NODES, rng = StableRNG(1))
let (points, boundary_nodes) = disjoint_domain(100, 3_000)
    SUITE["triangulate"]["disjoint, 3 discs, 9300 points"] = @benchmarkable triangulate($points; boundary_nodes = $boundary_nodes, rng = StableRNG(1))
end
let (points, curve) = curve_domain()
    # The points get mutated for curve-bounded domains, so they are copied
    SUITE["triangulate"]["curve-bounded"] = @benchmarkable triangulate(_points; boundary_nodes = $curve, rng = StableRNG(1)) setup = (_points = copy($points))
end
SUITE["triangulate"]["triangulate_rectangle, 300 × 300"] = @benchmarkable triangulate_rectangle(0, 1, 0, 1, 300, 300)
let points = circle_points(1_000)
    S = collect(eachindex(points))
    SUITE["triangulate"]["triangulate_convex, 1000 vertices"] = @benchmarkable triangulate_convex($points, _S; rng = StableRNG(1)) setup = (_S = copy($S))
end

## Refinement
# refine! adds points to the points, boundary nodes and segments of the triangulation, which are the inputs, so these are copied
SUITE["refine"] = BenchmarkGroup()
let (points, boundary_nodes) = square_domain(10, 0)
    SUITE["refine"]["square, max_area = 1e-4"] = @benchmarkable refine!(tri; max_area = 1.0e-4, rng = StableRNG(1)) setup = (tri = triangulate(copy($points); boundary_nodes = deepcopy($boundary_nodes), rng = StableRNG(1)))
    SUITE["refine"]["square, min_angle = 33, max_area = 1e-3"] = @benchmarkable refine!(tri; min_angle = 33.0, max_area = 1.0e-3, rng = StableRNG(1)) setup = (tri = triangulate(copy($points); boundary_nodes = deepcopy($boundary_nodes), rng = StableRNG(1)))
end
let (points, boundary_nodes) = holes_domain(20, 0)
    SUITE["refine"]["multiply connected, max_area = 1e-4"] = @benchmarkable refine!(tri; max_area = 1.0e-4, rng = StableRNG(1)) setup = (tri = triangulate(copy($points); boundary_nodes = deepcopy($boundary_nodes), rng = StableRNG(1)))
end
let (points, segments) = segments_input(1_000, 5)
    SUITE["refine"]["segments, max_area = 1e-4"] = @benchmarkable refine!(tri; max_area = 1.0e-4, rng = StableRNG(1)) setup = (tri = triangulate(copy($points); segments = copy($segments), rng = StableRNG(1)))
end
let (points, curve) = curve_domain()
    SUITE["refine"]["curve-bounded, max_area = 1e-3"] = @benchmarkable refine!(tri; max_area = 1.0e-3, rng = StableRNG(1)) setup = (tri = triangulate(copy($points); boundary_nodes = $curve, rng = StableRNG(1)))
end

## Operations on an existing triangulation
SUITE["operations"] = BenchmarkGroup()
let points = uniform_points(1_000; seed = 4)
    SUITE["operations"]["add_point!, 1000 points into 10000"] = @benchmarkable add_points!(tri, $points, StableRNG(1)) setup = (tri = copy($TRI_10K))
end
let vertices = filter(v -> !DT.is_boundary_node(TRI_10K, v)[1], 1:2:2_000)
    SUITE["operations"]["delete_point!, 1000 of 10000 points"] = @benchmarkable delete_points!(tri, $vertices, StableRNG(1)) setup = (tri = copy($TRI_10K))
end
let (points, segments) = segments_input(2_000, 2)
    tri = triangulate(points; rng = StableRNG(1))
    SUITE["operations"]["add_segment!, 98 segments"] = @benchmarkable add_segments!(_tri, $(collect(segments)), StableRNG(1)) setup = (_tri = copy($tri))
end
SUITE["operations"]["lock_convex_hull!, 10000 points"] = @benchmarkable lock_convex_hull!(tri; rng = StableRNG(1)) setup = (tri = copy($TRI_10K))
SUITE["operations"]["copy, 10000 points"] = @benchmarkable copy($TRI_10K)
SUITE["operations"]["each_solid_triangle, 100000 points"] = @benchmarkable count(_ -> true, each_solid_triangle($TRI_100K))
SUITE["operations"]["each_solid_edge, 100000 points"] = @benchmarkable count(_ -> true, each_solid_edge($TRI_100K))

## Point location
SUITE["point location"] = BenchmarkGroup()
SUITE["point location"]["find_triangle, 1000 queries in 10000 points"] = @benchmarkable find_triangles($TRI_10K, $QUERIES, StableRNG(1))
let queries = uniform_points(1_000; seed = 5, a = -0.5, b = 1.5)
    SUITE["point location"]["find_triangle, 1000 queries in and outside 10000 points"] = @benchmarkable find_triangles($TRI_10K, $queries, StableRNG(1))
end
SUITE["point location"]["find_triangle, 1000 queries in a multiply connected domain"] = @benchmarkable find_triangles($TRI_HOLES, $QUERIES, StableRNG(1))
SUITE["point location"]["find_polygon, 1000 queries in a multiply connected domain"] = @benchmarkable find_polygons($TRI_HOLES, $QUERIES)
SUITE["point location"]["get_nearest_neighbour, 1000 queries in 10000 points"] = @benchmarkable get_nearest_neighbours($TRI_10K, $QUERIES, StableRNG(1))

## Voronoi tessellations and power diagrams
SUITE["voronoi"] = BenchmarkGroup()
SUITE["voronoi"]["unclipped, 10000 points"] = @benchmarkable voronoi(tri; rng = StableRNG(1)) setup = (tri = copy($TRI_10K))
SUITE["voronoi"]["clipped, 10000 points"] = @benchmarkable voronoi(tri; clip = true, rng = StableRNG(1)) setup = (tri = copy($TRI_10K))
let clip_points = circle_points(50; centre = (0.5, 0.5), radius = 0.4)
    clip_polygon = (clip_points, [collect(eachindex(clip_points)); 1])
    SUITE["voronoi"]["clipped to a polygon, 10000 points"] = @benchmarkable voronoi(tri; clip = true, clip_polygon = $clip_polygon, rng = StableRNG(1)) setup = (tri = copy($TRI_10K))
end
let tri = triangulate(uniform_points(500); rng = StableRNG(1))
    SUITE["voronoi"]["centroidal, 500 points"] = @benchmarkable voronoi(_tri; clip = true, smooth = true, rng = StableRNG(1)) setup = (_tri = copy($tri))
end
let tri = triangulate(POINTS_10K; weights = 1.0e-5randn(StableRNG(3), 10_000), rng = StableRNG(1))
    SUITE["voronoi"]["power diagram, 10000 points"] = @benchmarkable voronoi(_tri; rng = StableRNG(1)) setup = (_tri = copy($tri))
end

## Geometry
SUITE["geometry"] = BenchmarkGroup()
SUITE["geometry"]["convex_hull, 100000 points"] = @benchmarkable convex_hull($POINTS_100K)
let points = circle_points(10_000)
    SUITE["geometry"]["convex_hull, 10000 points on a circle"] = @benchmarkable convex_hull($points)
end
SUITE["geometry"]["pole_of_inaccessibility, multiply connected domain"] = @benchmarkable DT.pole_of_inaccessibility($HOLES_POINTS, $HOLES_BOUNDARY_NODES)
SUITE["geometry"]["statistics, 10000 points"] = @benchmarkable statistics($TRI_10K)
SUITE["geometry"]["get_area, multiply connected domain"] = @benchmarkable get_area($TRI_HOLES)
