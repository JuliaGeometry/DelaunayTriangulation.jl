using DelaunayTriangulation, CairoMakie, ReferenceTests, StableRNGs

# Unconstrained 
rng = StableRNG(123)
points = rand(rng, 2, 50)
tri1 = triangulate(points; rng) # default predicate kernel is AdaptiveKernel()

# Voronoi example 
vorn2 = voronoi(tri1; rng)

# Clipped Voronoi 
vorn3 = voronoi(tri1, clip = true; rng, predicates = ExactKernel()) # you can change the predicate kernel

# Smoothed Voronoi 
vorn4 = centroidal_smooth(vorn3; rng, predicates = FastKernel())

# Constrained example with refinement 
boundary_points = [
    (0.0, 0.0), (1.0, 0.0), (1.0, 0.3), (0.5, 0.3),
    (0.3, 0.7), (0.1, 1.0), (0.0, 1.0), (0.0, 0.0),
]
boundary_nodes, points = convert_boundary_points_to_indices(boundary_points)
tri5 = triangulate(points; boundary_nodes, rng)
refine!(tri5; max_area = 1.0e-2get_area(tri5), rng)

# Disjoint constrained example with refinement 
boundary_points = [
    [[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]],
    [[(0.3, 0.3), (0.3, 0.7), (0.7, 0.7), (0.7, 0.3), (0.3, 0.3)]],
    [
        [
            (1.2, 0.0), (1.4, 0.0), (1.4, 1.2), (0.0, 1.2), (0.0, 1.1),
            (1.1, 1.1), (1.1, 0.0), (1.2, 0.0),
        ],
    ],
]
boundary_nodes, points = convert_boundary_points_to_indices(boundary_points)
tri6 = triangulate(points; boundary_nodes, rng)
refine!(tri6; max_area = 1.0e-2get_area(tri6), rng)

# Curve-bounded example
using DelaunayTriangulation: EllipticalArc
ellipse = EllipticalArc((1.0, 0.0), (1.0, 0.0), (0.0, 0.0), 1.0, 2.0, 0.0)
tri7 = triangulate(NTuple{2, Float64}[]; boundary_nodes = [ellipse], rng)
refine!(tri7; max_area = 1.0e-2get_area(tri7), rng)

# Disjoint curve-bounded example 
ellipse = EllipticalArc((7.0, 3.5), (7.0, 3.5), (0.0, 3.5), 7.0, 10.0, 0.0)
catmull_cp = [
    (0.0, 0.0), (-2.0, -1.0), (-4.0, 0.0), (-5.0, 2.0), (-1.0, 4.0), (0.0, 3.0),
    (1.0, 4.0), (5.0, 2.0), (4.0, 0.0), (2.0, -1.0), (0.0, 0.0),
]
catmull_spl = CatmullRomSpline(catmull_cp)
circle = CircularArc((0.5, 1.5), (0.5, 1.5), (0.0, 1.0))
circle2 = CircularArc((0.1, 1.5), (0.1, 1.5), (0.0, 1.0), positive = false)
points = [(-4.0, -10.0), (4.0, -10.0), (4.0, -7.0), (-4.0, -7.0)]
square = [1, 2, 3, 4, 1]
boundary_nodes = [[square], [[ellipse]], [[catmull_spl]], [[circle]], [[circle2]]]
tri8 = triangulate(points; boundary_nodes, rng)
refine!(tri8; max_area = 1.0e-2get_area(tri8), rng) # could also use find_polygon to help define a custom refinement function for each shape

# Plotting 
fig = Figure(fontsize = 42, size = (2800, 1480))
ax = Axis(fig[1, 1], title = "Unconstrained", width = 600, height = 600);            triplot!(ax, tri1)
ax = Axis(fig[1, 2], title = "Voronoi", width = 600, height = 600);                  voronoiplot!(ax, vorn2)
ax = Axis(fig[1, 3], title = "Clipped Voronoi", width = 600, height = 600);          voronoiplot!(ax, vorn3)
ax = Axis(fig[1, 4], title = "Centroidal Voronoi", width = 600, height = 600);       voronoiplot!(ax, vorn4)
ax = Axis(fig[2, 1], title = "Constrained", width = 600, height = 600);              triplot!(ax, tri5)
ax = Axis(fig[2, 2], title = "Disjoint Constrained", width = 600, height = 600);     triplot!(ax, tri6)
ax = Axis(fig[2, 3], title = "Curve-Bounded", width = 600, height = 600);            triplot!(ax, tri7)
ax = Axis(fig[2, 4], title = "Disjoint Curve-Bounded", width = 600, height = 600);   triplot!(ax, tri8)

readme_img = joinpath(dirname(dirname(pathof(DelaunayTriangulation))), "readme.png")
@test_reference readme_img fig by = psnr_equality(10)
