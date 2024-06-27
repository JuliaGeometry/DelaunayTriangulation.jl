using DelaunayTriangulation, CairoMakie

## Domain definition 
θ = 5π / 64
cs = θ -> (cos(θ), sin(θ))
p₁, q₁ = cs(-π / 2 - θ), cs(θ) # Absorbing 
p₂, q₂ = q₁, cs(π / 2 - θ)     # Reflecting 
p₃, q₃ = q₂, cs(π + θ)         # Absorbing 
p₄, q₄ = q₃, p₁                # Reflecting
c₀ = (0.0, 0.0)
𝒞₀₁ = CircularArc(p₁, q₁, c₀) # first, last, center
𝒞₀₂ = CircularArc(p₂, q₂, c₀)
𝒞₀₃ = CircularArc(p₃, q₃, c₀)
𝒞₀₄ = CircularArc(p₄, q₄, c₀)
c₁, p₅ = (-0.4, -0.4), (-0.65, -0.65)
c₂, p₆ = (0.4, 0.4), (0.65, 0.65)
𝒞₁ = CircularArc(p₅, p₅, c₁, positive=false) # Reflecting
𝒞₂ = CircularArc(p₆, p₆, c₂, positive=false) # Reflecting
sink = (0.0, 0.0)
fig = Figure()
ax = Axis(fig[1, 1])
t = LinRange(0, 1, 1000)
[lines!(ax, 𝒞.(t), color=:red, linewidth=4) for 𝒞 in (𝒞₀₁, 𝒞₀₃)]
[lines!(ax, 𝒞.(t), color=:blue, linewidth=4) for 𝒞 in (𝒞₀₂, 𝒞₀₄, 𝒞₁, 𝒞₂)]
scatter!(ax, [sink], color=:black, markersize=14)
text!(ax, [(0.0, -0.2),
        (-0.7, -0.5), (0.8, 0.6),
        (0.7, -0.5), (-0.7, 0.5)],
    text=[L"(x_s, y_s)",
        L"\Gamma_r", L"\Gamma_r",
        L"\Gamma_a", L"\Gamma_a"],
    fontsize=26)
fig
save("paper/figure0.png", fig)

## Example 1: Mean exit time
# The outer circle
θ = 5π / 64
cs = θ -> (cos(θ), sin(θ))
p₁, q₁ = cs(-π / 2 - θ), cs(θ) # Absorbing 
p₂, q₂ = q₁, cs(π / 2 - θ)     # Reflecting 
p₃, q₃ = q₂, cs(π + θ)         # Absorbing 
p₄, q₄ = q₃, p₁                # Reflecting
c₀ = (0.0, 0.0)
𝒞₀₁ = CircularArc(p₁, q₁, c₀) # first, last, center
𝒞₀₂ = CircularArc(p₂, q₂, c₀)
𝒞₀₃ = CircularArc(p₃, q₃, c₀)
𝒞₀₄ = CircularArc(p₄, q₄, c₀)
𝒞₀ = [[𝒞₀₁], [𝒞₀₂], [𝒞₀₃], [𝒞₀₄]]
# Inner circles
c₁, p₅ = (-0.4, -0.4), (-0.65, -0.65)
c₂, p₆ = (0.4, 0.4), (0.65, 0.65)
𝒞₁ = CircularArc(p₅, p₅, c₁, positive=false) # Reflecting
𝒞₂ = CircularArc(p₆, p₆, c₂, positive=false) # Reflecting
# Triangulate and refine
sink = (0.0, 0.0)
tri = triangulate([sink], boundary_nodes=[𝒞₀, [[𝒞₁]], [[𝒞₂]]])
refine!(tri; max_area=1e-3get_area(tri))
# Plot
fig, ax, sc = triplot(tri,
    axis=(width=500, height=500, title="(a): Triangulation"),
    figure=(fontsize=30,))
colors = [:red, :blue, :darkgreen, :magenta, :black, :purple]
section = (tri, i) -> [get_point(tri, i) for i in get_neighbours(tri, i)]
for i in each_ghost_vertex(tri)
    scatter!(ax, section(tri, i), color=colors[-i], markersize=15) # ghost vertices are negative
end
fig

# Solve
using FiniteVolumeMethod, OrdinaryDiffEq, SteadyStateDiffEq
D, zerof = 6.25e-4, Returns(0.0)
geo = FVMGeometry(tri)
zero_bc = ntuple(i -> zerof, Val(num_ghost_vertices(tri))) # all the boundary conditions are homogeneous
bc_types = (Dirichlet, Neumann, Dirichlet, Neumann, Neumann, Neumann)
BCs = BoundaryConditions(geo, zero_bc, bc_types)
ICs = InternalConditions(zerof, dirichlet_nodes=Dict(1 => findfirst(==(sink), get_points(tri))))
initial_guess = map(DelaunayTriangulation.each_point(tri)) do p
    p == sink && return 0.0
    inside = DelaunayTriangulation.dist(tri, p) > 0.0
    return !inside ? 0.0 : (1 - sqrt(getx(p)^2 + gety(p)^2)) / (4D) # (R^2 - r^2)/(4D) is solution in a disk of radius R
end
prob = FVMProblem(geo, BCs, ICs;
    diffusion_function=Returns(D), source_function=Returns(1.0),
    initial_condition=initial_guess, final_time=Inf)
met_prob = SteadyFVMProblem(prob)
sol = solve(met_prob, DynamicSS(Rosenbrock23()))
ax = Axis(fig[1, 2], width=500, height=500, title="(b): Mean exit time")
hm = tricontourf!(ax, tri, sol.u, levels=0:10:250, extendhigh=:auto)
poly!(ax, 𝒞₁.(LinRange(0, 1, 2500)), color=:white)
poly!(ax, 𝒞₂.(LinRange(0, 1, 2500)), color=:white)
ax.yticklabelsvisible = false
Colorbar(fig[1,3], hm, label = L"T(x, y)")
resize_to_layout!(fig)
save("paper/figure1.png", fig)

## Example 2: K-means clustering (assuming a dense data set)
using Random
Random.seed!(123)
k = 7
clusters = [(rand(), rand()) for _ in 1:k]
# Assume data live in [0, 1]²
push!(clusters, (0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0))
# Tessellate and smooth 
tri = triangulate(clusters)
vor = voronoi(tri, clip=true) # don't want unbounded
cvor = centroidal_smooth(vor)
# Generate data and assign 
data = rand(2, 2500)
label = p -> get_nearest_neighbour(cvor, p)
labels = label.(eachcol(data))
fig, ax, sc = voronoiplot(cvor, color=:white, strokewidth=4,
    show_generators=false)
colors = cgrad(:jet, num_polygons(cvor), categorical=true)
scatter!(ax, data, color=colors[labels])
fig
save("paper/figure2.png", fig)