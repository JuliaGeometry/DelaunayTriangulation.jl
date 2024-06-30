# # Cellular Biology 
# Our next application concerns the simulation of cell dynamics in two dimensions.[^1] We only consider a very basic 
# model here, allowing only diffusion and proliferation. The point here is to just demonstrate how we can make use of the neighbourhood graph 
# in a `Triangulation` to perform this type of simulation. A good paper discussing these types of simulations is the paper 
# [_Comparing individual-based approaches to modelling the self-organization of multicellular tissues_](https://doi.org/10.1371/journal.pcbi.1005387) by Osborne et al. (2017).
#
# [^1]: If you are interested in hwot ehse ideas are applied in one dimension, see also [EpithelialDynamics1D.jl](https://github.com/DanielVandH/EpithelialDynamics1D.jl).
#
# ## Cell model 
# Let us consider a domain $\Omega$, and suppose we have an initial 
# set of points $\mathcal P = \{\vb r_1, \ldots, \vb r_N\}$ in the plane. We use a spring based model based on Hooke's law and Newton's laws,
# writing 
#
# ```math 
# \dv{\vb r_i}{t} = \alpha \sum_{j \in \mathcal N_i} \left(\|\vb r_{ji}\| - s\right)\hat{\vb r}_{ji}, \quad i=1,2,\ldots,N,
# ```
# where $\mathcal N_i$ is the set of points neighbouring $\vb r_i$ in the triangulation $\mathcal D\mathcal T(\mathcal P)$, 
# $\vb r_{ji} = \vb r_j - \vb r_i$, $\hat{\vb r}_{ji} = \vb r_{ji} / \|\vb r_{ji}\|$, $s$ is the resting spring length, 
# and $\alpha$ is the ratio of the spring constant $k$ to the damping constant $\eta$.
# We will allow the boundary nodes to move freely rather than consider some types of boundary conditions.
#
# Now consider proliferation. We interpret the cells as being the Voronoi polygons associated with each $\vb r_i$, i.e. the cell associated 
# with $\vb r_i$ is $\mathcal V_i$. We allow only one 
# cell to proliferate over a time interval $[t, t + \Delta t)$ and, given that a cell does proliferate, the probability that $\mathcal V_i$ proliferates 
# is $G_i\Delta t$ with $G_i = \max\{0, \beta(1 - 1/(KA_i))\}$, where $\beta$ is the intrinsic proliferation rate, $K$ is the cell carrying 
# capacity density, and $A_i$ is the area of $\mathcal V_i$. For cells on the boundary, $G_i = 0$ so that they do not proliferate; this is done 
# to avoid issues with unbounded cells. The probability any proliferation event occurs is given by $\sum_i G_i\Delta t$. When a cell $\mathcal V_i$
# is selected to proliferate, we select a random angle $\theta \in [0, 2\pi)$ and place a new point $\vb r_i'$ at $\vb r_i' = \vb r_i + \delta \epsilon \vb d$,
# where $\vb d = (\cos\theta,\sin\theta)$, $\delta = \operatorname{dist}(\vb r_i, \mathcal V_i)$, and $\epsilon$ is a small separation constant in $[0, 1)$.
# 
# For our simulation, we will first determine if any proliferation event occurs and, if so, update the triangulation with the new $\vb r_i'$. To then 
# integrate forward in time, we use Euler's method, writing 
# 
# ```math 
# \vb r_i(t + \Delta t) = \vb r_i(t) + \alpha \Delta t \sum_{j \in \mathcal N_i} \left(\|\vb r_{ji}\| - s\right)\hat{\vb r}_{ji}, \quad i=1,2,\ldots,N.
# ```
#
# We update the triangulation after each step.[^2]

# [^2]: There are efficient methods for updating the triangulation after small perturbations, e.g. the paper [_Star splaying: an algorithm for repairing Delaunay triangulations and convex hulls_](https://doi.org/10.1145/1064092.1064129) by Shewchuk (2005). For this implementation, we do not concern ourselves with being overly efficient, and simply retriangulate.

# ## Implementation 
# Let us now implement this model. First, we define a struct for storing the parameters of our model. 
using DelaunayTriangulation
using StableRNGs
using LinearAlgebra
using StatsBase
using CairoMakie
@kwdef mutable struct CellModel{P}
    tri::Triangulation{P}
    new_r_cache::Vector{NTuple{2,Float64}} # for r(t + Δt)
    const α::Float64
    const s::Float64
    const Δt::Float64
    const rng::StableRNGs.LehmerRNG
    const final_time::Float64
    const β::Float64
    const K::Float64
    const ϵ::Float64
end

# Let's now write functions for performing the migration step.
function migrate_cells!(cells::CellModel) # a more efficient way would be to loop over edges rather than vertices
    tri = cells.tri
    for i in each_solid_vertex(tri)
        rᵢ = get_point(tri, i)
        F = 0.0
        for j in get_neighbours(tri, i)
            DelaunayTriangulation.is_ghost_vertex(j) && continue
            rⱼ = get_point(tri, j)
            rᵢⱼ = rⱼ .- rᵢ
            F = F .+ cells.α * (norm(rᵢⱼ) - cells.s) .* rᵢⱼ ./ norm(rᵢⱼ)
        end
        cells.new_r_cache[i] = rᵢ .+ cells.Δt .* F
    end
    for i in each_solid_vertex(tri)
        DelaunayTriangulation.set_point!(tri, i, cells.new_r_cache[i])
    end
    cells.tri = retriangulate(tri; cells.rng)
    return nothing
end

# Now we can write the proliferation functions. First, let us write a function that computes the Voronoi areas. If we had the 
# `VoronoiTessellation` computed, we would just use `get_area`, but we are aiming to avoid having to compute $\mathcal V\mathcal T(\mathcal P)$
# directly. 
function polygon_area(points) # this is the same function from the Interpolation section
    n = DelaunayTriangulation.num_points(points)
    p, q, r, s = get_point(points, 1, 2, n, n - 1)
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    sx, sy = getxy(s)
    area = px * (qy - ry) + rx * (py - sy)
    for i in 2:(n-1)
        p, q, r = get_point(points, i, i + 1, i - 1)
        px, py = getxy(p)
        qx, qy = getxy(q)
        rx, ry = getxy(r)
        area += px * (qy - ry)
    end
    return area / 2
end
function get_voronoi_area(tri::Triangulation, i)
    points = NTuple{2,Float64}[]
    !DelaunayTriangulation.has_vertex(tri, i) && return (0.0, points) # might not be included anymore due to retriangulation
    DelaunayTriangulation.is_boundary_node(tri, i)[1] && return (0.0, points) # to prevent boundary cells from proliferating
    N = get_neighbours(tri, i)
    N₁ = first(N)
    j = N₁
    k = get_adjacent(tri, i, j)
    Ttype = DelaunayTriangulation.triangle_type(tri)
    for _ in 1:num_neighbours(tri, i)
        T = DelaunayTriangulation.construct_triangle(Ttype, i, j, k)
        c = DelaunayTriangulation.triangle_circumcenter(tri, T)
        push!(points, c)
        j = k
        k = get_adjacent(tri, i, j)
    end
    push!(points, points[begin])
    return polygon_area(points), points
end

# The function `get_voronoi_area` above returns `0` if the cell is on the boundary. Finally, our function for performing the 
# proliferation step is below.
function proliferate_cells!(cells::CellModel)
    E = 0.0
    Δt = cells.Δt
    tri = cells.tri
    areas = [get_voronoi_area(tri, i)[1] for i in DelaunayTriangulation.each_point_index(tri)] # don't use each_solid_vertex since each_solid_vertex is not ordered
    for i in DelaunayTriangulation.each_point_index(tri)
        Gᵢ = iszero(areas[i]) ? 0.0 : max(0.0, cells.β * (1 - 1 / (cells.K * areas[i])))
        E += Gᵢ * Δt
    end
    u = rand(cells.rng)
    if u < E
        areas ./= sum(areas)
        weights = Weights(areas)
        i = sample(cells.rng, weights) # select cell to proliferate randomly based on area
        θ = 2π * rand(cells.rng)
        p = get_point(tri, i)
        poly = get_voronoi_area(tri, i)[2]
        δ = DelaunayTriangulation.distance_to_polygon(p, get_points(tri), poly)
        s, c = sincos(θ)
        q = p .+ (δ * cells.ϵ) .* (c, s)
        add_point!(tri, q; rng=cells.rng)
        push!(cells.new_r_cache, q)
    end
    return nothing
end

# Finally, our simulation function is below.
function perform_step!(cells::CellModel)
    proliferate_cells!(cells)
    migrate_cells!(cells)
    return nothing
end
function simulate_cells(cells::CellModel)
    t = 0.0
    all_points = Vector{Vector{NTuple{2,Float64}}}()
    push!(all_points, deepcopy(get_points(cells.tri)))
    while t < cells.final_time
        perform_step!(cells)
        t += cells.Δt
        push!(all_points, deepcopy(get_points(cells.tri)))
    end
    return all_points
end

# ## Example
# Let us now give an example. Our initial set of points will be randomly chosen inside 
# the rectangle $[-2, 2] \times [-5, 5]$. We use $\alpha = 5$, $s = 2$, $\Delta t = 10^{-3}$,
# $\beta = 0.25$, $K = 100^2$, and $\epsilon = 0.5$. 
rng = StableRNG(123444)
a, b, c, d = -2.0, 2.0, -5.0, 5.0
points = [(a + (b - a) * rand(rng), c + (d - c) * rand(rng)) for _ in 1:10]
tri = triangulate(points; rng=rng)
cells = CellModel(; tri=tri, new_r_cache=similar(points), α=5.0, s=2.0, Δt=1e-3,
    β=0.25, K=100.0^2, rng, final_time=25.0, ϵ=0.5)
results = simulate_cells(cells);

fig = Figure(fontsize=26)
title_obs = Observable(L"t = %$(0.0)")
ax1 = Axis(fig[1, 1], width=1200, height=400, title=title_obs, titlealign=:left)
Δt = cells.Δt
i = Observable(1)
voronoiplot!(ax1, @lift(voronoi(triangulate(results[$i]; rng), clip=true, rng=rng)),
    color=:darkgreen, strokecolor=:black, strokewidth=2, show_generators=false)
xlims!(ax1, -12, 12)
ylims!(ax1, -12, 12)
resize_to_layout!(fig)
t = 0:Δt:cells.final_time
record(fig, "cell_simulation.mp4", 1:10:length(t); framerate=60) do ii
    i[] = ii
    title_obs[] = L"t = %$(((ii-1) * Δt))"
end;

# ![](cell_simulation.mp4)
