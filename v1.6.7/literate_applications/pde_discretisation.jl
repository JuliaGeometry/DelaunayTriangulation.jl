# # Solving PDEs
# For our next application, we consider solving PDEs. In particular, 
# we focus on the discretisation of a domain on which a PDE can be solved.
# The method we discuss here is the basis of the package [FiniteVolumeMethod.jl](https://github.com/SciML/FiniteVolumeMethod.jl).
#
# ## Problem Statement 
# For this discussion, we will focus solely on solving a mean exit time problem of the form 
# ```math
# \begin{aligned}
# D\grad^2 T(\vb x) &= -1, \quad \vb x \in \Omega, \\
# T(\vb x) &= 0, \quad \vb x \in \partial \Omega,
# \end{aligned}
# ```
# where $T$ is the _mean exit time_, meaning the average amount of time it takes a particle to leave $\Omega$ through $\partial\Omega$ starting at 
# $\vb x$; $D$ is the diffusivity; and $\Omega$ is the domain of interest.[^1] For this discussion, we will let $\Omega$ be an annulus so that 
# $\Omega = \{(r, \theta) : 1 < r < 2,\, 0 < \theta < 2\pi\}$. We will exploit the linearity in this problem to express its solution as 
# being approximated by some linear system of the form $\vb A\vb T = \vb b$, but we note that for most PDEs nonlinear solvers are required.
#
# [^1]: See [this example](https://sciml.github.io/FiniteVolumeMethod.jl/stable/tutorials/mean_exit_time/) for more detail.
#
# ## Discretisation of the Domain
# Now let's discuss the discretisation of $\Omega$. The aim is to decompose $\Omega$ into a collection of triangles $\mathcal T_i$ such that 
# $\Omega = \cup_i\mathcal T_i$. We will then approximate the solution $T$ on each triangle $\mathcal T_i$ by a piecewise linear function.
# For this decomposition, we simply use $\mathcal D\mathcal T(\Omega)$, where $\mathcal D\mathcal T(\Omega)$ is some Delaunay triangulation of points in 
# $\Omega$ with $\partial\mathcal D\mathcal T(\Omega) = \partial\Omega$. Here is our discretisation of our annulus.
using DelaunayTriangulation
using CairoMakie
using StableRNGs
rng = StableRNG(123)
R₁ = 1.0
R₂ = 2.0
outer_circle = CircularArc((R₂, 0.0), (R₂, 0.0), (0.0, 0.0))
inner_circle = CircularArc((R₁, 0.0), (R₁, 0.0), (0.0, 0.0), positive = false)
points = NTuple{2, Float64}[]
tri = triangulate(points; rng, boundary_nodes = [[[outer_circle]], [[inner_circle]]])
A = 2π * (R₂^2 - R₁^2)
refine!(tri; max_area = 2.0e-3A, min_angle = 33.0, rng)
fig, ax, sc = triplot(tri)
fig

# ## Discretisation of the PDE 
# Now let's determine how we can use the discretisation above to solve the PDE. Central to this approach is the idea of a _control volume_ around 
# each point. In particular, connect the centroids of each triangle to the midpoints of the edges of the triangle. This defines a collection of polygons 
# $\Omega_i$ around each point $\vb x_i$, as shown below in blue.
points = NTuple{2, Float64}[] #hide
for T in each_solid_triangle(tri) #hide
    u, v, w = triangle_vertices(T) #hide
    p, q, r = get_point(tri, u, v, w) #hide
    c = DelaunayTriangulation.triangle_centroid(p, q, r) #hide
    push!(points, c, (p .+ q) ./ 2, c, (q .+ r) ./ 2, c, (r .+ p) ./ 2) #hide
end #hide
linesegments!(ax, points, color = :blue) #hide
fig #hide

# Consider a particular control volume $\Omega_i$. We integrate our PDE over this domain, writing 
# $\grad^2 = \grad \vdot \grad$ and use the divergence theorem:
# ```math 
# \begin{aligned}
# D\iint_{\Omega_i} \grad^2 T(\vb x)\, \dd A &= -\iint_{\Omega_i}\, \dd A \\
# D\oint_{\partial\Omega_i} \grad T(\vb x) \vdot \vu{n}\, \dd s &= -A_i,
# \end{aligned}
# ```
# where $\vu{n}$ is the outward normal to $\Omega_i$ and $A_i$ is the area of $\Omega_i$. Now, 
# let $\mathcal E_i$ be the set of edges defining $\partial\Omega_i$ so that $\partial\Omega_i = \cup_{\sigma\in\mathcal E_i} \sigma$. Thus, 
# ```math 
# D\sum_{\sigma\in\mathcal E_i}\int_\sigma \grad T(\vb x) \vdot \vu{n}\, \dd s = -A_i. 
# ```
# Now we want to approximate each $\int_\sigma \grad T(\vb x) \vdot \vu{n} \, \dd s$. First, we use a midpoint 
# approximation to give 
# ```math
# \int_{\sigma} \grad T(\vb x) \vdot \vu{n}\, \dd s \approx \left[\grad T(\vb x_\sigma) \vdot \vu{n}_\sigma\right] L_\sigma, 
# ```
# where $\vb x_\sigma$ is the midpoint of $\sigma$, and $\vu{n}_\sigma$ is the normal to $\sigma$. Now, 
# to approximate $\grad T(\vb x_\sigma)$, we use a piecewise linear approximation to $T$ on each triangle. In particular, 
# suppose $\vb x_\sigma$ is inside some triangle $\mathcal T_k$. We use a linear shape function inside $\mathcal T_k$ for 
# approximating $T$, writing 
# ```math 
# T(\vb x) = \alpha_k x + \beta_k y + \gamma_k, \quad \vb x \in \mathcal T_k. 
# ```
# Using Cramer's rule, we can determine $\alpha_k$, $\beta_k$, and $\gamma_k$ in terms of the values of $T$ at the vertices of $\mathcal T_k$. In particular, 
# ```math 
# \begin{aligned}
# \alpha_k &= s_{k, 11} T_{k1} + s_{k, 12} T_{k2} + s_{k, 13} T_{k3}, \\
# \beta_k &= s_{k, 21} T_{k1} + s_{k, 22} T_{k2} + s_{k, 23} T_{k3}, \\
# \gamma_k &= s_{k, 31} T_{k1} + s_{k, 32} T_{k2} + s_{k, 33} T_{k3},
# \end{aligned}
# ```
# where 
# ```math
# \begin{aligned}
# \vb S_k &= \frac{1}{\Delta_k} \begin{bmatrix} y_{k2} - y_{k3} & y_{k3} - y_{k1} & y_{k1} - y_{k2} \\ x_{k3} - x_{k2} & x_{k1} - x_{k3} & x_{k2} - x_{k1} \\ x_{k2}y_{k3} - x_{k3}y_{k2} & x_{k3}y_{k1} - x_{k1}y_{k3} & x_{k1}y_{k2} - x_{k2}y_{k1} \end{bmatrix}, \\
# \Delta_k &= x_{k1}y_{k2}-x_{k2}y_{k1}-x_{k1}y_{k3}+x_{k3}y_{k1}+x_{k2}y_{k3}-x_{k3}y_{k2}.
# \end{aligned}
# ```
# With this approximation, $\grad T(\vb x_\sigma) \approx (\alpha_k, \beta_k)^{\mathsf T\mkern-1.5mu}$. Thus, 
# have 
# ```math 
# \frac{D}{A_i}\sum_{\sigma \in \mathcal E_i} \left[\left(s_{k, 11}n_\sigma^x + s_{k, 21}n_\sigma^y\right)T_{k1} + \left(s_{k, 12}n_\sigma^x + s_{k, 22}n_\sigma^y\right)T_{k2} + \left(s_{k, 13}n_\sigma^x + s_{k,23}n_\sigma^y\right)T_{k3}\right]L_\sigma = -1.
# ```
#
# ## Building a Matrix System 
# Notice that the equation above is linear in the unknowns $T_{ki}$. Thus, we can write this as a matrix system $\vb A\vb T = \vb b$. Consider 
# the $i$th row of $\vb A$. This row corresponds to the control volume $\Omega_i$ and is given by $\vb a_i^{\mathsf T\mkern-1.5mu}$, where 
# the non-zero terms of $\vb a_i$ correspond to the vertices of the triangles that make up $\Omega_i$. The corresponding value of $b_i$ in $\vb b$ 
# is simply $-1$. For handling the boundaries, we will instead let $\vb a_i = \vb e_i$ for any nodes $\vb x_i$ on the boundary, where $\vb e_i = (0, \ldots, 1, \ldots, 0)$
# is the $i$th standard basis vector, and $b_i = 0$. The solution to this system will give us the mean exit time $T$.
# 
# ## Implementation
# Let's now implement these ideas. Some details of this implementation, like how we efficient loop over the mesh for building the system by using 
# edges rather than vertices, have been skipped and are described in more detail [here](https://sciml.github.io/FiniteVolumeMethod.jl/dev/math/).
using LinearAlgebra
using SparseArrays
function solve_met_problem(tri::Triangulation, D)
    ## To start, we need to build a map that takes the vertices from tri 
    ## into a range of consecutive integers, since not all vertices are used. 
    vertex_map = Dict{Int, Int}()
    inverse_vertex_map = Dict{Int, Int}()
    cur_idx = 1
    for i in DelaunayTriangulation.each_point_index(tri)
        if DelaunayTriangulation.has_vertex(tri, i)
            vertex_map[i] = cur_idx
            inverse_vertex_map[cur_idx] = i
            cur_idx += 1
        end
    end
    ## Next, we need to build up what we need from the geometry.
    nv = num_solid_vertices(tri)
    nt = num_solid_triangles(tri)
    cv_volumes = zeros(nv)
    Ttype = DelaunayTriangulation.triangle_type(tri)
    shape_function_coefficients = Dict{Ttype, NTuple{9, Float64}}()
    cv_edge_midpoints = Dict{Ttype, NTuple{3, NTuple{2, Float64}}}()
    cv_edge_normals = Dict{Ttype, NTuple{3, NTuple{2, Float64}}}()
    cv_edge_lengths = Dict{Ttype, NTuple{3, Float64}}()
    sizehint!.((cv_volumes, shape_function_coefficients, cv_edge_midpoints, cv_edge_normals, cv_edge_lengths), nt)
    for T in each_solid_triangle(tri)
        u, v, w = triangle_vertices(T)
        p, q, r = get_point(tri, u, v, w)
        px, py = getxy(p)
        qx, qy = getxy(q)
        rx, ry = getxy(r)
        cx, cy = (px + qx + rx) / 3, (py + qy + ry) / 3
        m₁x, m₁y = (px + qx) / 2, (py + qy) / 2
        m₂x, m₂y = (qx + rx) / 2, (qy + ry) / 2
        m₃x, m₃y = (rx + px) / 2, (ry + py) / 2
        ## Connect the centroid to each vertex, and all the midpoints to each other 
        pcx, pcy = cx - px, cy - py
        qcx, qcy = cx - qx, cy - qy
        rcx, rcy = cx - rx, cy - ry
        m₁₃x, m₁₃y = m₁x - m₃x, m₁y - m₃y
        m₂₁x, m₂₁y = m₂x - m₁x, m₂y - m₁y
        m₃₂x, m₃₂y = m₃x - m₂x, m₃y - m₂y
        ## Get the sub-control volume areas 
        S₁ = 1 / 2 * abs(pcx * m₁₃y - pcy * m₁₃x)
        S₂ = 1 / 2 * abs(qcx * m₂₁y - qcy * m₂₁x)
        S₃ = 1 / 2 * abs(rcx * m₃₂y - rcy * m₃₂x)
        cv_volumes[vertex_map[u]] += S₁
        cv_volumes[vertex_map[v]] += S₂
        cv_volumes[vertex_map[w]] += S₃
        ## Now get the shape function coefficients 
        Δ = qx * ry - qy * rx - px * ry + rx * py + px * qy - qx * py
        s₁ = (qy - ry) / Δ
        s₂ = (ry - py) / Δ
        s₃ = (py - qy) / Δ
        s₄ = (rx - qx) / Δ
        s₅ = (px - rx) / Δ
        s₆ = (qx - px) / Δ
        s₇ = (qx * ry - rx * qy) / Δ
        s₈ = (rx * py - px * ry) / Δ
        s₉ = (px * qy - qx * py) / Δ
        shape_function_coefficients[T] = (s₁, s₂, s₃, s₄, s₅, s₆, s₇, s₈, s₉)
        ## Get the midpoints 
        m₁c = (m₁x + cx) / 2, (m₁y + cy) / 2
        m₂c = (m₂x + cx) / 2, (m₂y + cy) / 2
        m₃c = (m₃x + cx) / 2, (m₃y + cy) / 2
        cv_edge_midpoints[T] = (m₁c, m₂c, m₃c)
        ## Now get the normal vectors on the control volume edges
        e₁x, e₁y = cx - m₁x, cy - m₁y
        e₂x, e₂y = cx - m₂x, cy - m₂y
        e₃x, e₃y = cx - m₃x, cy - m₃y
        ℓ₁ = norm((e₁x, e₁y))
        ℓ₂ = norm((e₂x, e₂y))
        ℓ₃ = norm((e₃x, e₃y))
        cv_edge_lengths[T] = (ℓ₁, ℓ₂, ℓ₃)
        n₁ = e₁y / ℓ₁, -e₁x / ℓ₁
        n₂ = e₂y / ℓ₂, -e₂x / ℓ₂
        n₃ = e₃y / ℓ₃, -e₃x / ℓ₃
        cv_edge_normals[T] = (n₁, n₂, n₃)
    end
    ## Now we can build A 
    A = zeros(nv, nv)
    for T in each_solid_triangle(tri)
        u, v, w = triangle_vertices(T)
        p, q, r = get_point(tri, u, v, w)
        s₁₁, s₁₂, s₁₃, s₂₁, s₂₂, s₂₃, s₃₁, s₃₂, s₃₃ = shape_function_coefficients[T]
        for (edge_index, e₁₂) in enumerate(DelaunayTriangulation.triangle_edges(T))
            e₁, e₂ = edge_vertices(e₁₂)
            x, y = cv_edge_midpoints[T][edge_index]
            nx, ny = cv_edge_normals[T][edge_index]
            ℓ = cv_edge_lengths[T][edge_index]
            Dℓ = D * ℓ
            a123 = (
                Dℓ * (s₁₁ * nx + s₂₁ * ny),
                Dℓ * (s₁₂ * nx + s₂₂ * ny),
                Dℓ * (s₁₃ * nx + s₂₃ * ny),
            )
            e1_is_bnd = DelaunayTriangulation.is_boundary_node(tri, e₁)[1]
            e2_is_bnd = DelaunayTriangulation.is_boundary_node(tri, e₂)[1]
            for vert in 1:3
                e1_is_bnd || (A[vertex_map[e₁], vertex_map[(u, v, w)[vert]]] += a123[vert] / cv_volumes[vertex_map[e₁]])
                e2_is_bnd || (A[vertex_map[e₂], vertex_map[(u, v, w)[vert]]] -= a123[vert] / cv_volumes[vertex_map[e₂]])
            end
        end
    end
    Asp = sparse(A)
    ## Now we can build b
    b = zeros(nv)
    for i in each_solid_vertex(tri)
        if !DelaunayTriangulation.is_boundary_node(tri, i)[1]
            b[vertex_map[i]] = -1.0
        else
            A[vertex_map[i], vertex_map[i]] = 1.0 # b[i] is already 0 
        end
    end
    ## Now solve and return the solution
    T = A \ b
    filled_out_T = zeros(DelaunayTriangulation.num_points(tri)) # make sure that T[i] actually refers to the vertex i
    for i in 1:nv
        filled_out_T[inverse_vertex_map[i]] = T[i]
    end
    return filled_out_T
end

# ## Solving the System
# Let's now solve this problem, taking $D = 6.25 \times 10^{-4}$.
D = 6.25e-4
T = solve_met_problem(tri, D)
fig, ax, sc = tricontourf(tri, T, levels = 0:5:200, extendhigh = :auto)
fig
