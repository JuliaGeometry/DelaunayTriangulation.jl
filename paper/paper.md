---
title: 'DelaunayTriangulation.jl: A Julia package for Delaunay triangulations and Voronoi tessellations in the plane'
tags:
    - Julia
    - geometry
    - visualisation
    - computational geometry
    - triangulation
    - tessellation
authors:
    - name: Daniel J. VandenHeuvel
    - orcid: https://orcid.org/0000-0001-6462-0135
    - affiliation: 1
    - corresponding: true
affiliations:
    - name: Department of Mathematics, Imperial College London, UK
      index: 1
date: 26 June 2024
bibliography: paper.bib
---

# Summary 

Given a set of points $\mathcal P$, edges $\mathcal E$, and piecewise linear boundaries $\mathcal B$ together defining some domain $\Omega$, a _Delaunay triangulation_ is a subdivision of this domain into triangles. The vertices of the triangles come from $\mathcal P$, and each of the edges in $\mathcal E$ and $\mathcal B$ is present as an edge of some triangle [@cheng2013delaunay; @aurenhammer2013voronoi]. A related geometric structure is the _Voronoi tessellation_ that partitions the plane into convex polygons for each $p \in \mathcal P$ such that, for a given polygon, each point in that polygon is closer to the associated polygon's point than to any other $q \in \mathcal P$ [@cheng2013delaunay; @aurenhammer2013voronoi]. Both of these objects have many applications, such as point location [@mucke1999fast], solving differential equations [@cheng2013delaunay], and geospatial interpolation [@bobach2009natural; @ledoux2005efficient].

DelaunayTriangulation.jl is a feature-rich Julia [@bezanson2017julia] package for computing Delaunay triangulations and Voronoi tessellations. The package, amongst many other features, supports constrained triangulations, mesh refinement, centroidal Voronoi tessellations, and dynamic updates. Thanks to the speed and genericity of Julia, the package is both performant and robust---making use of ExactPredicats.jl for computing predicates with exact arithmetic [@lairez2024exact; @devillers2002efficient; @melquiond2007formally; @meyer2008fpg]---while still allowing for generic representations of geometric primitives.

# Statement of Need 

Delaunay triangulations and Voronoi tessellations have applications in a myriad of fields. Delaunay triangulations have been used for point location [@mucke1999fast], solving differential equations [@cheng2013delaunay; @golias1997delaunay; @ju2006adaptive; @cendes1983magnetic], route planning [@chen2010enhanced;  @yan2008path; @sakthivel2022solving; @zhihai2021dynamic],etc. Voronoi tessellations are typically useful when there is some notion of _influence_ associated with a point, and have been applied to problems such as geospatial interpolation [@bobach2009natural; @ledoux2005efficient], image processing [@du2006centroidal; @wang2009edge; @du1999centroidal], visualisation [@pinho2006voromap; @balzer2005voronoi], clustering [@reddy2012initialisation; @du1999centroidal; @schrieber1991voronoi], cell biology [@hermann2008delaunay; @wang2024calibration; @osborne2017comparing].

Several software packages with support for computing Delaunay triangulations and Voronoi tessellations in two dimensions already exist, such as _Triangle_ [@shewchuk1996triangle], _MATLAB_ [@MATLAB], _SciPy_ [@SciPy], _CGAL_ [@CGAL], and _Gmsh_ [@GMSH]. DelaunayTriangulation.jl is the most feature-rich of these and benefits from the high-performance of Julia to efficiently support many operations. Julia's multiple dispatch [@bezanson2017julia] 
is leveraged to allow for complete customisation in how a user wishes to represent geometric primitives such as points and domain boundaries, a useful feature for allowing users to represent primitives in a way that suits their application without needing to sacrfice performance. The [documentation](https://juliageometry.github.io/DelaunayTriangulation.jl/stable/) lists many more features, including its ability a wide range of domains, even those that are disjoint and with holes. 

DelaunayTriangulation.jl has already seen use in several areas. DelaunayTriangulation.jl was used for mesh generation in [@vandenheuvel2023computational] and is used for the `tricontourf`, `triplot`, and `voronoiplot` routines inside Makie.jl [@danisch2021makie]. The packages [FiniteVolumeMethod.jl](https://github.com/SciML/FiniteVolumeMethod.jl) [@vandenheuvel2024finite] and [NaturalNeighbours.jl](https://github.com/DanielVandH/NaturalNeighbours.jl) [@vandenheuvel2024natural] are also built directly on top of DelaunayTriangulation.jl. The design of boundaries in DelaunayTriangulation.jl has been motivated especially for the efficient representation of boundary conditions along different parts of a boundary for solving differential equations, and this is heavily utilised by FiniteVolumeMethod.jl.  

# Examples

We give two examples of how the package can be used. Many more examples are given in the [documentation](https://juliageometry.github.io/DelaunayTriangulation.jl/stable/), including [tutorials](https://juliageometry.github.io/DelaunayTriangulation.jl/stable/tutorials/overview/) and [fully detailed applications](https://juliageometry.github.io/DelaunayTriangulation.jl/stable/applications/overview/) such as cell simulations. To fully demonstrate the utility of the package, our examples follow realistic applications. We omit code used for plotting with Makie.jl [@danisch2021makie] in the examples below for space reasons. The complete code is available [here](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/blob/paper/paper/paper.jl).

For our first example, we consider a domain motivated by mean exit time. In particular, consider the problem
```math 
\begin{align*}
\begin{array}{rcll}
D\nabla^2 T(x, y) & = & -1 & (x, y) \in \Omega, \\
T(x, y) & = & 0 & (x, y) \in \Gamma_a, \\
T(x, y) & = & 0 & (x, y) = (x_s, y_s), \\
\nabla T(x, y) \cdot \hat{\boldsymbol n}(x, y) & = & 0 & (x, y) \in \Gamma_r. 
\end{array}
\end{align*}
```
Here, $T(x, y)$ denotes the mean exit time of a particle exiting $\Omega$ with diffusivity $D$ starting at $(x, y)$ [@redner2001guide; @carr2022mean], $\hat{\boldsymbol n}(x, y)$ is the unit normal vector field on $\Gamma_r$, $(x_s, y_s) = (0, 0)$, and the domain $\Omega$ with boundary $\partial\Omega = \Gamma_a \cup \Gamma_r$ is shown in Figure ref{fig0}. This setup defines a mean exit time where the particle can only exit through $\Gamma_a$ or through the sink $(x_s, y_s)$, and it gets reflected off of $\Gamma_r$.

![The domain $\Omega$. The red part of the boundary defines the absorbing boundary $\Gamma_a$, and the blue part defines the reflecting boundary $\Gamma_r$.\label{fig0}](figure0.png)

The code to generate a mesh of the domain in Figure \ref{fig0} is given below. We use curves to define the boundary so that curve-bounded refinement can be applied [@gosselin2009delaunay]. The resulting mesh is shown in Figure \ref{fig1}, together with a solution of the mean exit time problem with $D = 6.25 \times 10^{-4}$; FiniteVolumeMethod.jl [@vandenheuvel2024finite] is used to solve this problem, and the code for this can be found [here](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/blob/paper/paper/paper.jl).

```julia
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
```

![(a) The generated mesh using DelaunayTriangulation.jl for the mean exit time domain. The different parts of the boundary are shown with different coloured dots. (b) The solution to the mean exit time problem using the mesh from (a) together with FiniteVolumeMethod.jl[@vandenheuvel2024finite].](figure1.png)

We now give an example using Voronoi tessellations. Our example is motivated from Lloyd's algorithm for $k$-means clustering [@kanugo2002efficient; @du1999centroidal]. We generate $k$ random points and compute their centroidal Voronoi tessellation. We then generate data and label them according to which Voronoi cell they belong to.[^1] The code is given below.

[^1]: This example is somewhat contrived. Our procedure only corresponds to $k$-means clustering when the data set is dense [@kanugo2002efficient; @du1999centroidal]. To exactly obtain $k$-means clustering, the density function used for computing the centroids would need to be adjusted. See Section 2.3 of [@du1999centroidal] for more details.

```julia
using Random
Random.seed!(123)
k = 7
clusters = [(rand(), rand()) for _ in 1:k]
# Assume data live in [0, 1]²
push!(clusters, (0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0))
# Tessellate and smooth 
tri = triangulate(clusters)
vor = voronoi(tri, clip=true) # clips to [0, 1]²
cvor = centroidal_smooth(vor)
# Generate data and assign 
data = rand(2, 2500)
label = p -> get_nearest_neighbour(cvor, p)
labels = label.(eachcol(data))
```

![Example of $k$-means clustering. The polygons are the clusters, and each point is coloured according to which cluster it belongs to, computed using `get_nearest_neighbour`.](figure2.png)

# Extensions

There are still several features that are intended to eventually be implemented, some of these being:
1. Weighted triangulations and Voronoi treemaps, using the algorithms described in [@cheng2013delaunay; arlind2012computing].
2. Support for maximum angle constraints and generalised Steiner points, using algorithms and ideas described in [@ungor2009off; hale2009quality; hale2009computing].
3. Clipped Voronoi tessellations to arbitrary boundaries, possibly using the VoroCrust algorithm [@ahmed2020vorocrust].
4. Centroidal tessellations with inhomogeneous mass densities, as described in [@du1999centroidal].
5. Inserting curves into an existing triangulation [@gosselin2009delaunay; zaide2014inserting].
6. Delaunay repair algorithms for retriangulating perturbed point sets, using ideas from [@shewchuk2005star; yuanfeng2010fast]. 

There is no intention to support three-dimensional geometries within DelaunayTriangulation.jl. For this, the best option is TetGen [@hang2015tetgen].

# References