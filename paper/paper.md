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
      orcid: 0000-0001-6462-0135
      affiliation: 1
      corresponding: true
affiliations:
    - name: Department of Mathematics, Imperial College London, UK
      index: 1
date: 10 September 2024
bibliography: paper.bib
---

# Summary 

DelaunayTriangulation.jl is a feature-rich Julia [@bezanson2017julia] package for computing Delaunay triangulations and Voronoi tessellations. The package, amongst many other features, supports unconstrained and constrained triangulations, mesh refinement, clipped and centroidal Voronoi tessellations, power diagrams, and dynamic updates. Thanks to the speed and genericity of Julia, the package is both performant and robust---making use of [AdaptivePredicates.jl](https://github.com/JuliaGeometry/AdaptivePredicates.jl) [@churavy2024adaptive; @shewchuk1997adaptive] and [ExactPredicates.jl](https://github.com/lairez/ExactPredicates.jl) [@lairez2024exact] for computing predicates with robust arithmetic---while still allowing for generic representations of geometric primitives.

Given a set of points $\mathcal P$, edges $\mathcal E$, and piecewise linear boundaries $\mathcal B$ together defining some domain $\Omega$, a _Delaunay triangulation_ is a subdivision of this domain into triangles. The vertices of the triangles come from $\mathcal P$, and each of the edges in $\mathcal E$ and $\mathcal B$ are present as an edge of some triangle [@cheng2013delaunay; @aurenhammer2013voronoi]. The boundaries $\mathcal B$ may also be given as parametric curves, in which case they are discretised until they accurately approximate the curved boundary [@gosselin2009delaunay]. A related geometric structure is the _Voronoi tessellation_ that partitions the plane into convex polygons for each $p \in \mathcal P$ such that, for a given polygon, each point in that polygon is closer to the associated polygon's point than to any other $q \in \mathcal P$ [@cheng2013delaunay; @aurenhammer2013voronoi]. Weighted triangulations and power diagrams are generalisations of these structures that allow for the inclusion of weights associated with the points [@cheng2013delaunay].

# Statement of Need 

Delaunay triangulations and Voronoi tessellations have applications in a myriad of fields. Delaunay triangulations have been used for point location [@mucke1999fast], solving differential equations [@golias1997delaunay; @ju2006adaptive], path planning [@yan2008path], etc. Voronoi tessellations are typically useful when there is some notion of _influence_ associated with a point, and have been applied to problems such as geospatial interpolation [@bobach2009natural], image processing [@du1999centroidal], and cell biology [@hermann2008delaunay; @wang2024calibration].

Several software packages with support for computing Delaunay triangulations and Voronoi tessellations in two dimensions already exist, such as [_Triangle_](https://www.cs.cmu.edu/~quake/triangle.html) [@shewchuk1996triangle], [_MATLAB_](https://uk.mathworks.com/help/matlab/computational-geometry.html?s_tid=CRUX_lftnav) [@MATLAB], [_SciPy_](https://docs.scipy.org/doc/scipy/tutorial/spatial.html) [@SciPy], [_CGAL_](https://www.cgal.org/) [@CGAL], and [_Gmsh_](https://gmsh.info/) [@GMSH]. There are also other Julia packages supporting some of these features, although none are as developed as DelaunayTriangulation.jl; a comparison with these other packages is given in DelaunayTriangulation.jl's [README](https://github.com/JuliaGeometry/DelaunayTriangulation.jl?tab=readme-ov-file#similar-packages). DelaunayTriangulation.jl supports a larger set of features than most of these other software packages, such as triangulation of curve-bounded domains and power diagrams, and benefits from the high-performance of Julia to efficiently support many operations. Julia's multiple dispatch [@bezanson2017julia] 
is leveraged to allow for complete customisation in how a user wishes to represent geometric primitives such as points and domain boundaries, a useful feature for allowing users to represent primitives in a way that suits their application without needing to sacrfice performance. The [documentation](https://juliageometry.github.io/DelaunayTriangulation.jl/stable/) lists many more features, including its ability to represent a wide range of domains, even those that are disjoint and with holes.

DelaunayTriangulation.jl has already seen use in several areas. DelaunayTriangulation.jl was used for mesh generation in @vandenheuvel2023computational and is used for the `tricontourf`, `triplot`, and `voronoiplot` routines inside [Makie.jl](https://github.com/MakieOrg/Makie.jl) [@danisch2021makie]. The packages [FiniteVolumeMethod.jl](https://github.com/SciML/FiniteVolumeMethod.jl) [@vandenheuvel2024finite] and [NaturalNeighbours.jl](https://github.com/DanielVandH/NaturalNeighbours.jl) [@vandenheuvel2024natural] are also built directly on top of DelaunayTriangulation.jl. 

# Example

We give one example of how the package can be used, focusing on Delaunay triangulations rather than Voronoi tessellations. Many more examples are given in the [documentation](https://juliageometry.github.io/DelaunayTriangulation.jl/stable/), including [tutorials](https://juliageometry.github.io/DelaunayTriangulation.jl/stable/tutorials/overview/) and [fully detailed applications](https://juliageometry.github.io/DelaunayTriangulation.jl/stable/applications/overview/) such as cell simulations. To fully demonstrate the utility of the package, we consider a realistic application. We omit code used for plotting with Makie.jl [@danisch2021makie] in the example below for space reasons. The complete code is available [here](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/blob/paper/paper/paper.jl).

We consider a domain motivated by mean exit time. In particular, consider the problem
$$
\begin{array}{rcll}
D\nabla^2 T(x, y) & = & -1 & (x, y) \in \Omega, \\
T(x, y) & = & 0 & (x, y) \in \Gamma_a, \\
T(x, y) & = & 0 & (x, y) = (x_s, y_s), \\
\nabla T(x, y) \cdot \hat{\boldsymbol n}(x, y) & = & 0 & (x, y) \in \Gamma_r. 
\end{array}
$$
Here, $T(x, y)$ denotes the mean exit time of a particle exiting $\Omega$ with diffusivity $D$ starting at $(x, y)$ [@redner2001guide; @carr2022mean], $\hat{\boldsymbol n}(x, y)$ is the unit normal vector field on $\Gamma_r$, $(x_s, y_s) = (0, 0)$, and the domain $\Omega$ with boundary $\partial\Omega = \Gamma_a \cup \Gamma_r$ is shown in \autoref{fig:1}(a). This setup defines a mean exit time where the particle can only exit through $\Gamma_a$ or through the sink $(x_s, y_s)$, and it gets reflected off of $\Gamma_r$.

The code to generate a mesh of the domain is given below. We use `CircularArc`s to define the boundary so that curve-bounded refinement can be applied using the algorithm of @gosselin2009delaunay. The resulting mesh is shown in \autoref{fig:1}, together with a solution of the mean exit time problem with $D = 6.25 \times 10^{-4}$; FiniteVolumeMethod.jl [@vandenheuvel2024finite] is used to solve this problem, and the code for this can be found [here](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/blob/paper/paper/paper.jl).

```julia
# The outer circle
θ = 5π / 64
cs = θ -> (cos(θ), sin(θ))
p1, q1 = cs(-π / 2 - θ), cs(θ) # Absorbing 
p2, q2 = q1, cs(π / 2 - θ)     # Reflecting 
p3, q3 = q2, cs(π + θ)         # Absorbing 
p4, q4 = q3, p1                # Reflecting
c0 = (0.0, 0.0)
C01 = CircularArc(p1, q1, c0) # first, last, center
C02 = CircularArc(p2, q2, c0)
C03 = CircularArc(p3, q3, c0)
C04 = CircularArc(p4, q4, c0)
C0 = [[C01], [C02], [C03], [C04]]
# Inner circles
c1, p5 = (-0.4, -0.4), (-0.65, -0.65)
c2, p6 = (0.4, 0.4), (0.65, 0.65)
C1 = CircularArc(p5, p5, c1, positive=false) # Reflecting
C2 = CircularArc(p6, p6, c2, positive=false) # Reflecting
# Triangulate and refine
sink = (0.0, 0.0)
tri = triangulate([sink], boundary_nodes=[C0, [[C1]], [[C2]]])
refine!(tri; max_area=1e-3get_area(tri))
```

![(a) The generated mesh using DelaunayTriangulation.jl for the mean exit time domain. The red dots along the boundary define the absorbing part of the boundary, $\Gamma_a$, and the blue dots define the reflecting part, $\Gamma_r$. (b) The solution to the mean exit time problem using the mesh from (a) together with FiniteVolumeMethod.jl [@vandenheuvel2024finite].\label{fig:1}](figure1.png)

# References