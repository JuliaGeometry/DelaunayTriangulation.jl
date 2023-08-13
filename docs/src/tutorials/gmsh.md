```@meta 
EditURL = "gmsh.md"
```

# Gmsh Integration 

This package provides some support for using Gmsh to compute constrained Delaunay triangulations. This is not so relevant with the latest versions of DelaunayTriangulation.jl now that we have native support for constrained triangulations, but this might be useful if, for example, you want to try [different mesh algorithms from Gmsh](https://gmsh.info/doc/texinfo/gmsh.html#Choosing-the-right-unstructured-algorithm). The support is also limited, in that we only support line segments rather than, for example, splines for curved boundaries, and interior holes and disjoint domains are not supported. Note also that we may remove Gmsh support once curved boundaries are supported natively.

To actually execute the Gmsh code, you need to have installed Gmsh. You can download it from [here](https://gmsh.info/#Download). You then define the path. In this tutorial, we use 

```julia 
GMSH_PATH = "./gmsh-4.11.1-Windows64/gmsh.exe"
```

We could have designed support based on, for example, [Gmsh.jl](https://github.com/JuliaFEM/Gmsh.jl), but we do not.

## Contiguous boundary 

Here we mesh a domain with a single non-segmented boundary curve. We show a comparison of the results with different refinement parameters. The default mesh algorithm used in `generate_mesh` is a 

```julia 
using DelaunayTriangulation
using CairoMakie
a = 4 / 5
t = LinRange(0, 2π, 100)
x = @. a * (2cos(t) + cos(2t))
y = @. a * (2sin(t) - sin(2t))
tri = generate_mesh(x, y, 0.1; gmsh_path=GMSH_PATH)
tri2 = generate_mesh(x, y, 1.0; gmsh_path=GMSH_PATH)
fig = Figure()
ax = Axis(fig[1, 1], xlabel="x", ylabel="y", width=300, height=300,
    title="Dense mesh", titlealign=:left)
triplot!(ax, tri, show_convex_hull=true, show_constrained_edges=true)
ax = Axis(fig[1, 2], xlabel="x", ylabel="y", width=300, height=300,
    title="Coarse mesh", titlealign=:left)
triplot!(ax, tri2, show_convex_hull=true, show_constrained_edges=true)
resize_to_layout!(fig)
fig
```

```@raw html
<figure>
    <img src='/gmsh_figures/gmsh_example_1.png', alt='Triangulation'><br>
</figure>
```