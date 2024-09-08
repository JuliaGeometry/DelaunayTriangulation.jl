# # Power Diagrams 
# 

# In this tutorial, we demonstrate how we can construct 
# power diagrams (also called weighted Voronoi tessellations).
# These are dual to the weighted Delaunay triangulation and, 
# instead of being based on the Euclidean metric like Voronoi tessellations,
# the power distance is used to define the Voronoi tiles. The power distance 
# is defined by $\pi(p, q) = d(p, q)^2 - w_p - w_q$, where $d(p, q)$
# is the Euclidean distance between points $p$ and $q$, and $w_p$ and $w_q$
# are _weights_ associated with $p$ and $q$. See [this page](../math/power.md)
# for more details. To start with the tutorial, we load in the packages we'll need. 
using DelaunayTriangulation
using CairoMakie
using StableRNGs
using StatsBase #src
using ReferenceTests #src
using Test #src
fig_path = joinpath(@__DIR__, "../figures") #src