# # Example Applications 
# ## Guibas et al.'s Algorithm 

# In this application, we implement an alternative method 
# for constructing unconstrained triangulations to what 
# we use in the package (namely, the Bowyer-Watson algorithm).
# This method, which we call Guibas et al.'s algorithm (see the paper
# *Randomized Incremental Construction of Delaunay and Voronoi Diagrams* 
# by Guibas et al., 1992), is described by de Berg et al. in the 
# book *Computational Geometry: Algorithms and Applications* (3 ed., 2008). 
# (We modify the algorithm in the book slightly to use three points at infinity rather 
# than two.)
#
# Guibas et al.'s algorithm is based on edge flipping, as described in the previous application, 
# but uses efficient point location and incremental construction 
# to make the procedure significantly faster and competitive with the 
# Bowyer-Watson method. (In fact, the package originally used this algorithm, 
# but it is not suitable for constrained triangulations with holes and cannot easily handle 
# vertex deletion.)
#
# At a high level, the algorithm is as follows. (This algorithm ignores degeneracies and collinear points.)
# ```julia
# Let p₋₁, p₋₂, p₋₃, be three points at infinity such that the triangle 
#   T₋₁₋₂₋₃ contains the input point set P = {p₁, ..., pₙ}.
# Initialise the triangulation with T₋₁₋₂₋₃
# for k in 1:n 
#     Sample pᵣ ∈ P randomly such that pᵣ is not yet in the triangulation.
#     Find the triangle Tᵢⱼₖ containing pᵣ.
#     Connect the vertices of Tᵢⱼₖ to pᵣ, thereby replacing Tᵢⱼₖ with three 
#       new triangles Tᵣₖᵢ, Tᵣⱼₖ, Tᵣᵢⱼ. 
#     For all the new edges, flip them if they are not legal. 
#       Apply this step recursively to all the new edges produced from the flipping.
# end 
# Remove the triangles containing p₋₁, p₋₂, p₋₃.
# Return the triangulation.
# ```
# 
# Of course, these lines of pseudocode hide a lot of details.
#
# The main complication is point location. Since each triangulation is on 
# top of the other, the idea is to somehow create a representation of the triangulation's 
# _history_, which can be used to find the triangle containing a point starting from 
# the initial triangle. In particular, we maintain a directed acyclic graph (DAG), 
# denoted $\mathcal D$, whose leaves give the current triangulation and the internal 
# nodes give past triangulations. Individual leaves correspond to individual triangles, 
# and proceeding up the DAG takes us to the original triangle that came before it, noting 
# that at each stage of the above algorithm a new point is used to split a triangle into 
# three. For example, suppose we split a triangle $T_{ijk}$ into three triangles. We 
# add three new leaves into $\mathcal D$, one for each new triangle, stemming from 
# $T_{ijk}$. Similar rules apply for edge flipping.
#
# To use $\mathcal D$ for point location, we proceed as follows. Starting at the root $T_{-1,-2,-3}$,
# we check the three children of the root to see if the point $p$ is inside it, and then 
# we enter the triangle that contains the point. This gives us a new node to continue 
# traversing in this manner. Eventually, we reach a leaf of $\mathcal D$, and this leaf is the 
# triangle in the current triangulation that contains $p$.
#
# The other issue is the choice of $p_{-1}$, $p_{-2}$, $p_{-3}$. While it is possible to 
# devise symbolic rules for treating these points without assigning actual values for 
# them, we will compute some actual values for them, remembering that (1) they must 
# define a counter-clockwise triangle, and (2) they must be large enough such that the triangle 
# $T_{-1, -2, -3}$ contains all the points.
#
### Implementation
#
# Let us now implement this algorithm. First, let us define the data structure we will be using 
# for the DAG. We call this the history graph, and use [SimpleGraphs.jl](https://github.com/scheinerman/SimpleGraphs.jl).
using SimpleGraphs
struct HistoryGraph 
    graph::DirectedGraph{NTuple{3,Int}}
    function HistoryGraph{NTuple{3,Int}}() 
        G = DirectedGraph{T}()
        forbid_loops!(G)
        return new{T}(G)
    end
    HistoryGraph(HG::DirectedGraph{T}) where {T} = new{T}(HG)
end
# Now we define a triangulation structure that wraps this graph together with a `Triangulation`.
struct GuibasTriangulation{T<:Triangulation}
    tri::T 
    history::HistoryGraph
end
# Next, let us define the point location function.
function dag_search(tri::GuibasTriangulation, r, init = (-1, -2, -3))
    if out_deg(tri, init) == 0 
        
end