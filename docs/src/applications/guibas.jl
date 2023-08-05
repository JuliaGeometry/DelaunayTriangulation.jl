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
# Let p₋₃, p₋₂, p₋₁, be three points at infinity such that the triangle 
#   T₋₃₋₂₋₁ contains the input point set P = {p₁, ..., pₙ}.
# Initialise the triangulation with T₋₃₋₂₋₁
# for k in 1:n 
#     Sample pᵣ ∈ P randomly such that pᵣ is not yet in the triangulation.
#     Find the triangle Tᵢⱼₖ containing pᵣ.
#     Connect the vertices of Tᵢⱼₖ to pᵣ, thereby replacing Tᵢⱼₖ with three 
#       new triangles Tᵣₖᵢ, Tᵣⱼₖ, Tᵣᵢⱼ. 
#     For all the new edges, flip them if they are not legal. 
#       Apply this step recursively to all the new edges produced from the flipping.
# end 
# Remove the triangles containing p₋₃, p₋₂, p₋₁.
# Return the triangulation.
# ```
# 
# Of course, these lines of pseudocode hide a lot of details.