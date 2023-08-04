# # Triangulation Operations 
# ## Adding or Clearing Ghost Triangles 

# An important component to be aware of when considering 
# dynamic updates to triangulations are the ghost triangles. 
# As we discussed in the [vertex insertion/deletion example](operations_vertex_insertion_deletion.md),
# ghost triangles are needed when we are making updates outside of the boundary 
# of the current triangulation. 
using DelaunayTriangulation
using CairoMakie
using ReferenceTests #src
using Test #src
fig_path = joinpath(@__DIR__, "../figures") #src

# Let us take an example triangulation. 
points = [(-1.0, -1.0), (1.0, -1.0), (0.0, 1.0)]
tri = triangulate(points)

#-
fig, ax, sc = triplot(tri, show_ghost_edges = true)
fig

# The ghost triangles are represented by the convex regions 
# bounded by the blue lines. Now, as shown by the output above, 
# we can see that there are no ghost triangles. To add them, 
# we can do 
add_ghost_triangles!(tri)
tri
@test DelaunayTriangulation.has_ghost_triangles(tri) #src

# or we can just use `triangulate` with `delete_ghosts = false`:
tri = triangulate(points; delete_ghosts = false)
@test DelaunayTriangulation.has_ghost_triangles(tri) #src

# If you do need to query whether your triangulation already 
# has ghost triangles, use 
DelaunayTriangulation.has_ghost_triangles(tri)
@test DelaunayTriangulation.has_ghost_triangles(tri) #src

# To clear the ghost triangles, use
delete_ghost_triangles!(tri)
DelaunayTriangulation.has_ghost_triangles(tri)
@test !DelaunayTriangulation.has_ghost_triangles(tri) #src

# An important note for us to make is that ghost triangles 
# are not just there as a concept, but they are actually 
# physically stored.
add_ghost_triangles!(tri)
get_triangles(tri)

# See that there is not just the triangle `(1, 2, 3)`, but
# also `(3, 2, -1)`, `(1, 3, -1)`, and `(2, 1, -1)` (where 
# the ghost triangles are also oriented counter-clockwise). For example, 
get_adjacent(tri, 3, 2)
@test get_adjacent(tri, 3, 2) == -1 #src

#-
get_adjacent(tri, 3, -1)
@test get_adjacent(tri, 3, -1) == 1 #src

#-
get_adjacent(tri, -1, 2)
@test get_adjacent(tri, -1, 2) == 1 #src

#-
get_neighbours(tri, -1)

#-
get_adjacent2vertex(tri, -1)

# If we delete them, they are no longer there. 
delete_ghost_triangles!(tri)
get_triangles(tri)

# As a last note, we remark that the boundary vertices that defines the 
# vertex of these ghost triangles is still there regardless of whether 
# the triangulation has ghost triangles. Thus, for example, the 
# following still work 
get_neighbours(tri, -1)

#-
get_adjacent2vertex(tri, -1)

#-
get_adjacent(tri, 3, 2)
@test get_adjacent(tri, 3, 2) == -1 #src


# You can remove them from the graph, using 
DelaunayTriangulation.delete_boundary_vertices_from_graph!(tri)

# so that e.g. `get_neighbours(tri, -1)` is then an error. This 
# will still not remove them from the `adjacent` and `adjacent2vertex` 
# maps, but it does mean for example that 
collect(each_solid_vertex(tri)) == collect(each_vertex(tri))
