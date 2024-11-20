mutable struct TriangleCounts
    num_solid_triangles::Int
    num_ghost_triangles::Int
    num_solid_edges::Int
    num_ghost_edges::Int
    num_solid_vertices::Int
    num_ghost_vertices::Int
end
TriangleCounts() = TriangleCounts(0, 0, 0, 0, 0, 0)

@inline add_solid_triangle!(tc::TriangleCounts) = (tc.num_solid_triangles += 1)
@inline add_ghost_triangle!(tc::TriangleCounts) = (tc.num_ghost_triangles += 1)
@inline add_solid_edge!(tc::TriangleCounts) = (tc.num_solid_edges += 1)
@inline add_ghost_edge!(tc::TriangleCounts) = (tc.num_ghost_edges += 1)
@inline add_solid_vertex!(tc::TriangleCounts) = (tc.num_solid_vertices += 1)
@inline add_ghost_vertex!(tc::TriangleCounts) = (tc.num_ghost_vertices += 1)

@inline num_solid_triangles(tc::TriangleCounts) = tc.num_solid_triangles
@inline num_ghost_triangles(tc::TriangleCounts) = tc.num_ghost_triangles
@inline num_triangles(tc::TriangleCounts) = num_solid_triangles(tc) + num_ghost_triangles(tc)
@inline num_solid_edges(tc::TriangleCounts) = tc.num_solid_edges
@inline num_ghost_edges(tc::TriangleCounts) = tc.num_ghost_edges
@inline num_edges(tc::TriangleCounts) = num_solid_edges(tc) + num_ghost_edges(tc)
@inline num_solid_vertices(tc::TriangleCounts) = tc.num_solid_vertices
@inline num_ghost_vertices(tc::TriangleCounts) = tc.num_ghost_vertices
@inline num_vertices(tc::TriangleCounts) = num_solid_vertices(tc) + num_ghost_vertices(tc)
