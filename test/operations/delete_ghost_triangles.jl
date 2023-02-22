using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using StructEquality

@struct_equal Triangulation
@struct_equal DT.Graph
@struct_equal DT.Adjacent2Vertex
@struct_equal DT.Adjacent

include("../helper_functions.jl")

tri, label_map, index_map = simple_geometry()
_tri = deepcopy(tri)
DT.add_ghost_triangles!(tri)
DT.delete_ghost_triangles!(tri)

@test tri == _tri
