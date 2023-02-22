using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using LinearAlgebra
using StatsBase
using Random

include("../helper_functions.jl")

tri, label_map, index_map = simple_geometry()
add_ghost_triangles!(tri)
pts = get_points(tri)
adj = get_adjacent(tri)
adj2v = get_adjacent2vertex(tri)
boundary_index_ranges = get_boundary_index_ranges(tri)
boundary_nodes = get_boundary_nodes(tri)
boundary_map = get_boundary_map(tri)
DT.compute_representative_points!(tri)
DT.RepresentativePointList[1].x = 10.0
DT.RepresentativePointList[1].y = 10.0
_pts = tri.points[[12, 11, 10, 9]]
DT.RepresentativePointList[2].x = mean([8.0, 8.0, 4.0, 4.0])
DT.RepresentativePointList[2].y = mean([16.0, 6.0, 6.0, 16.0])
_pts = tri.points[[18, 17, 16, 15, 14, 13]]
DT.RepresentativePointList[3].x = mean([18.0, 18.0, 14.0, 12.0, 14.0, 14.0])
DT.RepresentativePointList[3].y = mean([12.0, 6.0, 2.0, 4.0, 6.0, 10.0])
graph = get_graph(tri)

k = index_map["a"]
q = (6.0, 0.0)
dir, dir_cert, id = DT.check_for_intersections_with_adjacent_boundary_edges(pts, adj,
                                                                            boundary_index_ranges,
                                                                            boundary_map, k,
                                                                            q)
@test dir == Certificate.Right && dir_cert == Certificate.On && id == index_map["b"]
@inferred DT.check_for_intersections_with_adjacent_boundary_edges(pts, adj,
                                                                  boundary_index_ranges,
                                                                  boundary_map, k, q)
q = (14.0, 0.0)
dir, dir_cert, id = DT.check_for_intersections_with_adjacent_boundary_edges(pts, adj,
                                                                            boundary_index_ranges,
                                                                            boundary_map, k,
                                                                            q)
@test dir == Certificate.Right && dir_cert == Certificate.Right && id == index_map["b"]
q = (0.0, 4.0)
dir, dir_cert, id = DT.check_for_intersections_with_adjacent_boundary_edges(pts, adj,
                                                                            boundary_index_ranges,
                                                                            boundary_map, k,
                                                                            q)
@test dir == Certificate.Left && dir_cert == Certificate.On && id == index_map["h"]
q = (0.0, 14.0)
dir, dir_cert, id = DT.check_for_intersections_with_adjacent_boundary_edges(pts, adj,
                                                                            boundary_index_ranges,
                                                                            boundary_map, k,
                                                                            q)
@test dir == Certificate.Left && dir_cert == Certificate.Left && id == index_map["h"]
q = (-2.0, 0.0)
dir, dir_cert, id = DT.check_for_intersections_with_adjacent_boundary_edges(pts, adj,
                                                                            boundary_index_ranges,
                                                                            boundary_map, k,
                                                                            q)
@test dir == Certificate.Outside && dir_cert == Certificate.Outside && id == k
q = (0.0, -2.0)
dir, dir_cert, id = DT.check_for_intersections_with_adjacent_boundary_edges(pts, adj,
                                                                            boundary_index_ranges,
                                                                            boundary_map, k,
                                                                            q)
@test dir == Certificate.Outside && dir_cert == Certificate.Outside && id == k
q = (2.0, -2.0)
dir, dir_cert, id = DT.check_for_intersections_with_adjacent_boundary_edges(pts, adj,
                                                                            boundary_index_ranges,
                                                                            boundary_map, k,
                                                                            q)
@test dir == Certificate.Outside && dir_cert == Certificate.Outside && id == k
q = (4.0, 4.0)
dir, dir_cert, id = DT.check_for_intersections_with_adjacent_boundary_edges(pts, adj,
                                                                            boundary_index_ranges,
                                                                            boundary_map, k,
                                                                            q)
@test dir == Certificate.Outside && dir_cert == Certificate.Outside && id == k
k = index_map["d"]
q = (20.0, 2.0)
dir, dir_cert, id = DT.check_for_intersections_with_adjacent_boundary_edges(pts, adj,
                                                                            boundary_index_ranges,
                                                                            boundary_map, k,
                                                                            q)
@test dir == Certificate.Left && dir_cert == Certificate.On && id == index_map["c"]
q = (20.0, 17.0)
dir, dir_cert, id = DT.check_for_intersections_with_adjacent_boundary_edges(pts, adj,
                                                                            boundary_index_ranges,
                                                                            boundary_map, k,
                                                                            q)
@test dir == Certificate.Right && dir_cert == Certificate.On && id == index_map["e"]
q = (20.0, -5.0)
dir, dir_cert, id = DT.check_for_intersections_with_adjacent_boundary_edges(pts, adj,
                                                                            boundary_index_ranges,
                                                                            boundary_map, k,
                                                                            q)
@test dir == Certificate.Left && dir_cert == Certificate.Left && id == index_map["c"]
q = (20.0, 30.0)
dir, dir_cert, id = DT.check_for_intersections_with_adjacent_boundary_edges(pts, adj,
                                                                            boundary_index_ranges,
                                                                            boundary_map, k,
                                                                            q)
@test dir == Certificate.Right && dir_cert == Certificate.Right && id == index_map["e"]
q = (10.0, 10.0)
dir, dir_cert, id = DT.check_for_intersections_with_adjacent_boundary_edges(pts, adj,
                                                                            boundary_index_ranges,
                                                                            boundary_map, k,
                                                                            q)
@inferred DT.check_for_intersections_with_adjacent_boundary_edges(pts, adj,
                                                                  boundary_index_ranges,
                                                                  boundary_map, k, q)
@test dir == Certificate.Outside && dir_cert == Certificate.Outside && id == k

k = index_map["g"]
right = index_map["h"]
left = index_map["f"]
q = (3.8544234210286, 19.2195032844691)
p = pts[k]
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
dir, dir_cert, id, rc, ℓc = DT.check_for_intersections_with_adjacent_boundary_edges(pts,
                                                                                    adj,
                                                                                    boundary_index_ranges,
                                                                                    boundary_map,
                                                                                    k, q)
@test rc == right_cert && ℓc == left_cert
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@inferred DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts, adj,
                                                                                   graph,
                                                                                   boundary_index_ranges,
                                                                                   boundary_map,
                                                                                   k, q,
                                                                                   right_cert,
                                                                                   left_cert)
@test (i, j, k) == (index_map["a1"], index_map["f"], index_map["g"]) &&
      DT.is_none(edge_cert) && DT.is_inside(tri_cert) && k == k
q = (6.0, 18.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["f"], index_map["a1"], index_map["g"]) &&
      DT.is_single(edge_cert) && DT.is_outside(tri_cert) && k == k
q = (4.2913632954055, 16.864882850327)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["a1"], index_map["i"], index_map["g"]) &&
      DT.is_single(edge_cert) && DT.is_outside(tri_cert) && k == k
q = (3.368934671721, 17.3989204745654)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["i"], index_map["a1"], index_map["g"]) &&
      DT.is_none(edge_cert) && DT.is_inside(tri_cert) && k == k
q = (4.0, 17.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["i"], index_map["a1"], index_map["g"]) &&
      DT.is_on(edge_cert) && DT.is_inside(tri_cert) && k == k
q = (5.0, 15.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@inferred DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts, adj,
                                                                                   graph,
                                                                                   boundary_index_ranges,
                                                                                   boundary_map,
                                                                                   k, q,
                                                                                   right_cert,
                                                                                   left_cert)
@test (j, i, k) == (index_map["b1"], index_map["i"], index_map["g"]) &&
      DT.is_right(edge_cert) && DT.is_outside(tri_cert)
@test DT.is_positively_oriented(DT.triangle_orientation(tri, i, j, q))
q = (5.0, 12.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["i"], index_map["b1"], index_map["g"]) &&
      DT.is_single(edge_cert) && DT.is_outside(tri_cert)
q = (3.0, 12.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["i"], index_map["b1"], index_map["g"]) &&
      DT.is_single(edge_cert) && DT.is_outside(tri_cert)
q = (3.0, 14.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["b1"], index_map["i"], index_map["g"]) &&
      DT.is_on(edge_cert) && DT.is_inside(tri_cert)
q = (2.0, 14.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["b1"], index_map["i"], index_map["g"]) &&
      DT.is_none(edge_cert) && DT.is_inside(tri_cert)
q = (0.5, 14.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["h"], index_map["b1"], index_map["g"]) &&
      DT.is_none(edge_cert) && DT.is_inside(tri_cert)
q = (1.0, 16.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["h"], index_map["b1"], index_map["g"]) &&
      DT.is_on(edge_cert) && DT.is_inside(tri_cert)
q = (1.0, 19.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["b1"], index_map["i"], index_map["g"]) &&
      DT.is_on(edge_cert) && DT.is_inside(tri_cert)
q = (2.0, 19.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["i"], index_map["a1"], index_map["g"]) &&
      DT.is_on(edge_cert) && DT.is_inside(tri_cert)
q = (3.0, 19.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["a1"], index_map["f"], index_map["g"]) &&
      DT.is_none(edge_cert) && DT.is_inside(tri_cert)
q = (18.0, 19.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["f"], index_map["a1"], index_map["g"]) &&
      DT.is_single(edge_cert) && DT.is_outside(tri_cert)
q = (4.0, 19.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@inferred DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts, adj,
                                                                                   graph,
                                                                                   boundary_index_ranges,
                                                                                   boundary_map,
                                                                                   k, q,
                                                                                   right_cert,
                                                                                   left_cert)
@test (i, j, k) == (index_map["a1"], index_map["f"], index_map["g"]) &&
      DT.is_none(edge_cert) && DT.is_inside(tri_cert)
q = (4.0, 12.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["i"], index_map["b1"], index_map["g"]) &&
      DT.is_single(edge_cert) && DT.is_outside(tri_cert)
q = (-1.0, 19.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (0, 0, index_map["g"]) &&
      DT.is_none(edge_cert) && DT.is_outside(tri_cert)
q = (-1.0, 20.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (0, 0, index_map["g"]) &&
      DT.is_none(edge_cert) && DT.is_outside(tri_cert)
q = (0.0, 21.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (0, 0, index_map["g"]) &&
      DT.is_none(edge_cert) && DT.is_outside(tri_cert)
q = (-1.0, 21.0)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (0, 0, index_map["g"]) &&
      DT.is_none(edge_cert) && DT.is_outside(tri_cert)
q = (3.0, 21.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (0, 0, index_map["g"]) &&
      DT.is_none(edge_cert) && DT.is_outside(tri_cert)
q = (-2.0, -2.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (0, 0, index_map["g"]) &&
      DT.is_none(edge_cert) && DT.is_outside(tri_cert)
k = index_map["b"]
right = index_map["c"]
left = index_map["a"]
q = (10.0, 1.0)
p = pts[k]
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["p"], index_map["t"], index_map["b"]) &&
      DT.is_on(edge_cert) && DT.is_inside(tri_cert)
q = (16.27931, 0.4487)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["c"], index_map["p"], index_map["b"]) &&
      DT.is_none(edge_cert) && DT.is_inside(tri_cert)
q = (17.0, 1.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["c"], index_map["p"], index_map["b"]) &&
      DT.is_on(edge_cert) && DT.is_inside(tri_cert)
q = (15.0, -1.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (0, 0, index_map["b"]) &&
      DT.is_none(edge_cert) && DT.is_outside(tri_cert)
q = (23.0, -1.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (0, 0, index_map["b"]) &&
      DT.is_none(edge_cert) && DT.is_outside(tri_cert)
q = (10.0, -1.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (0, 0, index_map["b"]) &&
      DT.is_none(edge_cert) && DT.is_outside(tri_cert)
q = (8.0, -1.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (0, 0, index_map["b"]) &&
      DT.is_none(edge_cert) && DT.is_outside(tri_cert)
q = (18.0, 3.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["p"], index_map["c"], index_map["b"]) &&
      DT.is_single(edge_cert) && DT.is_outside(tri_cert)
q = (11.4965, 2.53)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["t"], index_map["p"], index_map["b"]) &&
      DT.is_single(edge_cert) && DT.is_outside(tri_cert)
q = (5.0, 2.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["a"], index_map["t"], index_map["b"]) &&
      DT.is_single(edge_cert) && DT.is_outside(tri_cert)
q = (4.0, 1.0)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["a"], index_map["t"], index_map["b"]) &&
      DT.is_single(edge_cert) && DT.is_outside(tri_cert)
q = (5.2, 0.6)
right_cert = DT.point_position_relative_to_line(p, pts[right], q)
left_cert = DT.point_position_relative_to_line(p, pts[left], q)
i, j, edge_cert, tri_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                     adj,
                                                                                                     graph,
                                                                                                     boundary_index_ranges,
                                                                                                     boundary_map,
                                                                                                     k,
                                                                                                     q,
                                                                                                     right_cert,
                                                                                                     left_cert)
@test (i, j, k) == (index_map["t"], index_map["a"], index_map["b"]) &&
      DT.is_none(edge_cert) && DT.is_inside(tri_cert)
@inferred DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts, adj,
                                                                                   graph,
                                                                                   boundary_index_ranges,
                                                                                   boundary_map,
                                                                                   k, q,
                                                                                   right_cert,
                                                                                   left_cert)

k = index_map["a"]
q = (6.0, 0.0)
direction, q_pos, next_vertex = DT.check_for_intersections_with_adjacent_boundary_edges(pts,
                                                                                        adj,
                                                                                        boundary_index_ranges,
                                                                                        boundary_map,
                                                                                        k,
                                                                                        q)
q_cert, u, v, w = DT.search_down_adjacent_boundary_edges(pts, adj, boundary_index_ranges,
                                                         boundary_map, k, q, direction,
                                                         q_pos, next_vertex)
@test DT.is_on(q_cert) && (u, v, w) == (index_map["a"], index_map["b"], index_map["t"])
q = (10.0, 0.0)
direction, q_pos, next_vertex = DT.check_for_intersections_with_adjacent_boundary_edges(pts,
                                                                                        adj,
                                                                                        boundary_index_ranges,
                                                                                        boundary_map,
                                                                                        k,
                                                                                        q)
@test DT.is_degenerate(q_pos)
direction, q_pos, next_vertex = DT.check_for_intersections_with_adjacent_boundary_edges(pts,
                                                                                        adj,
                                                                                        boundary_index_ranges,
                                                                                        boundary_map,
                                                                                        k,
                                                                                        q)
q_cert, u, v, w = DT.search_down_adjacent_boundary_edges(pts, adj, boundary_index_ranges,
                                                         boundary_map, k, q, direction,
                                                         q_pos, next_vertex)
@test DT.is_on(q_cert) && (u, v, w) == (index_map["a"], index_map["b"], index_map["t"])
q = (12.0, 0.0)
direction, q_pos, next_vertex = DT.check_for_intersections_with_adjacent_boundary_edges(pts,
                                                                                        adj,
                                                                                        boundary_index_ranges,
                                                                                        boundary_map,
                                                                                        k,
                                                                                        q)
q_cert, u, v, w = DT.search_down_adjacent_boundary_edges(pts, adj, boundary_index_ranges,
                                                         boundary_map, k, q, direction,
                                                         q_pos, next_vertex)
@inferred DT.search_down_adjacent_boundary_edges(pts, adj, boundary_index_ranges,
                                                 boundary_map, k, q, direction, q_pos,
                                                 next_vertex)
@test DT.is_on(q_cert) && (u, v, w) == (index_map["b"], index_map["c"], index_map["p"])
q = (20.0, 0.0)
direction, q_pos, next_vertex = DT.check_for_intersections_with_adjacent_boundary_edges(pts,
                                                                                        adj,
                                                                                        boundary_index_ranges,
                                                                                        boundary_map,
                                                                                        k,
                                                                                        q)
q_cert, u, v, w = DT.search_down_adjacent_boundary_edges(pts, adj, boundary_index_ranges,
                                                         boundary_map, k, q, direction,
                                                         q_pos, next_vertex)
@test DT.is_on(q_cert) && (u, v, w) == (index_map["b"], index_map["c"], index_map["p"])
q = (22.0, 0.0)
direction, q_pos, next_vertex = DT.check_for_intersections_with_adjacent_boundary_edges(pts,
                                                                                        adj,
                                                                                        boundary_index_ranges,
                                                                                        boundary_map,
                                                                                        k,
                                                                                        q)
q_cert, u, v, w = DT.search_down_adjacent_boundary_edges(pts, adj, boundary_index_ranges,
                                                         boundary_map, k, q, direction,
                                                         q_pos, next_vertex)
@test DT.is_outside(q_cert) &&
      (u, v, w) == (index_map["d"], index_map["c"], DT.BoundaryIndex)
q = (0.0, 2.0)
direction, q_pos, next_vertex = DT.check_for_intersections_with_adjacent_boundary_edges(pts,
                                                                                        adj,
                                                                                        boundary_index_ranges,
                                                                                        boundary_map,
                                                                                        k,
                                                                                        q)
q_cert, u, v, w = DT.search_down_adjacent_boundary_edges(pts, adj, boundary_index_ranges,
                                                         boundary_map, k, q, direction,
                                                         q_pos, next_vertex)
@test DT.is_on(q_cert) && (u, v, w) == (index_map["h"], index_map["a"], index_map["s"])
q = (0.0, 10.0)
direction, q_pos, next_vertex = DT.check_for_intersections_with_adjacent_boundary_edges(pts,
                                                                                        adj,
                                                                                        boundary_index_ranges,
                                                                                        boundary_map,
                                                                                        k,
                                                                                        q)
q_cert, u, v, w = DT.search_down_adjacent_boundary_edges(pts, adj, boundary_index_ranges,
                                                         boundary_map, k, q, direction,
                                                         q_pos, next_vertex)
@test DT.is_on(q_cert) && (u, v, w) == (index_map["h"], index_map["a"], index_map["s"])
q = (0.0, 12.0)
direction, q_pos, next_vertex = DT.check_for_intersections_with_adjacent_boundary_edges(pts,
                                                                                        adj,
                                                                                        boundary_index_ranges,
                                                                                        boundary_map,
                                                                                        k,
                                                                                        q)
q_cert, u, v, w = DT.search_down_adjacent_boundary_edges(pts, adj, boundary_index_ranges,
                                                         boundary_map, k, q, direction,
                                                         q_pos, next_vertex)
@test DT.is_on(q_cert) && (u, v, w) == (index_map["g"], index_map["h"], index_map["b1"])
q = (0.0, 20.0)
direction, q_pos, next_vertex = DT.check_for_intersections_with_adjacent_boundary_edges(pts,
                                                                                        adj,
                                                                                        boundary_index_ranges,
                                                                                        boundary_map,
                                                                                        k,
                                                                                        q)
q_cert, u, v, w = DT.search_down_adjacent_boundary_edges(pts, adj, boundary_index_ranges,
                                                         boundary_map, k, q, direction,
                                                         q_pos, next_vertex)
@test DT.is_on(q_cert) && (u, v, w) == (index_map["g"], index_map["h"], index_map["b1"])
q = (0.0, 22.0)
direction, q_pos, next_vertex = DT.check_for_intersections_with_adjacent_boundary_edges(pts,
                                                                                        adj,
                                                                                        boundary_index_ranges,
                                                                                        boundary_map,
                                                                                        k,
                                                                                        q)
q_cert, u, v, w = DT.search_down_adjacent_boundary_edges(pts, adj, boundary_index_ranges,
                                                         boundary_map, k, q, direction,
                                                         q_pos, next_vertex)
@test DT.is_outside(q_cert) &&
      (u, v, w) == (index_map["g"], index_map["f"], DT.BoundaryIndex)
@inferred DT.search_down_adjacent_boundary_edges(pts, adj, boundary_index_ranges,
                                                 boundary_map, k, q, direction, q_pos,
                                                 next_vertex)
q = (0.0, 0.0)
direction, q_pos, next_vertex = DT.check_for_intersections_with_adjacent_boundary_edges(pts,
                                                                                        adj,
                                                                                        boundary_index_ranges,
                                                                                        boundary_map,
                                                                                        k,
                                                                                        q)
q_cert, u, v, w = DT.search_down_adjacent_boundary_edges(pts, adj, boundary_index_ranges,
                                                         boundary_map, k, q, direction,
                                                         q_pos, next_vertex)
@test DT.is_on(q_cert) && (u, v, w) == (index_map["a"], index_map["b"], index_map["t"])

k = index_map["b"]
q = [2.0, 4.0]
direction, q_pos, next_vertex, right_cert, left_cert = DT.check_for_intersections_with_adjacent_boundary_edges(pts,
                                                                                                               adj,
                                                                                                               boundary_index_ranges,
                                                                                                               boundary_map,
                                                                                                               k,
                                                                                                               q)
@test DT.is_outside(direction)
i, j, edge_cert, triangle_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                          adj,
                                                                                                          graph,
                                                                                                          boundary_index_ranges,
                                                                                                          boundary_map,
                                                                                                          k,
                                                                                                          q,
                                                                                                          right_cert,
                                                                                                          left_cert)
@test i == index_map["a"] &&
      j == index_map["t"] &&
      DT.has_one_intersection(edge_cert) &&
      DT.is_outside(triangle_cert)

Random.seed!(191919)
q = (2.0, 4.0)
k = 6
i, j, edge_cert, triangle_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                          adj,
                                                                                                          graph,
                                                                                                          boundary_index_ranges,
                                                                                                          boundary_map,
                                                                                                          k,
                                                                                                          q,
                                                                                                          right_cert,
                                                                                                          left_cert)

for k in each_point_index(pts)
    local dir, rc, ℓc
    if DT.is_outer_boundary_node(tri, k)
        dir, qp, _k, rc, ℓc = DT.check_for_intersections_with_adjacent_boundary_edges(tri.points,
                                                                                      tri.adjacent,
                                                                                      tri.boundary_index_ranges,
                                                                                      boundary_map,
                                                                                      k,
                                                                                      get_point(tri.points,
                                                                                                boundary_map,
                                                                                                k))
        @test DT.is_degenerate(qp)
    end
end

k = 6
q = get_point(pts, 9)
direction, q_pos, next_vertex, right_cert, left_cert = DT.check_for_intersections_with_adjacent_boundary_edges(pts,
                                                                                                               adj,
                                                                                                               boundary_index_ranges,
                                                                                                               boundary_map,
                                                                                                               k,
                                                                                                               q)
i, j, edge_cert, triangle_cert = DT.check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
                                                                                                          adj,
                                                                                                          graph,
                                                                                                          boundary_index_ranges,
                                                                                                          boundary_map,
                                                                                                          k,
                                                                                                          q,
                                                                                                          right_cert,
                                                                                                          left_cert)
@test i == 25 && j == 9 && DT.is_on(edge_cert) && DT.is_inside(triangle_cert)
