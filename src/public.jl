@static if VERSION >= v"1.11.0-DEV.469"
    eval(Meta.parse(
        """
        public Adjacent,
            Adjacent2Vertex,
            ConvexHull,
            Graph,
            InsertionEventHistory, 
            Triangulation,
            VoronoiTessellation,
            ZeroWeight,
            add_weight!,
            all_ghost_vertices,
            check_args,
            compute_representative_points!,
            contains_boundary_edge,
            convert_boundary_points_to_indices,
            delete_ghost_vertices_from_graph!,
            dist,
            get_adjacent,
            get_adjacent2vertex,
            get_all_boundary_nodes,
            get_all_segments,
            get_area,
            get_boundary_curves,
            get_boundary_edge_map,
            get_boundary_nodes,
            get_convex_hull,
            get_convex_hull_vertices,
            get_curve_index,
            get_edges,
            get_ghost_vertex_map,
            get_ghost_vertex_range,
            get_ghost_vertex_ranges,
            get_graph,
            get_insertion_order,
            get_interior_segments,
            get_left_boundary_node,
            get_neighbours,
            get_point,
            get_points,
            get_representative_point_coordinates,
            get_representative_point_list,
            get_right_boundary_node,
            get_section_index,
            get_triangles,
            get_vertices,
            get_weight,
            get_weights,
            has_ghost_vertices,
            has_vertex,
            is_exterior_ghost_vertex,
            is_interior_ghost_vertex,
            iterated_neighbourhood,
            map_ghost_vertex,
            num_curves,
            num_neighbours,
            num_points,
            num_sections,
            refine!,
            retriangulate,
            triangulate,
            triangulate_convex,
            triangulate_rectangle,
            validate_triangulation,
            add_boundary_information!,
            add_ghost_triangles!,
            add_point!,
            add_segment!,
            add_triangle!,
            clear_empty_features!,
            complete_split_edge_and_legalise!,
            complete_split_triangle_and_legalise!,
            delete_ghost_triangles!,
            delete_holes!,
            delete_point!,
            delete_triangle!,
            flip_edge!,
            get_surrounding_polygon,
            legalise_edge!,
            lock_convex_hull!,
            split_edge!,
            split_triangle!,
            unlock_convex_hull!,
            centroidal_smooth,
            get_boundary_polygons,
            get_centroid,
            get_circumcenter_to_triangle,
            get_cocircular_circumcenters,
            get_generator,
            get_generators,
            get_polygon,
            get_polygon_coordinates,
            get_polygon_point,
            get_polygon_points,
            get_polygons,
            get_triangle_to_circumcenter,
            get_triangulation,
            get_unbounded_polygons,
            num_generators,
            num_polygon_vertices,
            num_polygons,
            polygon_bounds,
            polygon_features,
            toggle_inf_warn!,
            voronoi,
            convex_hull,
            convex_hull!,
            AbstractParametricCurve,
            BSpline,
            BezierCurve,
            CatmullRomSpline,
            CircularArc,
            EllipticalArc,
            LineSegment,
            angle_between,
            arc_length,
            curvature,
            differentiate,
            get_circle_intersection,
            get_closest_point,
            get_equidistant_split,
            get_equivariation_split,
            get_inverse,
            marked_total_variation,
            orientation_markers,
            point_position_relative_to_curve,
            thrice_differentiate,
            total_variation,
            twice_differentiate,
            each_boundary_edge,
            each_edge,
            each_generator,
            each_ghost_edge,
            each_ghost_triangle,
            each_ghost_vertex,
            each_point,
            each_point_index,
            each_polygon,
            each_polygon_index,
            each_polygon_vertex,
            each_segment,
            each_solid_edge,
            each_solid_triangle,
            each_solid_vertex,
            each_triangle,
            each_unbounded_polygon,
            each_vertex,
            num_edges,
            num_ghost_edges,
            num_ghost_triangles,
            num_ghost_vertices,
            num_solid_edges,
            num_solid_triangles,
            num_solid_vertices,
            num_triangles,
            num_vertices,
            brute_force_search,
            find_polygon,
            find_triangle,
            get_nearest_neighbour,
            Certificate,
            contains_segment,
            contains_triangle,
            edge_exists,
            find_edge,
            has_boundary_nodes,
            has_ghost_triangles,
            has_multiple_curves,
            has_multiple_intersections,
            has_multiple_sections,
            has_no_intersections,
            has_one_intersection,
            is_above,
            is_acute,
            is_below,
            is_boundary_edge,
            is_boundary_node,
            is_boundary_triangle,
            is_closer,
            is_collinear,
            is_constrained,
            is_degenerate,
            is_equidistant,
            is_further,
            is_ghost_edge,
            is_ghost_triangle,
            is_ghost_vertex,
            is_illegal,
            is_inside,
            is_left,
            is_legal,
            is_multiple,
            is_negatively_oriented,
            is_negativelyoriented,
            is_none,
            is_obtuse,
            is_on,
            is_outside,
            is_positively_oriented,
            is_positivelyoriented,
            is_right,
            is_single,
            is_touching,
            is_weighted,
            line_segment_intersection_type,
            opposite_angle,
            point_closest_to_line,
            point_position_on_line_segment,
            point_position_relative_to_circle,
            point_position_relative_to_circumcircle,
            point_position_relative_to_diametral_circle,
            point_position_relative_to_diametral_lens,
            point_position_relative_to_line,
            point_position_relative_to_oriented_outer_halfplane,
            point_position_relative_to_triangle,
            point_position_relative_to_witness_plane,
            triangle_line_segment_intersection,
            triangle_orientation,
            unoriented_edge_exists,
            IndividualTriangleStatistics,
            TriangulationStatistics,
            get_all_stat,
            get_angles,
            get_aspect_ratio,
            get_circumcenter,
            get_circumradius,
            get_edge_midpoints,
            get_individual_statistics,
            get_inradius,
            get_largest_angle,
            get_largest_area,
            get_largest_radius_edge_ratio,
            get_lengths,
            get_maximum_angle,
            get_median_angle,
            get_minimum_angle,
            get_offcenter,
            get_perimeter,
            get_radius_edge_ratio,
            get_sink,
            get_smallest_angle,
            get_smallest_area,
            get_smallest_radius_edge_ratio,
            num_boundary_segments,
            num_convex_hull_vertices,
            num_interior_segments,
            num_segments,
            statistics,
            triangle_angles,
            triangle_area,
            triangle_aspect_ratio,
            triangle_centroid,
            triangle_circumcenter,
            triangle_circumradius,
            triangle_edge_midpoints,
            triangle_inradius,
            triangle_lengths,
            triangle_offcenter,
            triangle_perimeter,
            triangle_radius_edge_ratio,
            triangle_sink,
            clip_polygon,
            construct_polygon_hierarchy,
            distance_to_polygon,
            number_type,
            pole_of_inaccessibility,
            is_point2
    """))
end