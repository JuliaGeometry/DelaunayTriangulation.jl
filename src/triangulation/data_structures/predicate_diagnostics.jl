mutable struct PredicateDiagnostics
    orient2::Int
    orient3::Int
    incircle::Int
    parallelorder::Int
    angle_is_acute::Int
    sameside::Int
    meet::Int
    triangle_orientation::Int
    point_position_relative_to_circle::Int
    point_position_relative_to_line::Int
    point_closest_to_line::Int
    point_position_on_line_segment::Int
    line_segment_intersection_type::Int
    point_position_relative_to_triangle::Int
    point_position_relative_to_oriented_outer_halfplane::Int
    is_legal::Int
    triangle_line_segment_intersection::Int
    opposite_angle::Int
    point_position_relative_to_diametral_circle::Int
    point_position_relative_to_diametral_lens::Int
    find_edge::Int
    point_position_relative_to_circumcircle::Int
    point_position_relative_to_witness_plane::Int
    test_visiblity::Int
end
PredicateDiagnostics() = PredicateDiagnostics(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

for f in fieldnames(PredicateDiagnostics)
    addf = Symbol("add_", f, "!")
    numf = Symbol("num_", f)
    @eval begin 
        @inline $addf(pd::PredicateDiagnostics) = (pd.$f += 1)
        @inline $numf(pd::PredicateDiagnostics) = pd.$f
        @inline $addf(pd::Nothing) = pd
    end
end
