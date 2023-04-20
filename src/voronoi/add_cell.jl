function add_voronoi_cell!(vorn::VoronoiTessellation, i)
    I = integer_type(vorn)
    T = triangle_type(vorn)
    S = get_surrounding_polygon(vorn, i)
    B = VoronoiCell{I}()
    boundary_cells = get_boundary_cells(vorn)
    unbounded_cells = get_unbounded_cells(vorn)
    # Find the first triangle to start in
    j = S[begin]
    m = firstindex(S) + 1
    k = S[m]
    V = (rotate_triangle_to_standard_form ∘ construct_triangle)(T, i, j, k)
    ci = get_triangle_to_circumcenter(vorn, V)
    # cx, cy = get_point(cell_points, ci)
    add_vertex!(B, ci)
    # Now go over to the next triangles
    for m in (firstindex(S)+2):lastindex(S)
        j = k
        k = S[m]
        V = (rotate_triangle_to_standard_form ∘ construct_triangle)(T, i, j, k)
        ci = get_triangle_to_circumcenter(vorn, V)
        if is_boundary_index(ci)
            push!(boundary_cells, i)
        end
        # cx, cy = get_point(cell_points, ci)
        add_vertex!(B, ci)
        if is_boundary_index(ci)
            push!(unbounded_cells, i) # do this here in case it got converted previously
        end
    end
    add_vertex!(B, get_vertices(B, 1))
    cells = get_cells(vorn)
    cells[i] = B
    return B
end