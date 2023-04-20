function add_topological_information!(vorn::VoronoiTessellation)
    adj_cell = get_adjacent_cell(vorn)
    adj_vertex = get_adjacent_vertex(vorn)
    vert_to_cells = get_vertex_to_cells(vorn)
    vert_to_verts = get_vertex_to_vertices(vorn)
    for i in each_cell(vorn)
        C = get_cell(vorn, i)
        v = get_vertices(C)
        ne = num_boundary_edges(v)
        for j in 1:ne 
            a = get_boundary_nodes(v, j)
            b = get_boundary_nodes(v, j+1)
            c = nextindex_circular(v, b) # get_boundary_nodes(v, j+2)
            add_adjacent!(adj_cell, a, b, i)
            add_adjacent!(adj_vertex, a, b, c)
            adj2v_cells = get!(vert_to_cells, a)
            push!(adj2v_cells, a)
            add_neighbour!(vert_to_verts, a, b)
        end
    end
    cell_to_cells = deepcopy(get_graph(get_triangulation(vorn)))
    delete_boundary_vertices_from_graph!(cell_to_cells)
    return nothing
end