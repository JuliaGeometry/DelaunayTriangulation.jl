
function initialise_voronoi_tessellation(tri::Tr) where {Tr<:Triangulation}
    ## Setup
    I = integer_type(tri)
    T = triangle_type(tri)
    F = number_type(tri)
    P = NTuple{2,F}
    constrained_edges = get_all_constrained_edges(tri)
    Es = typeof(constrained_edges)
    E = edge_type(tri)
    cell_points = Vector{P}()
    circumcenter_to_triangle = Dict{I,T}()
    triangle_to_circumcenter = Dict{T,I}()
    sizehint!(cell_points, num_triangles(tri) + num_edges(constrained_edges))
    sizehint!(circumcenter_to_triangle, num_triangles(tri))
    sizehint!(triangle_to_circumcenter, num_triangles(tri))
    cur_ghost_idx = I(0)

    ## Build the circumcenter map
    for V in each_triangle(tri)
        V = rotate_triangle_to_standard_form(V)
        if !is_ghost_triangle(V)
            cx, cy = triangle_circumcenter(tri, V)
            push_point!(cell_points, cx, cy)
            circumcenter_to_triangle[num_points(cell_points)] = V
            triangle_to_circumcenter[V] = num_points(cell_points)
        else
            circumcenter_to_triangle[I(BoundaryIndex)-cur_ghost_idx] = V
            triangle_to_circumcenter[V] = I(BoundaryIndex) - cur_ghost_idx
            cur_ghost_idx += I(1)
        end
    end

    ## Add in the obstacles 
    obstacle_mapping = Dict{I,I}()
    obstacles = initialise_edges(Es)
    sizehint!(obstacle_mapping, num_edges(constrained_edges))
    sizehint!(obstacles, num_edges(constrained_edges))
    for e in each_edge(constrained_edges)
        u, v = edge_indices(e)
        if u ∉ keys(obstacle_mapping)
            p = get_point(tri, u)
            push_point!(cell_points, p)
            obstacle_mapping[u] = num_points(cell_points)
        end
        if v ∉ keys(obstacle_mapping)
            p = get_point(tri, v)
            push_point!(cell_points, p)
            obstacle_mapping[v] = num_points(cell_points)
        end
        i = obstacle_mapping[u]
        j = obstacle_mapping[v]
        e = construct_edge(E, i, j)
        add_edge!(obstacles, e)
    end

    ## Get the generators
    generators = Vector{P}()
    sizehint!(generators, num_solid_vertices(tri))
    for i in each_point_index(tri) # This can even include points that aren't in the triangulation, but still in the point set. It's just useful to keep the indices here matching up with the original triangulation. Note also that we are avoiding future mutations of points by just copying
        p = get_point(tri, i)
        push_point!(generators, p)
    end

    ## Initialise the cells 
    cells = Vector{VoronoiCell{I}}(undef, num_points(generators))

    ## Initialise the adjacent maps
    adjacent_cell = Adjacent{I,E}()
    adjacent_vertex = Adjacent{I,E}()

    ## Initialise the Adjacent2Vertex map 
    vertex_to_cells = Adjacent2Vertex{I,Set{I},I}()

    ## Initialise the graphs 
    vertex_to_vertices = Graph{I}()
    cell_to_cells = Graph{I}()

    ## Prepare the boundary cells 
    boundary_cells = Set{I}()
    unbounded_cells = Set{I}()

    ## Return the Voronoi tessellation
    return VoronoiTessellation{P,I,E,Es,T,Tr}(generators, cell_points, cells,
        adjacent_cell, adjacent_vertex, vertex_to_cells, vertex_to_vertices,
        cell_to_cells, obstacles, boundary_cells, unbounded_cells,
        circumcenter_to_triangle, triangle_to_circumcenter, tri)
end
