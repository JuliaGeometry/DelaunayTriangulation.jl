function Triangulation(points::P;
    IntegerType::Type{I}=Int64,
    EdgeType::Type{E}=NTuple{2,IntegerType},
    TriangleType::Type{V}=NTuple{3,IntegerType},
    EdgesType::Type{Es}=Set{EdgeType},
    TrianglesType::Type{Ts}=Set{TriangleType},
    boundary_nodes::BN=IntegerType[],
    constrained_edges=initialise_edges(EdgesType),
    representative_point_list=get_empty_representative_points(IntegerType, number_type(points))) where {P,Ts,I,E,Es,BN,V}
    T = initialise_triangles(Ts)
    adj = Adjacent{I,E}()
    adj2v = Adjacent2Vertex{I,Es,E}()
    graph = Graph{I}()
    bnn_map = construct_boundary_edge_map(boundary_nodes; IntegerType=I, EdgeType=E)
    bn_map = construct_boundary_map(boundary_nodes; IntegerType=I)
    bn_range = construct_boundary_index_ranges(boundary_nodes; IntegerType=I)
    ch = ConvexHull(points, I[])
    all_constrained_edges = initialise_edges(EdgesType)
    n = num_points(points)
    sizehint!(T, 2n - 5) # maximum number of triangles
    sizehint!(adj, 3n - 6) # maximum number of edges 
    sizehint!(adj2v, n)
    sizehint!(graph, 3n - 6, n, n)
    sizehint!(ch, n)
    tri = Triangulation(points, T, adj, adj2v, graph, boundary_nodes, bnn_map, bn_map, bn_range,
        constrained_edges, all_constrained_edges, ch, representative_point_list)
    return tri
end

function Triangulation(points::P, triangles::T, boundary_nodes::BN;
    IntegerType::Type{I}=Int64,
    EdgeType::Type{E}=NTuple{2,IntegerType},
    TriangleType::Type{V}=NTuple{3,IntegerType},
    EdgesType::Type{Es}=Set{EdgeType},
    TrianglesType::Type{Ts}=Set{TriangleType},
    add_ghost_triangles=false) where {P,BN,I,E,V,Es,Ts,T}
    tri = Triangulation(points; boundary_nodes, IntegerType, EdgeType, TriangleType, EdgesType, TrianglesType)
    adj = get_adjacent(tri)
    adj2v = get_adjacent2vertex(tri)
    graph = get_graph(tri)
    tris = get_triangles(tri)
    for τ in each_triangle(triangles)
        add_triangle!(adj, τ)
        add_triangle!(adj2v, τ)
        add_triangle!(graph, τ)
        add_triangle!(tris, τ)
    end
    add_boundary_information!(tri)
    add_ghost_triangles && add_ghost_triangles!(tri)
    convex_hull!(tri; reconstruct=true)
    constrained_edges = get_all_constrained_edges(tri)
    bn_map = get_boundary_map(tri)
    all_edges = merge_constrained_edges(bn_map, boundary_nodes, initialise_edges(Es))
    for edge in each_edge(all_edges)
        add_edge!(constrained_edges, edge)
    end
    return tri
end

function remake_triangulation_with_constraints(tri::Triangulation{P,Ts,I,E,Es,BN,B,BIR}, edges, boundary_nodes) where {P,Ts,I,E,Es,BN,B,BIR}
    points = get_points(tri)
    triangles = get_triangles(tri)
    adjacent = get_adjacent(tri)
    adjacent2vertex = get_adjacent2vertex(tri)
    graph = get_graph(tri)
    boundary_nodes = @something boundary_nodes get_boundary_nodes(tri)
    boundary_edge_map = construct_boundary_edge_map(boundary_nodes; IntegerType=I, EdgeType=E)
    bn_map = construct_boundary_map(boundary_nodes; IntegerType=I)
    bn_range = construct_boundary_index_ranges(boundary_nodes; IntegerType=I)
    constrained_edges = @something edges get_constrained_edges(tri)
    all_constrained_edges = get_all_constrained_edges(tri)
    convex_hull = get_convex_hull(tri)
    representative_point_list = get_representative_point_list(tri)
    return bn_map, bn_range, Triangulation(
        points,
        triangles,
        adjacent,
        adjacent2vertex,
        graph,
        boundary_nodes,
        boundary_edge_map,
        get_boundary_map(tri), # Delay putting these in until we are done with the triangulation
        get_boundary_index_ranges(tri),
        constrained_edges,
        all_constrained_edges,
        convex_hull,
        representative_point_list)
end

function replace_boundary_dict_information(tri::Triangulation, bn_map, bn_range)
    points = get_points(tri)
    triangles = get_triangles(tri)
    adjacent = get_adjacent(tri)
    adjacent2vertex = get_adjacent2vertex(tri)
    graph = get_graph(tri)
    boundary_nodes = get_boundary_nodes(tri)
    constrained_edges = get_constrained_edges(tri)
    all_constrained_edges = get_all_constrained_edges(tri)
    bnn_map = get_boundary_edge_map(tri)
    convex_hull = get_convex_hull(tri)
    representative_point_list = get_representative_point_list(tri)
    return Triangulation(
        points,
        triangles,
        adjacent,
        adjacent2vertex,
        graph,
        boundary_nodes,
        bnn_map,
        bn_map,
        bn_range,
        constrained_edges,
        all_constrained_edges,
        convex_hull,
        representative_point_list
    )
end




