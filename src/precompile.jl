PrecompileTools.@setup_workload begin
    p1 = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    b1 = [1, 2, 3, 4, 1]

    p2 = [(0.25, 0.25), (0.25, 0.75), (0.75, 0.75), (0.75, 0.25)]
    p2 = [p1; p2]
    b2 = [[[1, 2, 3, 4, 1]], [[5, 6, 7, 8, 5]]]

    PrecompileTools.@compile_workload begin
        tri = triangulate(p1)
        vor = voronoi(tri, clip = true, smooth = true)
        refine!(tri)
        tri = triangulate(p1; boundary_nodes = b1)
        tri = triangulate(p2; boundary_nodes = b2)
        circ = CircularArc((1.0, 0.0), (1.0, 0.0), (0.0, 0.0))
        tri = triangulate(NTuple{2, Float64}[]; boundary_nodes = [[circ]])
        refine!(tri)
    end
end