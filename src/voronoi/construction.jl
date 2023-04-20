function voronoi(tri::Triangulation, add_topology=true)
    has_ghost = has_ghost_triangles(tri)
    !has_ghost && add_ghost_triangles!(tri)
    vorn = initialise_voronoi_tessellation(tri)
    for i in each_generator(vorn)
        add_voronoi_cell!(vorn, i)
    end
    add_topology && add_topological_information!(vorn)
    !has_ghost && delete_ghost_triangles!(tri)
    return vorn
end