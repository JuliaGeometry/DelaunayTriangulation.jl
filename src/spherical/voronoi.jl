const SphericalTessellation = VoronoiTessellation{<:SphericalTriangulation}

function spherical_voronoi(tri::SphericalTriangulation; kwargs...)
    vorn = voronoi(tri; kwargs...)
    circumcenter_to_triangle = get_circumcenter_to_triangle(vorn)
    triangle_to_circumcenter = get_triangle_to_circumcenter(vorn)
    polygon_points = get_polygon_points(vorn)
    I = integer_type(tri)
    mapped_ghosts = Dict{I, I}()
    for polygon_vertices in each_polygon(vorn)
        for (i, v) in pairs(polygon_vertices) 
            is_ghost_vertex(v) || continue 
            if haskey(mapped_ghosts, v)
                polygon_vertices[i] = mapped_ghosts[v]
            else
                T = get_circumcenter_to_triangle(vorn, v)
                c = triangle_circumcenter(tri, T)
                pidx = findfirst(==(c), polygon_points)
                if isnothing(pidx)
                    push_polygon_point!(vorn, c)
                    n = num_polygon_vertices(vorn)
                else
                    n = pidx
                end
                triangle_to_circumcenter[T] = n
                mapped_ghosts[v] = n
                polygon_vertices[i] = n
            end
        end
    end
    return vorn
end
