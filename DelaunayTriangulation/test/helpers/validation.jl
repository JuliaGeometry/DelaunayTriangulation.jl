# Validation of statistics, tessellations, refinement and insertion event histories, and checks for the constrained triangulation tests.

export validate_refinement, validate_statistics, validate_tessellation, is_conformal
export test_segment_triangle_intersections, test_split_edges, test_intersections, validate_insertion_event_history
export slow_triangle_assess_queue, slow_triangle_assess, compare_encroach_queues, compare_triangle_queues

function test_intersections(tri, e, allT, constrained_edges)
    for e in ((e[1], e[2]), (e[2], e[1]))
        intersecting_triangles, collinear_segments, left, right = DT.locate_intersecting_triangles(tri, e)
        @test all(T -> DT.is_positively_oriented(DT.triangle_orientation(tri, T)), intersecting_triangles)
        @test all(!DT.is_none, [DT.triangle_line_segment_intersection(tri, T..., e...) for T in intersecting_triangles])
        @test DT.compare_triangle_collections(allT, intersecting_triangles)
        @test allunique(intersecting_triangles)
        if typeof(constrained_edges) <: AbstractVector
            @test collinear_segments == constrained_edges
        else # Tuple of possibilities, in case the edge's endpoints have equal degree so that we could start at any point
            @test any(==(collinear_segments), constrained_edges)
        end
    end
end

function test_split_edges(tri, edge, current_constrained_edges)
    constrained_edges = get_interior_segments(tri)
    DT.add_edge!(constrained_edges, edge)
    _, collinear_segments = DT.locate_intersecting_triangles(tri, edge)
    DT.split_segment!(tri, edge, collinear_segments)
    @test any(==(constrained_edges), current_constrained_edges)
    return nothing
end

function test_segment_triangle_intersections(tri, edge, true_triangles, true_collinear_segments, current_constrained_edges)
    constrained_edges = get_interior_segments(tri)
    for edge in ((edge[1], edge[2]), (edge[2], edge[1]))
        intersecting_triangles, collinear_segments, left, right = DT.locate_intersecting_triangles(tri, edge)
        @test all(T -> DT.is_positively_oriented(DT.triangle_orientation(tri, T)), intersecting_triangles)
        @test all(!DT.is_none, [DT.triangle_line_segment_intersection(tri, T..., edge...) for T in intersecting_triangles])
        if typeof(true_triangles) <: AbstractVector || typeof(true_triangles) <: AbstractSet
            @test DT.compare_triangle_collections(true_triangles, intersecting_triangles)
        else
            @test any(V -> DT.compare_triangle_collections(intersecting_triangles, V), true_triangles)
        end
        @test allunique(intersecting_triangles)
        if typeof(true_collinear_segments) <: AbstractVector || typeof(true_collinear_segments) <: AbstractSet
            @test collinear_segments == true_collinear_segments
        else # Tuple of possibilities, in case the edge's endpoints have equal degree so that we could start at any point
            @test any(==(collinear_segments), true_collinear_segments)
        end
    end
    constrained_edges = get_interior_segments(tri)
    DT.add_edge!(constrained_edges, edge)
    _, collinear_segments = DT.locate_intersecting_triangles(tri, edge)
    DT.split_segment!(tri, edge, collinear_segments)
    !isnothing(current_constrained_edges) && @test compare_edge_vectors(constrained_edges, current_constrained_edges)
end

function validate_statistics(tri::DT.Triangulation, stats=statistics(tri); msg="")
    ## Build up the array (see also test_iterators)
    I = DT.integer_type(tri)
    T = NTuple{3,I}
    E = NTuple{2,I}
    solid_triangles = T[]
    ghost_triangles = T[]
    all_triangles = T[]
    solid_vertices = Set{I}()
    ghost_vertices = Set{I}()
    all_vertices = Set{I}()
    solid_edges = E[]
    ghost_edges = E[]
    all_edges = E[]
    for T in each_triangle(tri)
        i, j, k = DT.triangle_vertices(T)
        push!(all_triangles, (i, j, k))
        if DT.is_ghost_triangle(i, j, k)
            push!(ghost_triangles, (i, j, k))
        else
            push!(solid_triangles, (i, j, k))
        end
        for (u, v) in DT.triangle_edges(i, j, k)
            push!(all_edges, (u, v))
            if DT.is_ghost_edge(u, v)
                push!(ghost_edges, (u, v))
            else
                push!(solid_edges, (u, v))
            end
        end
        for k in (i, j, k)
            push!(all_vertices, k)
            if DT.is_ghost_vertex(k)
                push!(ghost_vertices, k)
            else
                push!(solid_vertices, k)
            end
        end
    end
    push!(all_vertices, DT.all_ghost_vertices(tri)...)
    push!(ghost_vertices, DT.all_ghost_vertices(tri)...)
    for e in each_ghost_edge(tri)
        u, v = DT.edge_vertices(e)
        push!(ghost_edges, (min(u, v), max(u, v)))
        push!(all_edges, (min(u, v), max(u, v)))
    end
    for (i, e) in enumerate(all_edges)
        u, v = DT.edge_vertices(e)
        all_edges[i] = (min(u, v), max(u, v))
    end
    for (i, e) in enumerate(solid_edges)
        u, v = DT.edge_vertices(e)
        solid_edges[i] = (min(u, v), max(u, v))
    end
    for (i, e) in enumerate(ghost_edges)
        u, v = DT.edge_vertices(e)
        ghost_edges[i] = (min(u, v), max(u, v))
    end
    unique!(all_edges)
    unique!(solid_edges)
    unique!(ghost_edges)
    all_vertices = collect(all_vertices)
    solid_vertices = collect(solid_vertices)
    ghost_vertices = collect(ghost_vertices)
    sort!(all_vertices)
    sort!(solid_vertices)
    sort!(ghost_vertices)
    NT = DT.number_type(tri)
    ## Build up the individual statistics 
    areas = Dict{T,NT}()
    lengths = Dict{T,NTuple{3,NT}}()
    circumcenters = Dict{T,NTuple{2,NT}}()
    circumradii = Dict{T,NT}()
    angles = Dict{T,NTuple{3,NT}}()
    radius_edge_ratio = Dict{T,NT}()
    edge_midpoints = Dict{T,NTuple{3,NTuple{2,NT}}}()
    aspect_ratio = Dict{T,NT}()
    inradius = Dict{T,NT}()
    perimeter = Dict{T,NT}()
    centroid = Dict{T,NTuple{2,NT}}()
    offcenters = Dict{T,NTuple{2,NT}}()
    sinks = Dict{T,NTuple{2,NT}}()
    total_A = zero(NT)
    # The checks for each triangle use @check rather than a @test each, which was by far the slowest part of the refinement
    # tests (about 30 tests per triangle). The first triangle that fails each check is printed below.
    failures = Dict{String, Any}()
    for T in each_solid_triangle(tri)
        u, v, w = DT.triangle_vertices(T)
        p, q, r = get_point(tri, u, v, w)
        p = [getx(p), gety(p)]
        q = [getx(q), gety(q)]
        r = [getx(r), gety(r)]
        areas[triangle_vertices(T)] = 0.5 * (p[1] * (q[2] - r[2]) + q[1] * (r[2] - p[2]) + r[1] * (p[2] - q[2]))
        total_A += areas[triangle_vertices(T)]
        ℓ1 = norm(q - p)
        ℓ2 = norm(r - q)
        ℓ3 = norm(r - p)
        ℓmin = min(ℓ1, ℓ2, ℓ3)
        ℓmax = max(ℓ1, ℓ2, ℓ3)
        ℓmed = ℓ1 + ℓ2 + ℓ3 - ℓmin - ℓmax
        ℓ1, ℓ2, ℓ3 = ℓmin, ℓmed, ℓmax
        lengths[triangle_vertices(T)] = (ℓ1, ℓ2, ℓ3)
        r′ = p - r
        s′ = q - r
        ox = r[1] + det([norm(r′)^2 r′[2]; norm(s′)^2 s′[2]]) / (4areas[triangle_vertices(T)])
        oy = r[2] + det([r′[1] norm(r′)^2; s′[1] norm(s′)^2]) / (4areas[triangle_vertices(T)])
        circumcenters[triangle_vertices(T)] = (ox, oy)
        circumradii[triangle_vertices(T)] = norm(r - collect(circumcenters[triangle_vertices(T)]))
        all_angles = [(norm(p - r)^2 + norm(q - r)^2 - norm(p - q)^2) / (2norm(p - r) * norm(q - r)) for (p, q, r) in ((p, q, r), (q, r, p), (r, p, q))]
        all_angles[all_angles.<-1.0] .= -1.0
        all_angles[all_angles.>1.0] .= 1.0
        all_angles = acos.(all_angles)
        sort!(all_angles)
        radius_edge_ratio[triangle_vertices(T)] = circumradii[triangle_vertices(T)] / ℓ1
        edge_midpoints[triangle_vertices(T)] = ((Tuple(0.5 * (p + q))), Tuple(0.5 * (q + r)), Tuple(0.5 * (r + p)))
        inradius[triangle_vertices(T)] = 2areas[triangle_vertices(T)] / (ℓ1 + ℓ2 + ℓ3)
        perimeter[triangle_vertices(T)] = ℓ1 + ℓ2 + ℓ3
        aspect_ratio[triangle_vertices(T)] = inradius[triangle_vertices(T)] / circumradii[triangle_vertices(T)]
        centroid[triangle_vertices(T)] = (1 / 3 * (p[1] + q[1] + r[1]), 1 / 3 * (p[2] + q[2] + r[2]))
        angles[triangle_vertices(T)] = Tuple(all_angles)
        offcenters[triangle_vertices(T)] = DT.triangle_offcenter(p, q, r)
        sinks[triangle_vertices(T)] = DT.triangle_sink(tri, T)
        @check failures T radius_edge_ratio[triangle_vertices(T)] ≥ 1 / sqrt(3) - 0.1
        @check failures T DT.get_radius_edge_ratio(stats, T) ≥ 1 / sqrt(3) - 0.1
        @check failures T angles[triangle_vertices(T)][1] ≤ deg2rad(60) + 0.01
        @check failures T DT.get_minimum_angle(stats, T) ≤ deg2rad(60) + 0.01
    end

    ## Now compare the statistics 
    for T in each_solid_triangle(tri)
        @check failures T areas[triangle_vertices(T)] ≈ DT.get_area(stats, T) rtol = 1.0e-6 atol = 1.0e-8
        @check failures T areas[triangle_vertices(T)] ≈ DT.triangle_area((lengths[triangle_vertices(T)] .^ 2)...) rtol = 1.0e-6 atol = 1.0e-8
        @check failures T collect(lengths[triangle_vertices(T)]) ≈ collect(DT.get_lengths(stats, T)) rtol = 1.0e-6 atol = 1.0e-8
        @check failures T collect(circumcenters[triangle_vertices(T)]) ≈ collect(DT.get_circumcenter(stats, T)) rtol = 1.0e-6 atol = 1.0e-8
        @check failures T circumradii[triangle_vertices(T)] ≈ DT.get_circumradius(stats, T) rtol = 1.0e-6 atol = 1.0e-8
        @check failures T radius_edge_ratio[triangle_vertices(T)] ≈ DT.get_radius_edge_ratio(stats, T) rtol = 1.0e-6 atol = 1.0e-8
        @check failures T collect(collect.(edge_midpoints[triangle_vertices(T)])) ≈ collect(collect.(DT.get_edge_midpoints(stats, T))) rtol = 1.0e-6 atol = 1.0e-8
        @check failures T aspect_ratio[triangle_vertices(T)] ≈ DT.get_aspect_ratio(stats, T) rtol = 1.0e-6 atol = 1.0e-8
        @check failures T aspect_ratio[triangle_vertices(T)] ≈ DT.triangle_aspect_ratio(get_point(tri, triangle_vertices(T)...)...) rtol = 1.0e-6 atol = 1.0e-8
        @check failures T inradius[triangle_vertices(T)] ≈ DT.get_inradius(stats, T) rtol = 1.0e-6 atol = 1.0e-8
        @check failures T inradius[triangle_vertices(T)] ≈ DT.triangle_inradius(get_point(tri, triangle_vertices(T)...)...) rtol = 1.0e-6 atol = 1.0e-8
        @check failures T perimeter[triangle_vertices(T)] ≈ DT.get_perimeter(stats, T) rtol = 1.0e-6 atol = 1.0e-8
        @check failures T radius_edge_ratio[triangle_vertices(T)] ≈ 1 / (2sin(minimum(angles[triangle_vertices(T)]))) rtol = 1.0e-6 atol = 1.0e-8
        @check failures T (2sin(DT.get_minimum_angle(stats, T) / 2)^2 - 0.1 ≤ DT.get_aspect_ratio(stats, T) ≤ 2tan(DT.get_minimum_angle(stats, T) / 2) + 0.1)
        @check failures T DT.get_radius_edge_ratio(stats, T) ≈ 1 / (2(sin(DT.get_minimum_angle(stats, T)))) rtol = 1.0e-6 atol = 1.0e-8
        @check failures T areas[triangle_vertices(T)] ≈ inradius[triangle_vertices(T)] * 0.5perimeter[triangle_vertices(T)] rtol = 1.0e-6 atol = 1.0e-8
        @check failures T DT.get_area(stats, T) ≈ DT.get_inradius(stats, T) * 0.5DT.get_perimeter(stats, T) rtol = 1.0e-6 atol = 1.0e-8
        @check failures T collect(centroid[triangle_vertices(T)]) ≈ collect(DT.get_centroid(stats, T)) rtol = 1.0e-6 atol = 1.0e-8
        @check failures T DT.get_angles(stats, T)[1] ≈ angles[triangle_vertices(T)][1] rtol = 1.0e-6 atol = 1.0e-8
        @check failures T DT.get_angles(stats, T)[2] ≈ angles[triangle_vertices(T)][2] rtol = 1.0e-6 atol = 1.0e-8
        @check failures T DT.get_angles(stats, T)[3] ≈ angles[triangle_vertices(T)][3] rtol = 1.0e-6 atol = 1.0e-8
        @check failures T sum(DT.get_angles(stats, T)) ≈ π rtol = 1.0e-6 atol = 1.0e-8
        @check failures T DT.get_minimum_angle(stats, T) ≈ angles[triangle_vertices(T)][1] rtol = 1.0e-6 atol = 1.0e-8
        @check failures T DT.get_maximum_angle(stats, T) ≈ angles[triangle_vertices(T)][3] rtol = 1.0e-6 atol = 1.0e-8
        @check failures T DT.get_minimum_angle(stats, T) ≈ DT.get_angles(stats, T)[1] rtol = 1.0e-6 atol = 1.0e-8
        @check failures T DT.get_maximum_angle(stats, T) ≈ DT.get_angles(stats, T)[3] rtol = 1.0e-6 atol = 1.0e-8
        @check failures T collect(offcenters[triangle_vertices(T)]) ≈ collect(DT.get_offcenter(stats, T)) rtol = 1.0e-6 atol = 1.0e-8
        @check failures T collect(sinks[triangle_vertices(T)]) ≈ collect(DT.get_sink(stats, T)) rtol = 1.0e-6 atol = 1.0e-8
    end
    @test_with_log msg stats.individual_statistics == DT.get_individual_statistics(stats)
    @test_with_log msg stats.area ≈ DT.get_area(stats) rtol = 1.0e-6 atol = 1.0e-8
    @test_with_log msg stats.area ≈ total_A rtol = 1.0e-6 atol = 1.0e-8

    ## Inference 
    rng = StableRNG(1)
    for _ in 1:5
        _T = rand(rng, DT.each_solid_triangle(tri))
        _p, _q, _r = @inferred DT.get_point(tri, DT.triangle_vertices(_T)...)
        _ℓmin², _ℓmed², _ℓmax² = @inferred DT.squared_triangle_lengths(_p, _q, _r)
        _ℓmin, _ℓmed, _ℓmax = sqrt(_ℓmin²), sqrt(_ℓmed²), sqrt(_ℓmax²)
        _A = @inferred DT.triangle_area(_p, _q, _r)
        _circumcenter = @inferred DT.triangle_circumcenter(_p, _q, _r, _A)
        _circumradius = @inferred DT.triangle_circumradius(_A, _ℓmin², _ℓmed², _ℓmax²)
        _radius_edge_ratio = @inferred DT.triangle_radius_edge_ratio(_circumradius, _ℓmin)
        _edge_midpoints = @inferred DT.triangle_edge_midpoints(_p, _q, _r)
        _perimeter = @inferred DT.triangle_perimeter(_ℓmin, _ℓmed, _ℓmax)
        _inradius = @inferred DT.triangle_inradius(_A, _perimeter)
        _aspect_ratio = @inferred DT.triangle_aspect_ratio(_inradius, _circumradius)
        _centroid = @inferred DT.triangle_centroid(_p, _q, _r)
        _angles = @inferred DT.triangle_angles(_p, _q, _r)
        _offcenter = @inferred DT.triangle_offcenter(_p, _q, _r, _circumcenter)
        _sink = @inferred DT.triangle_sink(tri, _T)
        @inferred DT.triangle_orthocenter(tri, _T)
        @inferred DT.triangle_orthoradius_squared(tri, _T)
    end

    ## Test the number statistics 
    @test_with_log msg DT.num_vertices(stats) == length(all_vertices) == stats.num_vertices
    @test_with_log msg DT.num_solid_vertices(stats) == length(solid_vertices) == stats.num_solid_vertices
    @test_with_log msg DT.num_ghost_vertices(stats) == length(DT.all_ghost_vertices(tri)) == length(ghost_vertices) == stats.num_ghost_vertices
    @test_with_log msg DT.num_triangles(stats) == length(all_triangles) == stats.num_triangles
    @test_with_log msg DT.num_solid_triangles(stats) == length(solid_triangles) == stats.num_solid_triangles
    @test_with_log msg DT.num_ghost_triangles(stats) == length(ghost_triangles) == stats.num_ghost_triangles
    @test_with_log msg DT.num_edges(stats) == length(all_edges) == stats.num_edges
    @test_with_log msg DT.num_solid_edges(stats) == length(solid_edges) == stats.num_solid_edges
    @test_with_log msg DT.num_ghost_edges(stats) == length(ghost_edges) == stats.num_ghost_edges
    @test_with_log msg DT.num_boundary_segments(stats) == length(keys(get_boundary_edge_map(tri))) == stats.num_boundary_segments
    @test_with_log msg DT.num_interior_segments(stats) == num_edges(get_interior_segments(tri)) == stats.num_interior_segments
    @test_with_log msg DT.num_segments(stats) == num_edges(get_all_segments(tri)) == stats.num_segments
    @test_with_log msg DT.num_convex_hull_vertices(stats) == length(get_convex_hull_vertices(tri)) - 1 == stats.num_convex_hull_vertices

    ## Global statistics  
    smallest_angle = minimum([angles[triangle_vertices(T)][1] for T in each_solid_triangle(tri)])
    largest_angle = maximum([angles[triangle_vertices(T)][3] for T in each_solid_triangle(tri)])
    smallest_area = minimum([areas[triangle_vertices(T)] for T in each_solid_triangle(tri)])
    largest_area = maximum([areas[triangle_vertices(T)] for T in each_solid_triangle(tri)])
    smallest_radius_edge_ratio = minimum([radius_edge_ratio[triangle_vertices(T)] for T in each_solid_triangle(tri)])
    largest_radius_edge_ratio = maximum([radius_edge_ratio[triangle_vertices(T)] for T in each_solid_triangle(tri)])
    @test_with_log msg DT.get_smallest_angle(stats) ≈ smallest_angle rtol = 1.0e-6 atol = 1.0e-8
    @test_with_log msg DT.get_largest_angle(stats) ≈ largest_angle rtol = 1.0e-6 atol = 1.0e-8
    @test_with_log msg DT.get_smallest_area(stats) ≈ smallest_area rtol = 1.0e-6 atol = 1.0e-8
    @test_with_log msg DT.get_largest_area(stats) ≈ largest_area rtol = 1.0e-6 atol = 1.0e-8
    @test_with_log msg DT.get_smallest_radius_edge_ratio(stats) ≈ smallest_radius_edge_ratio rtol = 1.0e-6 atol = 1.0e-8
    @test_with_log msg DT.get_largest_radius_edge_ratio(stats) ≈ largest_radius_edge_ratio rtol = 1.0e-6 atol = 1.0e-8
    @test_with_log msg DT.get_smallest_radius_edge_ratio(stats) ≥ 1 / sqrt(3) - 0.1
    @test_with_log msg DT.get_smallest_angle(stats) ≤ deg2rad(60) + 0.01

    ## Orthocenter and orthoradius
    for T in each_solid_triangle(tri)
        o = DT.triangle_orthocenter(tri, T)
        i, j, k = DT.triangle_vertices(T)
        r1 = DT.dist_sqr(get_point(tri, i, o)...) - DT.get_weight(tri, i)
        r2 = DT.dist_sqr(get_point(tri, j, o)...) - DT.get_weight(tri, j)
        r3 = DT.dist_sqr(get_point(tri, k, o)...) - DT.get_weight(tri, k)
        @check failures T r1 ≈ r2 rtol = 1.0e-6 atol = 1.0e-8
        @check failures T r2 ≈ r3 rtol = 1.0e-6 atol = 1.0e-8
        @check failures T r1 ≈ r3 rtol = 1.0e-6 atol = 1.0e-8
        if !DT.is_weighted(tri)
            c = DT.triangle_circumcenter(tri, T)
            @check failures T collect(o) ≈ collect(c) rtol = 1.0e-6 atol = 1.0e-8
        end

        r² = DT.triangle_orthoradius_squared(tri, T)
        if !DT.is_weighted(tri)
            cr2 = DT.triangle_circumradius(get_point(tri, T...)...)^2
            @check failures T r² ≈ cr2 rtol = 1.0e-6 atol = 1.0e-8
        end
        @check failures T r² ≈ r1
        @check failures T r² ≈ r2
        @check failures T r² ≈ r3
        @check failures T DT.dist_sqr(get_point(tri, i, o)...) - DT.get_weight(tri, i) - r² ≈ 0 atol = 1.0e-8 * max(1, abs(r²))
    end
    @test_with_log msg isempty(failures)
    isempty(failures) || println("  Failed checks in validate_statistics, with the first triangle that failed each: ", failures)
end

## TODO: Implement a brute-force DT.VoronoiTessellation that we can compare with
function validate_tessellation(vorn::DT.VoronoiTessellation; check_convex=true, check_area=true, check_adjacent=true, predicates::DT.AbstractPredicateKernel=DT.AdaptiveKernel())
    tri = DT.get_triangulation(vorn)
    for (i, p) in DT.get_generators(vorn)
        flag = get_point(tri, i) == get_generator(vorn, i) == p
        if !flag
            println("Generator $i is not correct, mapped to $p.")
            return false
        end
    end
    flag = DT.get_triangulation(vorn) == vorn.triangulation
    if !flag
        println("DT.Triangulation is not correct.")
        return false
    end
    circumcenter_to_triangle = DT.get_circumcenter_to_triangle(vorn)
    triangle_to_circumcenter = DT.get_triangle_to_circumcenter(vorn)
    for V in DT.each_solid_triangle(DT.get_triangulation(vorn))
        V = DT.sort_triangle(V)
        c = DT.get_triangle_to_circumcenter(vorn, V)
        c = get_polygon_point(vorn, c)
        i, j, k = triangle_vertices(V)
        p, q, r = get_point(DT.get_triangulation(vorn), i, j, k)
        if DT.is_weighted(tri)
            cx, cy = DT.triangle_orthocenter(tri, V)
        else
            cx, cy = DT.triangle_circumcenter(p, q, r)
        end
        flag = cx == c[1] && cy == c[2]
        if !flag
            println("Circumcenter of $V is not correct, mapped to $c.")
            return false
        end
    end
    for (c, V) in circumcenter_to_triangle
        flag = DT.get_circumcenter_to_triangle(vorn, c) == V
        if !flag
            println("Circumcenter $c is not correct, mapped to $V.")
            return false
        end
        flag = DT.get_triangle_to_circumcenter(vorn, V) == c
        if !flag
            println("Triangle $V is not correct, mapped to $c.")
            return false
        end
    end
    for (V, c) in triangle_to_circumcenter
        flag = DT.get_circumcenter_to_triangle(vorn, c) == V
        if !flag # cocircular points disrupt these dicts somewhat
            _flag = c ∈ vorn.cocircular_circumcenters
            if !_flag
                println("Circumcenter $c is duplicated from another triangle, but not in cocircular_circumcenters.")
                return false
            end
            i, j, k = triangle_vertices(V)
            for i in (i, j, k)
                for e in each_edge(get_adjacent2vertex(DT.get_triangulation(vorn), i))
                    v, w = DT.edge_vertices(e)
                    V′ = DT.sort_triangle(DT.construct_triangle(DT.triangle_type(DT.get_triangulation(vorn)), i, v, w))
                    flag = flag || DT.get_triangle_to_circumcenter(vorn, V′) == c
                end
            end
        end
        if !flag
            println("Triangle $V is not correct, mapped to $c.")
            return false
        end
        flag = DT.get_triangle_to_circumcenter(vorn, V) == c
        if !flag
            println("Circumcenter $c is not correct, mapped to $V.")
            return false
        end
    end
    for i in each_polygon_index(vorn)
        A = get_area(vorn, i)
        flag = A ≥ 0
        if !flag
            println("Polygon $i has area $A.")
            return false
        end
        if !isfinite(A)
            flag = i ∈ DT.get_unbounded_polygons(vorn)
            if !flag
                println("Polygon $i has infinite area, but is not in the unbounded polygons.")
                return false
            end
        end
    end
    if check_adjacent
        for i in each_polygon_index(vorn)
            C = get_polygon(vorn, i)
            ne = DT.num_boundary_edges(C)
            for j in 1:ne
                u = get_boundary_nodes(C, j)
                v = get_boundary_nodes(C, j + 1)
                flag = get_adjacent(vorn, u, v) == get_adjacent(vorn, DT.construct_edge(DT.edge_type(DT.get_triangulation(vorn)), u, v)) == i
                if !flag
                    println("Polygon $i is not adjacent to points $u and $v.")
                    return false
                end
            end
        end
    end
    for i in each_polygon_index(vorn)
        if i ∉ DT.get_unbounded_polygons(vorn)
            verts = get_polygon(vorn, i)
            poly_points = get_polygon_point.(Ref(vorn), verts)
            flag = @views allunique(poly_points[begin:(end-1)])
            if !flag
                println("Polygon $i has repeated vertices.")
                return false
            end
            if check_convex
                _pts = unique(poly_points)
                ch = convex_hull(_pts; predicates)
                _poly_points = _pts[DT.get_vertices(ch)]
                flag = DT.circular_equality(collect.(poly_points), collect.(_poly_points), ≈)
                if !flag
                    println("Polygon $i is not convex.")
                    return false
                end
            end
        end
    end
    if isempty(DT.get_unbounded_polygons(vorn)) && check_area
        A = 0.0
        for i in each_polygon_index(vorn)
            A += get_area(vorn, i)
        end
        @test isapprox(A, get_area(vorn.triangulation), rtol=1.0e-4)
    end
    return true
end

function manual_insertion_event_history(tri::Triangulation, orig_tri::Triangulation)
    # Initialise
    added_triangles = Set{DT.triangle_type(tri)}()
    deleted_triangles = Set{DT.triangle_type(tri)}()
    added_segments = Set{DT.edge_type(tri)}()
    deleted_segments = Set{DT.edge_type(tri)}()
    added_boundary_segments = Set{DT.edge_type(tri)}()
    deleted_boundary_segments = Set{DT.edge_type(tri)}()
    # Find added or deleted triangles
    for T in each_triangle(tri)
        _, flag = DT.contains_triangle(orig_tri, T)
        !flag && push!(added_triangles, T)
    end
    for T in each_triangle(orig_tri)
        _, flag = DT.contains_triangle(tri, T)
        !flag && push!(deleted_triangles, T)
    end
    # Find added or deleted segments 
    for e in DT.get_interior_segments(tri)
        flag = DT.contains_unoriented_edge(e, DT.get_interior_segments(orig_tri))
        !flag && push!(added_segments, e)
    end
    for e in DT.get_interior_segments(orig_tri)
        flag = DT.contains_unoriented_edge(e, DT.get_interior_segments(tri))
        !flag && push!(deleted_segments, e)
    end
    # Find added or delete boundary segments 
    for e in DT.each_boundary_edge(tri)
        flag = DT.contains_unoriented_edge(e, DT.each_boundary_edge(orig_tri))
        !flag && push!(added_boundary_segments, e)
    end
    for e in DT.each_boundary_edge(orig_tri)
        flag = DT.contains_unoriented_edge(e, DT.each_boundary_edge(tri))
        !flag && push!(deleted_boundary_segments, e)
    end
    # Done 
    return DT.InsertionEventHistory(added_triangles, deleted_triangles, added_segments, deleted_segments, added_boundary_segments, deleted_boundary_segments)
end

function validate_insertion_event_history(tri::Triangulation, orig_tri::Triangulation, history::DT.InsertionEventHistory)
    manual_history = manual_insertion_event_history(tri, orig_tri)
    @test length(manual_history.added_triangles) == length(history.added_triangles)
    @test length(manual_history.deleted_triangles) == length(history.deleted_triangles)
    @test length(manual_history.added_segments) == length(history.added_segments)
    @test length(manual_history.deleted_segments) == length(history.deleted_segments)
    @test length(manual_history.added_boundary_segments) == length(history.added_boundary_segments)
    @test length(manual_history.deleted_boundary_segments) == length(history.deleted_boundary_segments)
    @test DT.compare_triangle_collections(manual_history.added_triangles, history.added_triangles)
    @test DT.compare_unoriented_edge_collections(manual_history.added_segments, history.added_segments)
    @test DT.compare_unoriented_edge_collections(manual_history.added_boundary_segments, history.added_boundary_segments)
end

function compare_encroach_queues(args::DT.RefinementArguments, manual_enqueue)
    _manual_enqueue_pairs = collect(manual_enqueue)
    _args_queue_segments_pairs = collect(args.queue.segments)
    for pars in (_manual_enqueue_pairs, _args_queue_segments_pairs)
        for i in eachindex(pars)
            u, v = pars[i].first
            if u > v
                pars[i] = (v, u) => pars[i].second
            end
        end
    end
    _compare_pairs(_manual_enqueue_pairs, _args_queue_segments_pairs)
end

function compare_triangle_queues(args::DT.RefinementArguments, manual_enqueue)
    _manual_enqueue_pairs = collect(manual_enqueue)
    _args_queue_triangle_pairs = collect(args.queue.triangles)
    for pars in (_manual_enqueue_pairs, _args_queue_triangle_pairs)
        for i in eachindex(pars)
            T = pars[i].first
            pars[i] = DT.sort_triangle(T) => pars[i].second
        end
    end
    _compare_pairs(_manual_enqueue_pairs, _args_queue_triangle_pairs)
end

function is_conformal(tri::Triangulation; predicates::DT.AbstractPredicateKernel=DT.AdaptiveKernel())
    points = get_points(tri)
    segment_tree = DT.BoundaryRTree(points)
    for e in each_segment(tri)
        i, j = DT.edge_vertices(e)
        insert!(segment_tree, i, j)
    end
    for r in each_solid_vertex(tri)
        intersects = DT.get_intersections(segment_tree, r, cache_id=1) # can't use multithreading here
        c = get_point(tri, r)
        for box in intersects
            i, j = DT.get_edge(box)
            any(==(r), (i, j)) && continue
            a, b = get_point(tri, i, j)
            flag = DT.is_inside(DT.point_position_relative_to_diametral_circle(predicates, a, b, c))
            if flag
                #cert = test_visibility(tri, segment_tree, i, j, r)
                cert = DT.test_visibility(predicates, tri, i, j, r)
                flag = DT.is_visible(cert)
            end
            if flag
                println("Triangulation is not conformal, as vertex $r is in the diametral circle of the edge $((i, j))")
                return false
            end
        end
    end
    return true
end

function slow_triangle_assess(tri, args)
    good_T = NTuple{3,Int}[]
    bad_T = NTuple{3,Int}[]
    for T in each_solid_triangle(tri)
        u, v, w = T
        p, q, r = get_point(tri, u, v, w)
        t1, t2, t3 = DT.triangle_angles(p, q, r)
        t = min(t1, t2, t3) * 180 / π
        t′ = max(t1, t2, t3) * 180 / π
        A = DT.triangle_area(p, q, r)
        seditious = DT.is_triangle_seditious(tri, args, u, v, w, DT.squared_triangle_lengths_and_smallest_index(p, q, r)[4])
        nestled = DT.is_triangle_nestled(tri, T, DT.squared_triangle_lengths_and_smallest_index(p, q, r)[4])
        ρ, qual = DT.assess_triangle_quality(tri, args, T)
        @inferred DT.assess_triangle_quality(tri, args, T)
        condition = ((A > args.constraints.max_area) || (t < args.constraints.min_angle && ρ > args.constraints.max_radius_edge_ratio && !seditious && !nestled)) && A ≥ args.constraints.min_area && !args.constraints.custom_constraint(tri, T) && t′ ≤ args.constraints.max_angle
        @test ρ == DT.triangle_radius_edge_ratio(p, q, r)
        @test qual == condition
        if qual
            push!(bad_T, T)
        else
            push!(good_T, T)
        end
    end
    return good_T, bad_T
end

function slow_triangle_assess_queue(tri, args)
    good_T, bad_T = slow_triangle_assess(tri, args)
    queue = PriorityQueue{NTuple{3,Int},Float64}(Base.Order.Reverse)
    for T in bad_T
        queue[T] = DT.triangle_radius_edge_ratio(get_point(tri, T...)...)
    end
    return queue
end

function compute_φmin(tri)
    # φmin is the minimum angle between two adjoining segments 
    φmin = Inf
    locked = !DT.has_boundary_nodes(tri)
    if !DT.has_boundary_nodes(tri)
        lock_convex_hull!(tri; rng=StableRNG(1))
    end
    for e in each_segment(tri)
        for η in each_segment(tri)
            if e ≠ η
                u, v = e
                x, y = η
                w = DT.get_shared_vertex(e, η)
                w == DT.∅ && continue
                if w == u && w == x
                    e1 = (u, v)
                    e2 = (x, y)
                    i, j, k = w, v, y
                elseif w == u && w == y
                    e1 = (u, v)
                    e2 = (y, x)
                    i, j, k = w, v, x
                elseif w == v && w == x
                    e1 = (v, u)
                    e2 = (x, y)
                    i, j, k = w, u, y
                else # if w == v && w == y
                    e1 = (v, u)
                    e2 = (y, x)
                    i, j, k = w, u, x
                end
                p, q, r = get_point(tri, i, j, k)
                s1 = q .- p
                s2 = r .- p
                rat = dot(s1, s2) / (norm(s1) * norm(s2))
                if rat > 1
                    rat = 1
                elseif rat < -1
                    rat = -1
                end
                θ = acosd(rat)
                φmin = min(φmin, θ)
            end
        end
    end
    locked && unlock_convex_hull!(tri)
    return φmin
end

function validate_refinement(tri, args; check_conformal=true, warn=true)
    rng = StableRNG(1) # for find_triangle
    ## Things to check:
    ##   1. All angle constraints are met, except for seditious and nestled triangles.
    ##   2. If !use_lens, the triangulation is conformal. 
    ##   3. All area constraints are met.
    ##   4. All custom constraints are met.
    ##   5. The triangulation is Delaunay.
    ##   6. All circumcenters are inside the domain, if !use_lens. 
    ##   7. There are no more points than there are maximum points.
    ##   8. There are no more encroached edges.
    ##   9. There are no more bad triangles.
    minθ = Inf
    maxθ = -Inf
    all_segments = Int[]
    for e in each_segment(tri)
        push!(all_segments, e...)
    end
    for T in each_solid_triangle(tri)
        is_in_corner = all(∈(all_segments), triangle_vertices(T)) # triangles in corners sometimes can't have their circumcenters inserted if doing use_lens, so their centroid gets inserted, but this breaks the constraint conditions sometimes for minimum angles.
        # Get some stats
        t1, t2, t3 = DT.triangle_angles(get_point(tri, T...)...)
        t1, t2, t3 = (t1, t2, t3) .* 180 ./ π
        A = DT.triangle_area(get_point(tri, T...)...)
        seditious = DT.is_triangle_seditious(tri, args, T..., DT.squared_triangle_lengths_and_smallest_index(get_point(tri, T...)...)[4])
        nestled = DT.is_triangle_nestled(tri, T, DT.squared_triangle_lengths_and_smallest_index(get_point(tri, T...)...)[4])
        _flag, steiner_point = DT.get_steiner_point(tri, args, T)
        show_print = true
        if !DT.is_none(_flag)
            warn && @warn "Steiner point found for triangle $T has precision issues. Proceeding with the tests, but will not return any results for this triangle."
            show_print = warn
        else
            # Only update extrema for good triangles 
            minθ = min(minθ, t1, t2, t3)
            maxθ = max(maxθ, t1, t2, t3)
        end
        V = find_triangle(tri, steiner_point; predicates=DT.ExactKernel(), rng)
        flag = DT.point_position_relative_to_triangle(DT.ExactKernel(), tri, V, steiner_point)
        if DT.is_on(flag) && DT.is_ghost_triangle(V)
            V = DT.replace_ghost_triangle_with_boundary_triangle(tri, V)
        end
        if !(seditious || nestled || is_in_corner) && min(t1, t2, t3) < args.constraints.min_angle
            show_print && println("Triangle $T has minimum angle $(min(t1, t2, t3)), which is less than the minimum angle constraint $(args.constraints.min_angle).")
            DT.is_none(_flag) && return false
        end
        if max(t1, t2, t3) > args.constraints.max_angle
            show_print && println("Triangle $T has maximum angle $(max(t1, t2, t3)), which is greater than the maximum angle constraint $(args.constraints.max_angle).")
            return false
        end
        flag = args.constraints.min_area - 1.0e-16 ≤ A ≤ args.constraints.max_area
        if !flag
            show_print && println("Triangle $T has area $A, which is not in the range $(args.constraints.min_area) to $(args.constraints.max_area).")
            DT.is_none(_flag) && return false
        end
        flag = !args.constraints.custom_constraint(tri, T)
        if !flag
            show_print && println("Triangle $T violates the custom constraint.")
            DT.is_none(_flag) && return false
        end
        if !args.use_lens
            flag = DT.dist(tri, steiner_point) > -eps(Float64)
            if !flag
                show_print && println("The Steiner point associated with the triangle $T is a distance $(-DT.dist(tri, steiner_point)) away from the domain.")
                DT.is_none(_flag) && return false
            end
        end
        for (u, v) in DT.triangle_edges(T)
            if DT.contains_segment(tri, u, v)
                flag = DT.is_encroached(tri, args, (u, v))
                if flag
                    show_print && println("The edge ($u, $v) is encroached.")
                    DT.is_none(_flag) && return false
                end
            end
        end
        flag = DT.assess_triangle_quality(tri, args, T)[2]
        if flag && !seditious && !nestled && !is_in_corner
            show_print && println("Triangle $T is of poor quality.")
            DT.is_none(_flag) && return false
        end
        ρ = DT.triangle_radius_edge_ratio(get_point(tri, T...)...)
        flag_2 = flag == ((A > args.constraints.max_area) || (min(t1, t2, t3) < args.constraints.min_angle && ρ > args.constraints.max_radius_edge_ratio && !seditious && !nestled && !is_in_corner)) && A ≥ args.constraints.min_area && !args.constraints.custom_constraint(tri, T) && max(t1, t2, t3) ≤ args.constraints.max_angle
        if !flag_2 && !is_in_corner
            if flag
                show_print && println("Triangle $T was incorrectly marked as being of good quality.")
            else
                show_print && println("Triangle $T was incorrectly marked as being of poor quality.")
            end
            DT.is_none(_flag) && return false
        end
    end
    if !args.use_lens && check_conformal
        flag = is_conformal(tri; predicates=DT.ExactKernel())
        if !flag
            println("The triangulation is not conformal.")
            return false
        end
    end
    flag = validate_triangulation(tri)
    if !flag
        println("The triangulation is not Delaunay.")
        return false
    end
    flag = DT.num_solid_vertices(tri) > args.constraints.max_points
    if flag
        println("The number of points in the triangulation ($(DT.num_points(tri))) is greater than the maximum number of points allowed ($(args.constraints.max_points)).")
        return false
    end
    if !args.use_lens
        φmin = compute_φmin(tri)
        ρ = args.constraints.max_radius_edge_ratio
        τmin = asind(sind(φmin) / sqrt(5 - 4cosd(φmin)))
        τmax = 180 - 2asind(1 / (2ρ))
        flag = minθ ≥ min(τmin, args.constraints.min_angle)
        if !flag
            println("The minimum angle $(minθ) is not in the range $(τmin) to $(τmax).")
            return false
        end
        flag = maxθ ≤ τmax
        if !flag
            println("The maximum angle $(maxθ) is not in the range $(τmin) to $(τmax).")
            return false
        end
    end
    return true
end

validate_refinement(tri; check_conformal=true, warn=true, predicates=DT.ExactKernel(), kwargs...) = validate_refinement(tri, DT.RefinementArguments(tri; predicates, kwargs...); warn, check_conformal)
