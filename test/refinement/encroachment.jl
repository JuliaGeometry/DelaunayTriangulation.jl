using ..DelaunayTriangulation
const DT = DelaunayTriangulation

include("../helper_functions.jl")


@testset "Testing if point is in diametral circle" begin
    p = Vector{NTuple{2,Float64}}(undef, 7)
    q = Vector{NTuple{2,Float64}}(undef, 7)
    r = Vector{NTuple{2,Float64}}(undef, 7)

    p[1], q[1], r[1] = (2.6, 1.63), (5.68, 1.37), (4.5, 4.63)
    p[2], q[2], r[2] = p[1], q[1], (4.18, 1.96)
    p[3], q[3], r[3] = p[2], q[2], (1.5, 2.0)
    p[4], q[4], r[4] = p[3], q[3], (3.9, 1.324)
    p[5], q[5], r[5] = p[4], q[4], (3.6, -1.25)
    p[6], q[6], r[6] = (2.0, 0.0), (4.0, 0.0), (4.0, 2.0)
    p[7], q[7], r[7] = (2.0, 0.0), (4.0, 2.0), (4.0, 0.0)

    certs = [
        DT.is_outside,
        DT.is_inside,
        DT.is_outside,
        DT.is_inside,
        DT.is_outside,
        DT.is_outside,
        DT.is_on,
    ]

    for (p, q, r, c) in zip(p, q, r, certs)
        @test c(DT.point_position_relative_to_diametral_circle(p, q, r))
        @test c(DT.point_position_relative_to_diametral_circle(q, p, r))
        if c == DT.is_on || c == DT.is_inside
            @test DT.encroaches_upon(p, q, r)
        else
            @test !DT.encroaches_upon(p, q, r)
        end
    end
end

if !(get(ENV, "CI", "false") == "true")
    @testset "Testing if encroached edges are detected" begin
        _x, _y = complicated_geometry()
        x = _x
        y = _y
        tri_1 = generate_mesh(x, y, 0.1; convert_result=true, add_ghost_triangles=true)
        tri_2 = generate_mesh(x[1], y[1], 0.1; convert_result=true, add_ghost_triangles=true)
        tri_3 = generate_mesh([0.0, 2.0, 2.0, 0.0, 0.0], [0.0, 0.0, 2.0, 2.0, 0.0], 0.1;
            convert_result=true, add_ghost_triangles=true)
        tri_4 = generate_mesh(reverse(reverse.(x[2])), reverse(reverse.(y[2])), 0.1; convert_result=true, add_ghost_triangles=true)
        a, b = 0.0, 5.0
        c, d = 3.0, 7.0
        nx = 3
        ny = 3
        tri_5 = triangulate_rectangle(a, b, c, d, nx, ny; add_ghost_triangles=true, single_boundary=false)
        tri_6 = triangulate_rectangle(a, b, c, d, nx, ny; add_ghost_triangles=true, single_boundary=true)
        tri_7 = triangulate(rand(2, 500))
        tri_8 = triangulate(rand(2, 500), delete_ghosts=false)
        for tri in (tri_1, tri_2, tri_3, tri_4, tri_5, tri_6, tri_7, tri_8)
            in_dt_encroached_edges, not_in_dt_encroached_edges = slow_encroachment_test(tri)
            all_bn = DT.get_all_boundary_nodes(tri)
            for (e, (b, k)) in not_in_dt_encroached_edges
                if DT.initial(e) ∈ all_bn && DT.terminal(e) ∈ all_bn && !DT.contains_constrained_edge(tri, e) # e.g. if an edge crosses an interior
                    continue
                end
                @test DT.is_encroached(tri, e) == b
            end
            for (e, (b, k)) in in_dt_encroached_edges
                if DT.initial(e) ∈ all_bn && DT.terminal(e) ∈ all_bn && !DT.contains_constrained_edge(tri, e)
                    continue
                end
                @test DT.is_encroached(tri, e) == b
            end
        end
    end
end

@testset "Concentric shell splitting" begin
    for _ in 1:5000
        p = rand(2)
        q = rand(2)
        t = DT.compute_concentric_shell_ternary_split_position(p, q)
        @test 1 / 3 ≤ t ≤ 2 / 3 # Split should between 1/3 and 2/3 of the segment pq 
        t = DT.compute_concentric_shell_quarternary_split_position(p, q)
        @test 1 / 4 ≤ 1 / 2 # Split should between 1/4 and 3/4 of the segment pq
    end
end
