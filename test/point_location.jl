#=
@testset "Can we make sure we do not sample a ghost edge?" begin
    p1 = @SVector[-3.32, 3.53]
    p2 = @SVector[-5.98, 2.17]
    p3 = @SVector[-6.36, -1.55]
    p4 = @SVector[-2.26, -4.31]
    p5 = @SVector[6.34, -3.23]
    p6 = @SVector[-3.24, 1.01]
    p7 = @SVector[0.14, -1.51]
    p8 = @SVector[0.2, 1.25]
    p9 = @SVector[1.0, 4.0]
    p10 = @SVector[4.74, 2.21]
    p11 = @SVector[2.32, -0.27]
    pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11]
    T, adj, adj2v, DG, HG = DT.triangulate_berg(pts)
    [DT.delete_triangle!(i, j, k, T, adj, adj2v, DG) for (i, j, k) in (
        (1, 8, 9), (9, 8, 10), (10, 11, 5), (5, 7, 4), (4, 6, 3), (3, 6, 2), (2, 6, 1)
    )]
    [DT.add_triangle!(i, j, k, T, adj, adj2v, DG; update_ghost_edges=true) for (i, j, k) in (
        (1, 8, 9), (9, 8, 10), (10, 11, 5), (5, 7, 4), (4, 6, 3), (3, 6, 2), (2, 6, 1)
    )]
    DT.clear_empty_keys!(adj, DG)
    for k in [9, 10, 5, 4, 3, 2, 1]
        for _ in 1:508
            r = rand(1:20022991)
            Random.seed!(r)
            i, j = rand(DT.get_edge(adj2v, k))
            Random.seed!(r)
            _i, _j = DT.select_random_edge(k, adj2v; include_ghost_edges=true)
            @test (i, j) == (_i, _j)
        end
        for _ in 1:50409
            i, j = DT.select_random_edge(k, adj2v; include_ghost_edges=false)
            @test !DT.is_ghost_edge(i, j)
        end
    end
end
=#

@testset "Point location with DAG or jump-and-march" begin
    for n = 3:111
        r = 5sqrt.(rand(n))
        θ = 2π * rand(n)
        pts = [@SVector[r * cos(θ), r * sin(θ)] for (r, θ) in zip(r, θ)]
        T, adj, adj2v, DG, HG = DT.triangulate_berg(pts)
        _pts = deepcopy(pts)
        for k in rand(1:n, ceil(Int64, 0.45n^(1 / 3)))
            for _ in 1:length(T)
                sample_T = rand(T)
                a, b, c = DT.get_point(pts, sample_T...)
                q = (a .+ b .+ c) ./ 3
                push!(_pts, q)
                τ, flag = locate_triangle(HG, _pts, n + 1, DT.BoundingTriangle)
                @test DT.isintriangle(τ, _pts, n + 1) == flag == 1
                @test DT.isoriented(τ, pts) == 1
                for _ in 1:4
                    τ = jump_and_march(q, adj, adj2v, DG, pts; k=k)
                    @test DT.isintriangle(τ, _pts, n + 1) == 1
                    @test DT.isoriented(τ, pts) == 1
                    τ = jump_and_march(q, adj, adj2v, DG, pts)
                    @test DT.isintriangle(τ, _pts, n + 1) == 1
                    @test DT.isoriented(τ, pts) == 1
                    τ = jump_and_march(n + 1, adj, adj2v, DG, _pts; k=k, pt_idx=1:n)
                    @test DT.isintriangle(τ, _pts, n + 1) == 1
                    @test DT.isoriented(τ, pts) == 1
                    τ = jump_and_march(n + 1, adj, adj2v, DG, _pts; pt_idx=1:n)
                    @test DT.isintriangle(τ, _pts, n + 1) == 1
                    @test DT.isoriented(τ, pts) == 1

                    τ = jump_and_march(q, adj, adj2v, DG, pts; k=k, pt_idx=[1])
                    @test DT.isintriangle(τ, _pts, n + 1) == 1
                    @test DT.isoriented(τ, pts) == 1
                    τ = jump_and_march(q, adj, adj2v, DG, pts; pt_idx=[1])
                    @test DT.isintriangle(τ, _pts, n + 1) == 1
                    @test DT.isoriented(τ, pts) == 1
                    τ = jump_and_march(n + 1, adj, adj2v, DG, _pts; k=k, pt_idx=[1])
                    @test DT.isintriangle(τ, _pts, n + 1) == 1
                    @test DT.isoriented(τ, pts) == 1
                    τ = jump_and_march(n + 1, adj, adj2v, DG, _pts; pt_idx=[1])
                    @test DT.isintriangle(τ, _pts, n + 1) == 1
                    @test DT.isoriented(τ, pts) == 1
                end
                pop!(_pts)
                for q in (a, b, c)
                    push!(_pts, q)
                    τ, flag = locate_triangle(HG, _pts, n + 1, DT.BoundingTriangle)
                    @test DT.isintriangle(τ, _pts, n + 1) == flag == 0
                    @test DT.isoriented(τ, pts) == 1
                    for _ in 1:4
                        τ = jump_and_march(q, adj, adj2v, DG, pts; k=k)
                        @test DT.isintriangle(τ, _pts, n + 1) == 0
                        @test DT.isoriented(τ, pts) == 1
                        τ = jump_and_march(q, adj, adj2v, DG, pts)
                        @test DT.isintriangle(τ, _pts, n + 1) == 0
                        @test DT.isoriented(τ, pts) == 1
                    end
                    pop!(_pts)
                end
            end
        end
    end
end
