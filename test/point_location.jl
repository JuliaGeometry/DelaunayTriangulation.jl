@testset "Point location with DAG or jump-and-march" begin
    for n = 3:350
        r = 5sqrt.(rand(n))
        θ = 2π * rand(n)
        pts = [@SVector[r * cos(θ), r * sin(θ)] for (r, θ) in zip(r, θ)]
        T, adj, adj2v, DG, HG = DT.triangulate_berg(pts)
        _pts = deepcopy(pts)
        for k in rand(1:n, ceil(Int64, 0.45n^(1 / 3)))
            for _ in 1:length(T)
                sample_T = rand(T)
                a, b, c = DT.get_point(pts, sample_T...)
                q = (a + b + c) / 3
                push!(_pts, q)
                τ, flag = locate_triangle(HG, _pts, n + 1, DT.BoundingTriangle)
                @test DT.isintriangle(τ, _pts, n + 1) == flag == 1
                @test DT.isoriented(τ, pts) == 1
                for _ in 1:4
                    τ = jump_and_march(q, adj, adj2v, pts; k=k)
                    @test DT.isintriangle(τ, _pts, n + 1) == 1
                    @test DT.isoriented(τ, pts) == 1
                    τ = jump_and_march(q, adj, adj2v, pts)
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
                        τ = jump_and_march(q, adj, adj2v, pts; k=k)
                        @test DT.isintriangle(τ, _pts, n + 1) == 0
                        @test DT.isoriented(τ, pts) == 1
                        τ = jump_and_march(q, adj, adj2v, pts)
                        @test DT.isintriangle(τ, _pts, n + 1) == 0
                        @test DT.isoriented(τ, pts) == 1
                    end
                    pop!(_pts)
                end
            end
        end
    end
end
