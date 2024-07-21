using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using DataStructures
using StaticArrays

global def_adj = DT.∅
global default_1 = Dict{NTuple{2, Int}, Int}()
global default_2 = Dict{NTuple{2, Int32}, Int32}()
global default_3 = Dict{Vector{Int}, Int}()
global adj_1 = DT.Adjacent{Int, NTuple{2, Int}}()
global adj_2 = DT.Adjacent{Int32, NTuple{2, Int32}}()
global adj_3 = DT.Adjacent{Int, Vector{Int}}()

@testset "Constructors and getters" begin
    @test adj_1.adjacent == default_1
    @test adj_2.adjacent == default_2
    @test adj_3.adjacent == default_3
    @test adj_1 == DT.Adjacent(default_1)
    @test adj_2 == DT.Adjacent(default_2)
    @test adj_3 == DT.Adjacent(default_3)
    @test DT.get_adjacent(adj_1) == default_1
    @test DT.get_adjacent(adj_2) == default_2
    @test DT.get_adjacent(adj_3) == default_3
end

global dict_1 = Dict((1, 2) => 4, (2, 3) => 10, (5, 6) => 15, (20, 5) => 72)
global dict_2 = Dict(
    @SVector[1, 2] => 4, @SVector[2, 3] => 10, @SVector[5, 6] => 15,
    @SVector[20, 5] => 72,
)
global dict_3 = Dict{NTuple{2, Int32}, Int32}(
    (1, 2) => 4, (2, 3) => 10, (5, 6) => 15,
    (20, 5) => 72,
)
global ddict_1 = Dict(dict_1)
global ddict_2 = Dict(dict_2)
global ddict_3 = Dict(dict_3)
global adj_1 = DT.Adjacent(ddict_1)
global adj_2 = DT.Adjacent(ddict_2)
global adj_3 = DT.Adjacent(ddict_3)

@testset "Getting, adding, and deleting adjacencies" begin
    for adj in (adj_1, adj_2, adj_3)
        @test get_adjacent(adj, 1, 2) == 4
        @inferred get_adjacent(adj, 1, 2)
        @test get_adjacent(adj, 2, 3) == 10
        @test get_adjacent(adj, 5, 6) == 15
        @test get_adjacent(adj, 20, 5) == 72
        @test get_adjacent(adj, 20, 27) == def_adj

        DT.add_adjacent!(adj, 17, 23, 50)
        DT.add_adjacent!(adj, 28, 38, 173)
        @test get_adjacent(adj, 17, 23) == 50
        @test get_adjacent(adj, 28, 38) == 173

        DT.delete_adjacent!(adj, 2, 3)
        DT.delete_adjacent!(adj, 20, 5)
        @test get_adjacent(adj, 2, 3) == def_adj
        @test get_adjacent(adj, 20, 5) == def_adj

        DT.add_triangle!(adj, 51, 52, 53)
        @test get_adjacent(adj, 51, 52) == 53
        @test get_adjacent(adj, 52, 53) == 51
        @test get_adjacent(adj, 53, 51) == 52

        DT.add_triangle!(adj, (60, 61, 62))
        @test get_adjacent(adj, 60, 61) == 62
        @test get_adjacent(adj, 61, 62) == 60
        @test get_adjacent(adj, 62, 60) == 61

        T1 = (120, 125, 132)
        T2 = (501, 502, 591)
        T3 = (6019, 919, 821)
        DT.add_triangle!.(Ref(adj), (T1, T2, T3))
        for T in (T1, T2, T3)
            i, j, k = T
            @test get_adjacent(adj, i, j) == k
            @test get_adjacent(adj, j, k) == i
            @test get_adjacent(adj, k, i) == j
        end

        DT.delete_triangle!(adj, 60, 61, 62)
        @test get_adjacent(adj, 60, 61) == def_adj
        @inferred get_adjacent(adj, 60, 61)
        @test get_adjacent(adj, 61, 62) == def_adj
        @test get_adjacent(adj, 62, 60) == def_adj

        DT.delete_triangle!(adj, T2)
        @test get_adjacent(adj, 501, 502) == def_adj
        @test get_adjacent(adj, 502, 591) == def_adj
        @test get_adjacent(adj, 591, 501) == def_adj

        DT.delete_triangle!.(Ref(adj), (T1, T3))
        for T in (T1, T2)
            i, j, k = T
            @test get_adjacent(adj, i, j) == def_adj
            @test get_adjacent(adj, j, k) == def_adj
            @test get_adjacent(adj, k, i) == def_adj
        end
    end
end

@testset "Testing if edges exist" begin
    tri = triangulate(rand(2, 50))
    @test !DT.edge_exists(get_adjacent(tri, 122, -7))
end

@testset "empty!" begin
    tri = triangulate(rand(2, 50))
    adj = get_adjacent(tri)
    @test !isempty(get_adjacent(adj))
    empty!(adj)
    @test isempty(get_adjacent(adj))
end
