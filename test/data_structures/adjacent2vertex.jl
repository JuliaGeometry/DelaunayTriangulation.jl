using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using DataStructures
using StaticArrays

global dict_1 = Dict{Int,Set{NTuple{2,Int}}}()
global dict_2 = Dict{Int,Vector{NTuple{2,Int}}}()
global dict_3 = Dict{Int32,Set{SVector{2,Int32}}}()
global adj2v_1 = DT.Adjacent2Vertex{Int,Set{NTuple{2,Int}}}()
global adj2v_2 = DT.Adjacent2Vertex{Int,Vector{NTuple{2,Int}}}()
global adj2v_3 = DT.Adjacent2Vertex{Int32,Set{SVector{2,Int32}}}()

@testset "Constructors and getters" begin
    @test adj2v_1.adjacent2vertex == dict_1
    @test adj2v_2.adjacent2vertex == dict_2
    @test adj2v_3.adjacent2vertex == dict_3
    @test adj2v_1 == DT.Adjacent2Vertex(dict_1)
    @test adj2v_2 == DT.Adjacent2Vertex(dict_2)
    @test adj2v_3 == DT.Adjacent2Vertex(dict_3)
    @test get_adjacent2vertex(adj2v_1) == dict_1
    @test get_adjacent2vertex(adj2v_2) == dict_2
    @test get_adjacent2vertex(adj2v_3) == dict_3
end

global dict_1 = Dict(1 => Set(((1, 2), (3, 4), (10, 15), (2, 5))),
    2 => Set(((5, 7), (10, 14), (2, 3), (5, 9))),
    3 => Set(((10, 25), (23, 29))))
global dict_2 = Dict(1 => [(1, 2), (3, 4), (10, 15), (2, 5)],
    2 => [(5, 7), (10, 14), (2, 3), (5, 9)],
    3 => [(10, 25), (23, 29)])
global dict_3 = Dict{Int32,Set{SVector{2,Int32}}}(1 => Set{SVector{2,Int32}}((@SVector[1, 2],
        @SVector[3, 4],
        @SVector[10, 15],
        @SVector[2, 5])),
    2 => Set{SVector{2,Int32}}((@SVector[5, 7],
        @SVector[10, 14],
        @SVector[2, 3],
        @SVector[5, 9])),
    3 => Set{SVector{2,Int32}}((@SVector[10, 25],
        @SVector[23, 29])))
global adj2v_1 = DT.Adjacent2Vertex(dict_1)
global adj2v_2 = DT.Adjacent2Vertex(dict_2)
global adj2v_3 = DT.Adjacent2Vertex(dict_3)

@testset "Adding and deleting edges, triangles" begin
    for adj2v in (adj2v_1, adj2v_2, adj2v_3)
        w1 = get_adjacent2vertex(adj2v, 1)
        w2 = get_adjacent2vertex(adj2v, 2)
        w3 = get_adjacent2vertex(adj2v, 3)
        if adj2v === adj2v_1
            @test w1 == Set(((1, 2), (3, 4), (10, 15), (2, 5)))
            @test w2 == Set(((5, 7), (10, 14), (2, 3), (5, 9)))
            @test w3 == Set(((10, 25), (23, 29)))
        elseif adj2v === adj2v_2
            @test w1 == [(1, 2), (3, 4), (10, 15), (2, 5)]
            @test w2 == [(5, 7), (10, 14), (2, 3), (5, 9)]
            @test w3 == [(10, 25), (23, 29)]
        elseif adj2v === adj2v_3
            @test w1 ==
                  Set{SVector{2,Int32}}((@SVector[1, 2], @SVector[3, 4], @SVector[10, 15],
                @SVector[2, 5]))
            @test w2 ==
                  Set{SVector{2,Int32}}((@SVector[5, 7], @SVector[10, 14], @SVector[2, 3],
                @SVector[5, 9]))
            @test w3 == Set{SVector{2,Int32}}((@SVector[10, 25], @SVector[23, 29]))
        end

        DT.add_adjacent2vertex!(adj2v, 1, 23, 50)
        DT.add_adjacent2vertex!(adj2v, 3, 38, 173)
        if adj2v === adj2v_1
            @test w1 == Set(((1, 2), (3, 4), (10, 15), (2, 5), (23, 50)))
            @test w2 == Set(((5, 7), (10, 14), (2, 3), (5, 9)))
            @test w3 == Set(((10, 25), (23, 29), (38, 173)))
        elseif adj2v === adj2v_2
            @test w1 == [(1, 2), (3, 4), (10, 15), (2, 5), (23, 50)]
            @test w2 == [(5, 7), (10, 14), (2, 3), (5, 9)]
            @test w3 == [(10, 25), (23, 29), (38, 173)]
        elseif adj2v === adj2v_3
            @test w1 ==
                  Set{SVector{2,Int32}}((@SVector[1, 2], @SVector[3, 4], @SVector[10, 15],
                @SVector[2, 5], @SVector[23, 50]))
            @test w2 ==
                  Set{SVector{2,Int32}}((@SVector[5, 7], @SVector[10, 14], @SVector[2, 3],
                @SVector[5, 9]))
            @test w3 == Set{SVector{2,Int32}}((@SVector[10, 25], @SVector[23, 29],
                @SVector[38, 173]))
        end

        if !(adj2v === adj2v_2)
            DT.delete_adjacent2vertex!(adj2v, 1, 3, 4)
            DT.delete_adjacent2vertex!(adj2v, 2, 5, 9)
            if adj2v === adj2v_1
                @test w1 == Set(((1, 2), (10, 15), (2, 5), (23, 50)))
                @test w2 == Set(((5, 7), (10, 14), (2, 3)))
                @test w3 == Set(((10, 25), (23, 29), (38, 173)))
            elseif adj2v === adj2v_2
                @test w1 == [(1, 2), (10, 15), (2, 5), (23, 50)]
                @test w2 == [(5, 7), (10, 14), (2, 3)]
                @test w3 == [(10, 25), (23, 29), (38, 173)]
            elseif adj2v === adj2v_3
                @test w1 ==
                      Set{SVector{2,Int32}}((@SVector[1, 2], @SVector[10, 15], @SVector[2, 5],
                    @SVector[23, 50]))
                @test w2 ==
                      Set{SVector{2,Int32}}((@SVector[5, 7], @SVector[10, 14], @SVector[2, 3]))
                @test w3 == Set{SVector{2,Int32}}((@SVector[10, 25], @SVector[23, 29],
                    @SVector[38, 173]))
            end
        end

        DT.delete_adjacent2vertex!(adj2v, 3)
        @test 3 ∉ keys(adj2v.adjacent2vertex)
        if adj2v === adj2v_1
            @test w1 == Set(((1, 2), (10, 15), (2, 5), (23, 50)))
            @test w2 == Set(((5, 7), (10, 14), (2, 3)))
        elseif adj2v === adj2v_2
            @test w1 == [(1, 2), (3, 4), (10, 15), (2, 5), (23, 50)]
            @test w2 == [(5, 7), (10, 14), (2, 3), (5, 9)]
        elseif adj2v === adj2v_3
            @test w1 ==
                  Set{SVector{2,Int32}}((@SVector[1, 2], @SVector[10, 15], @SVector[2, 5],
                @SVector[23, 50]))
            @test w2 ==
                  Set{SVector{2,Int32}}((@SVector[5, 7], @SVector[10, 14], @SVector[2, 3]))
        end

        DT.add_triangle!(adj2v, 51, 52, 53)
        @test DT.contains_edge(51, 52, get_adjacent2vertex(adj2v, 53))
        @test DT.contains_edge(52, 53, get_adjacent2vertex(adj2v, 51))
        @test DT.contains_edge(53, 51, get_adjacent2vertex(adj2v, 52))

        DT.add_triangle!(adj2v, (60, 61, 62))
        @test DT.contains_edge(60, 61, get_adjacent2vertex(adj2v, 62))
        @test DT.contains_edge(61, 62, get_adjacent2vertex(adj2v, 60))
        @test DT.contains_edge(62, 60, get_adjacent2vertex(adj2v, 61))

        T1 = (120, 125, 132)
        T2 = (501, 502, 591)
        T3 = (6019, 919, 821)
        DT.add_triangle!.(Ref(adj2v), (T1, T2, T3))
        for T in (T1, T2, T3)
            i, j, k = T
            @test DT.contains_edge(i, j, get_adjacent2vertex(adj2v, k))
            @test DT.contains_edge(j, k, get_adjacent2vertex(adj2v, i))
            @test DT.contains_edge(k, i, get_adjacent2vertex(adj2v, j))
        end

        if !(adj2v === adj2v_2)
            DT.delete_triangle!(adj2v, 60, 61, 62)
            @test !DT.contains_edge(60, 61, get_adjacent2vertex(adj2v, 62))
            @test !DT.contains_edge(61, 62, get_adjacent2vertex(adj2v, 60))
            @test !DT.contains_edge(62, 60, get_adjacent2vertex(adj2v, 61))

            DT.delete_triangle!(adj2v, T2)
            @test !DT.contains_edge(501, 502, get_adjacent2vertex(adj2v, 591))
            @test !DT.contains_edge(502, 591, get_adjacent2vertex(adj2v, 501))
            @test !DT.contains_edge(591, 501, get_adjacent2vertex(adj2v, 502))

            DT.delete_triangle!.(Ref(adj2v), (T1, T3))
            for T in (T1, T2)
                i, j, k = T
                @test !DT.contains_edge(i, j, get_adjacent2vertex(adj2v, k))
                @test !DT.contains_edge(j, k, get_adjacent2vertex(adj2v, i))
                @test !DT.contains_edge(k, i, get_adjacent2vertex(adj2v, j))
            end
        end
    end
end

@testset "Seeing if Adjacent2Vertex is empty and clearing empty sets" begin
    adj2v = DT.Adjacent2Vertex{Int,Set{NTuple{2,Int}}}()
    DT.add_adjacent2vertex!(adj2v, 2, 5, 7)
    DT.add_adjacent2vertex!(adj2v, 2, 7, 13)
    DT.add_adjacent2vertex!(adj2v, 13, 5, 23)
    adj2v_clean = deepcopy(adj2v)
    DT.add_adjacent2vertex!(adj2v, 26, 2, 10)
    @test adj2v ≠ adj2v_clean
    DT.delete_adjacent2vertex!(adj2v, 26, 2, 10)
    @test adj2v ≠ adj2v_clean
    @test isempty(get_adjacent2vertex(adj2v, 26))
    DT.clear_empty_keys!(adj2v)
    @test adj2v == adj2v_clean
end

@testset "empty!" begin
    tri = triangulate(rand(2, 50))
    adj2v = get_adjacent2vertex(tri)
    @test !isempty(get_adjacent2vertex(adj2v))
    empty!(adj2v)
    @test isempty(get_adjacent2vertex(adj2v))
end