using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using Test
using StaticArraysCore
using Random

global i = 17
global j = 5
global e1 = (i, j)
global e2 = [i, j]
global e3 = SVector{2,Int32}((i, j))

throw_f = expr -> @static if VERSION < v"1.8"
    ErrorException(expr)
else
    expr
end

@testset "Individual edges" begin
    @testset "Constructing an edge" begin
        @test_throws throw_f("The construct_edge function has not been defined for the type String.") DT.construct_edge(String, 1, 2)
        for e in (e1, e2, e3)
            @test DT.construct_edge(typeof(e), i, j) == e
        end
    end

    @testset "Getting indices" begin
        @test_throws throw_f("The initial function has not been defined for the type DataType.") initial(String)
        for e in (e1, e2, e3)
            @test initial(e) == i
        end
        @test_throws throw_f("The terminal function has not been defined for the type DataType.") terminal(String)
        for e in (e1, e2, e3)
            @test terminal(e) == j
        end
        for e in (e1, e2, e3)
            @test DT.edge_indices(e) == (e[1], e[2])
        end
    end

    @testset "Reversing an edge" begin
        for e in (e1, e2, e3)
            @test DT.reverse_edge(e) == DT.construct_edge(typeof(e), reverse(e)...)
        end
    end

    @testset "Comparing unoriented edges" begin
        @test DT.compare_unoriented_edge((1, 2), (2, 1))
        @test !DT.compare_unoriented_edge((1, 5), (7, 3))
        @test DT.compare_unoriented_edge((1, 7), (1, 7))
    end
end

global es1 = Set{typeof(e1)}(((1, 3), (4, 1), (10, 1), (3, 9), (5, 3)))
global es2 = Set{typeof(e2)}(([1, 3], [4, 1], [10, 1], [3, 9], [5, 3]))
global es3 = Set{typeof(e3)}((SVector{2,Int32}((1, 3)),
    SVector{2,Int32}((4, 1)),
    SVector{2,Int32}((10, 1)),
    SVector{2,Int32}((3, 9)),
    SVector{2,Int32}((5, 3))))

@testset "Collection of edges" begin
    @testset "Initialising a collection of edges" begin
        @test_throws throw_f("The initialise_edges function has not been defined for the type String.") DT.initialise_edges(String)
        for es in (es1, es2, es3)
            F = eltype(es)
            @test DT.initialise_edges(Set{F}) == Set{F}()
        end
        @test DT.initialise_edges(Vector{NTuple{2,Int}}) == Vector{NTuple{2,Int}}()
    end

    @testset "Getting the type of edges in a collection" begin
        @test_throws throw_f("The edge_type function has not been defined for the type String.") DT.edge_type(String)
        for es in (es1, es2, es3)
            F = eltype(es)
            @test DT.edge_type(typeof(es)) == F
        end
        @test DT.edge_type(Vector{NTuple{2,Int}}) == NTuple{2,Int}
    end

    @testset "Number of edges" begin
        @test_throws throw_f("The num_edges function has not been defined for the type DataType.") DT.num_edges(String)
        for es in (es1, es2, es3)
            @test num_edges(es) == length(es) == 5
        end
        @test num_edges([(1, 2), (2, 3), (4, 5), (10, 11)]) == 4
    end

    @testset "Seeing if a collection contains an edge" begin
        @test_throws throw_f("The contains_edge function has not been defined for the type Tuple{Int64, Int64}.") DT.contains_edge(String, (1, 2))
        for es in (es1, es2, es3)
            for e in es
                @test DT.contains_edge(e, es)
                @test DT.contains_edge(e[1], e[2], es)
            end
            f1 = DT.construct_edge(eltype(es), 17, 29)
            f2 = DT.construct_edge(eltype(es), 30, 50)
            @test !DT.contains_edge(f1, es)
            @test !DT.contains_edge(f2, es)
            @test !DT.contains_edge(17, 29, es)
            @test !DT.contains_edge(30, 50, es)
        end
        @test DT.contains_edge(1, 7, [[2, 3], [1, 7]])
        @test !DT.contains_edge(10, 5, [[50, 10]])
    end

    @testset "Adding to a collection of edges" begin
        @test_throws throw_f("The add_to_edges! function has not been defined for the type DataType.") DT.add_to_edges!(String, (1, 2))
        for es in (es1, es2, es3)
            f1 = DT.construct_edge(eltype(es), 17, 29)
            f2 = DT.construct_edge(eltype(es), 30, 50)
            f3 = DT.construct_edge(eltype(es), 50, 20)
            DT.add_to_edges!(es, f1)
            DT.add_edge!(es, f2, f3)
            @test length(es) == 8
            @test f1 ∈ es && DT.contains_edge(f1, es)
            @test f2 ∈ es && DT.contains_edge(f2, es)
            @test f3 ∈ es && DT.contains_edge(f3, es)
        end
    end

    @testset "Deleting from a collection of edges" begin
        @test_throws throw_f("The delete_from_edges! function has not been defined for the type DataType.") DT.delete_from_edges!(String, (1, 2))
        for es in (es1, es2, es3)
            f1 = DT.construct_edge(eltype(es), 17, 29)
            f2 = DT.construct_edge(eltype(es), 30, 50)
            f3 = DT.construct_edge(eltype(es), 50, 20)
            DT.delete_from_edges!(es, f1)
            DT.delete_edge!(es, f2, f3)
            @test length(es) == 5
            @test f1 ∉ es && !DT.contains_edge(f1, es)
            @test f2 ∉ es && !DT.contains_edge(f2, es)
            @test f3 ∉ es && !DT.contains_edge(f3, es)
        end
        Es = [(1, 2), (5, 10), (17, 23)]
        DT.delete_edge!(Es, (17, 23))
        @test Es == [(1, 2), (5, 10)]
        Es = [(1, 2), (3, 4), (10, 15), (23, 10), (17, 2)]
        DT.delete_edge!(Es, (1, 2), (23, 10))
        @test Es == [(3, 4), (10, 15), (17, 2)]
    end

    @testset "Iterating over each edge in a collection" begin
        @test_throws throw_f("The each_edge function has not been defined for the type DataType.") each_edge(String)
        for es in (es1, es2, es3)
            @test each_edge(es) == es
        end
        E = [1 2; 10 15; 23 5]'
        @test each_edge(E) == eachcol(E)
        E = [(1, 2), (2, 3), (4, 5)]
        @test each_edge(E) == E
    end

    @testset "Getting a random edge from a collection" begin
        for es in (es1, es2, es3)
            local i
            i = rand(1:29991)
            Random.seed!(i)
            u = rand(es)
            Random.seed!(i)
            v = DT.random_edge(es)
            @test u == v
        end
    end

    @testset "Seeing if a collection of edges is empty" begin
        @test DT.is_empty([])
        @test !DT.is_empty([1, 2, 3])
        @test DT.is_empty(Set{NTuple{2,Int}}())
        @test !DT.is_empty(Set(((1, 2), (3, 4))))
    end
end
