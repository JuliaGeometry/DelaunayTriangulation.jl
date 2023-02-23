using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using Test
using StaticArraysCore
using Random

i, j = 17, 5
e1 = (i, j)
e2 = [i, j]
e3 = SVector{2,Int32}((i, j))

@test_throws "The" DT.construct_edge(String, 1, 2)
for e in (e1, e2, e3)
    @test DT.construct_edge(typeof(e), i, j) == e
end

@test_throws "The" initial(String)
for e in (e1, e2, e3)
    @test initial(e) == i
end

@test_throws "The" terminal(String)
for e in (e1, e2, e3)
    @test terminal(e) == j
end

for e in (e1, e2, e3)
    @test DT.edge_indices(e) == (e[1], e[2])
end

es1 = Set{typeof(e1)}(((1, 3), (4, 1), (10, 1), (3, 9), (5, 3)))
es2 = Set{typeof(e2)}(([1, 3], [4, 1], [10, 1], [3, 9], [5, 3]))
es3 = Set{typeof(e3)}((SVector{2,Int32}((1, 3)),
                       SVector{2,Int32}((4, 1)),
                       SVector{2,Int32}((10, 1)),
                       SVector{2,Int32}((3, 9)),
                       SVector{2,Int32}((5, 3))))

@test_throws "The" DT.initialise_edges(String)
for es in (es1, es2, es3)
    F = eltype(es)
    @test DT.initialise_edges(Set{F}) == Set{F}()
end
@test DT.initialise_edges(Vector{NTuple{2,Int64}}) == Vector{NTuple{2,Int64}}()

@test_throws "The" DT.edge_type(String)
for es in (es1, es2, es3)
    F = eltype(es)
    @test DT.edge_type(typeof(es)) == F
end
@test DT.edge_type(Vector{NTuple{2,Int64}}) == NTuple{2,Int64}

@test_throws "The" DT.num_edges(String)
for es in (es1, es2, es3)
    @test num_edges(es) == length(es) == 5
end

@test_throws "The" DT.contains_edge(String, (1, 2))
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

@test_throws "The" DT.add_to_edges!(String, (1, 2))
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

@test_throws "The" DT.delete_from_edges!(String, (1, 2))
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

@test_throws "The" each_edge(String)
for es in (es1, es2, es3)
    @test each_edge(es) == es
end
E = [1 2; 10 15; 23 5]'
@test each_edge(E) == eachcol(E)

for es in (es1, es2, es3)
    local i
    i = rand(1:29991)
    Random.seed!(i)
    u = rand(es)
    Random.seed!(i)
    v = DT.random_edge(es)
    @test u == v
end

@test DT.is_empty([])
@test !DT.is_empty([1, 2, 3])
@test DT.is_empty(Set{NTuple{2,Int64}}())
@test !DT.is_empty(Set(((1, 2), (3, 4))))