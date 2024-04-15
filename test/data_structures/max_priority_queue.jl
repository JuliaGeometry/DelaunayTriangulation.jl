using Test, DataStructures
using ..DelaunayTriangulation
const DT = DelaunayTriangulation
include("../helper_functions.jl")

@test DT.hparent(2) == 1
@test DT.hparent(3) == 1
@test DT.hparent(27) == 27 ÷ 2
@test DT.hleft(3) == 6
@test DT.hleft(4) == 8
@test DT.hright(3) == 7
@test DT.hright(4) == 9
@test DT.hchildren(5) == (10, 11)
@test DT.hchildren(10) == (20, 21)
@test !DT.is_root(2)
@test DT.is_root(1)

pmax = 1000

for n in [1:10; (11:50:10000)]
    if n ≠ 11
        r = rand(1:pmax, n)
        ks, vs = 1:n+1, rand(1:pmax, n)
    else
        r = [670, 438, 251, 527, 40, 428, 913, 922, 756, 344, 562]
        ks, vs = 1:n+1, [591, 242, 443, 24, 640, 173, 83, 409, 197, 125, 352]
    end
    priorities = Dict(zip(ks, vs))
    queue = DT.MaxPriorityQueue(priorities)
    dt_queue = PriorityQueue{Int,Int}(Base.Reverse, priorities)
    if n == 11
        @test sprint(show, MIME"text/plain"(), queue) == "DelaunayTriangulation.MaxPriorityQueue{Int64, Int64} with 11 entries:\n  5 => 640\n  1 => 591\n  3 => 443\n  8 => 409\n  11 => 352\n  2 => 242\n  9 => 197\n  6 => 173\n  10 => 125\n  7 => 83\n  4 => 24"
    end
    @test queue == dt_queue
    @test length(queue) == length(dt_queue)
    @test !isempty(queue)
    for k in keys(queue)
        @test haskey(queue, k) && haskey(dt_queue, k)
    end
    @test first(queue).second == first(dt_queue).second
    queue[27] = 5000.0
    dt_queue[27] = 5000.0
    @test queue == dt_queue
    for _ in 1:n
        k = rand(1:n)
        v = rand(1:pmax)
        queue[k] = v
        dt_queue[k] = v
        @test queue[k] == dt_queue[k]
    end
    @test queue == dt_queue
    pairs = Pair{Int,Int}[]
    dt_pairs = Pair{Int,Int}[]
    while !isempty(queue)
        pair = popfirst!(queue)
        dt_pair = dequeue_pair!(dt_queue)
        push!(pairs, pair)
        push!(dt_pairs, dt_pair)
    end
    @test isempty(queue) && isempty(dt_queue)
    _compare_pairs(pairs, dt_pairs)
end