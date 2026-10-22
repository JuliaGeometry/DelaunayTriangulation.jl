using Test, DataStructures
using ..DelaunayTriangulation
using StableRNGs
const DT = DelaunayTriangulation
queue = DT.Queue{Float64}()
@test isempty(queue)
@test length(queue) == 0
@test eltype(queue) == Float64
ds_queue = DataStructures.Queue{Float64}()
els = Float64[]
ds_els = Float64[]
rng = StableRNG(1)
for k in 1:100
    r = randn(rng)
    push!(queue, r)
    push!(ds_queue, r)
end
for k in 1:100
    push!(els, popfirst!(queue))
    push!(ds_els, popfirst!(ds_queue))
end
@test isempty(queue) && isempty(ds_queue)
@test els == ds_els

queue = DT.Queue{Int}()
rng = StableRNG(2)
points = rand(rng, 2, 50)
DT.enqueue_all!(queue, DT.each_point_index(points))
@test queue.data == 1:50

q1 = DT.Queue{Int}()
q2 = DT.Queue{Int}()
@test q1 == q2
push!(q1, 1)
@test q1 ≠ q2
push!(q2, 1)
@test q1 == q2

@testset "copy/deepcopy" begin
    queue = DT.Queue{Int}([1, 2, 3, 4])
    queue2 = copy(queue)
    @test queue == queue2 && !(queue === queue2)
    queue2 = deepcopy(queue)
    @test queue == queue2 && !(queue === queue2)
end