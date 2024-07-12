using Test, DataStructures
using ..DelaunayTriangulation
const DT = DelaunayTriangulation
queue = DT.Queue{Float64}()
@test isempty(queue)
@test length(queue) == 0
@test eltype(queue) == Float64
dt_queue = DT.Queue{Float64}()
els = Float64[]
dt_els = Float64[]
for k in 1:100
    r = randn()
    push!(queue, r)
    push!(dt_queue, r)
end
for k in 1:100
    push!(els, popfirst!(queue))
    push!(dt_els, popfirst!(dt_queue))
end
@test isempty(queue)
@test els == dt_els


queue = DT.Queue{Int}()
points = rand(2, 50)
DT.enqueue_all!(queue, DT.each_point_index(points))
@test queue.data == 1:50


q1 = DT.Queue{Int}()
q2 = DT.Queue{Int}()
@test q1 == q2
push!(q1, 1)
@test q1 ≠ q2
push!(q2, 1)
@test q1 == q2
