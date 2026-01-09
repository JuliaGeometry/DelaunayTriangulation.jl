using ..DelaunayTriangulation
using Test
using Random
using StableRNGs
const DT = DelaunayTriangulation

S = [1, 6, 15, 13, 23, 71]
rng = StableRNG(123)
_rng = StableRNG(123)
list = DT.ShuffledPolygonLinkedList(S; rng)
manual_next = [mod1(i + 1, length(S)) for i in eachindex(S)]
manual_prev = [mod1(i - 1, length(S)) for i in eachindex(S)]
manual_shuffled_indices = collect(eachindex(S))
shuffle!(_rng, manual_shuffled_indices)
@test list.next == manual_next
@test list.prev == manual_prev
@test list.shuffled_indices == manual_shuffled_indices
@test list.k == length(S)
@test list.S == S

DT.reset!(list; rng)
@test list.next == manual_next
@test list.prev == manual_prev
shuffle!(_rng, manual_shuffled_indices)
@test list.shuffled_indices == manual_shuffled_indices
@test list.k == length(S)
@test list.S == S == [1, 6, 15, 13, 23, 71] # careful about aliasing

π = manual_shuffled_indices
for i in eachindex(S)
    πᵢ = π[i]
    @test DT.get_triplet(list, i) == (S[πᵢ], S[manual_next[πᵢ]], S[manual_prev[πᵢ]])
end

@test_throws AssertionError DT.ShuffledPolygonLinkedList([1], [1, 2], [1, 2], 2, [1, 2])
@test_throws AssertionError DT.ShuffledPolygonLinkedList([1, 2], [1], [1, 2], 2, [1, 2])
@test_throws AssertionError DT.ShuffledPolygonLinkedList([1, 2], [1, 2], [1], 2, [1, 2])
@test_throws AssertionError DT.ShuffledPolygonLinkedList([1, 2], [1, 2], [1, 2], 1, [1, 2])
@test_throws AssertionError DT.ShuffledPolygonLinkedList([1, 2], [1, 2], [1, 2], 2, [1, 2, 3])

next = [2, 3, 4, 5, 6, 1]
prev = [6, 1, 2, 3, 4, 5]
shuffled_indices = [3, 1, 5, 6, 2, 4]
k = 6
S = [1, 6, 15, 13, 23, 71]
list = DT.ShuffledPolygonLinkedList(next, prev, shuffled_indices, k, S)
π₂ = shuffled_indices[2]
DT.delete_vertex!(list, 2)
@test next[prev[π₂]] == next[π₂]
@test prev[next[π₂]] == prev[π₂]

π = shuffled_indices
_π = copy(π)
DT.swap_permutation!(list, 1, 2)
@test π[1] == _π[2]
@test π[2] == _π[1]
DT.swap_permutation!(list, 3, 5)
@test π[3] == _π[5]
@test π[5] == _π[3]
