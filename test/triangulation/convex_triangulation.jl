using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using Random
using Test
using DataStructures
using CairoMakie

save_path = basename(pwd()) == "test" ? "figures" : "test/figures"

a = [4.0, 8.0]
b = [7.0, 7.0]
c = [8.0, 3.0]
d = [5.0, 1.0]
e = [2.0, 4.0]
f = [1.0, 7.0]
pts = [a, b, c, d, e, f] |> reverse # counterclockwise
rng = Random.default_rng()
S = convert(Vector{I}, collect(each_point_index(pts)))
perm_S = shuffle(rng, S)
k = length(S)
next = [mod1(i + 1, k) for i in eachindex(S)]
prev = [mod1(i - 1, k) for i in eachindex(S)]
for i in k:-1:4
    next[prev[perm_S[i]]] = next[perm_S[i]]
    prev[next[perm_S[i]]] = prev[perm_S[i]]
end


a = [4.0, 8.0]
b = [7.0, 7.0]
c = [8.0, 3.0]
d = [5.0, 1.0]
e = [2.0, 4.0]
f = [1.0, 7.0]
pts = [a, b, c, d, e, f] |> reverse # counterclockwise
rng = Random.default_rng()
I = Int64
S = convert(Vector{I}, collect(each_point_index(pts)))
k = length(S)
next = zeros(I, k)
prev = zeros(I, k)
