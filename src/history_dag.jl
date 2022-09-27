"""
    HistoryDAG <: AbstractSimpleGraph

Data structure for the Delaunay triangulation. This is a directed acyclic graph 
that is useful for point location.
"""
struct HistoryDAG <: AbstractSimpleGraph
    graph::DirectedGraph{TriangleType}
    function HistoryDAG()
        G = DirectedGraph{TriangleType}()
        forbid_loops!(G)
        TDAG = new(G)
        return TDAG
    end
end

graph(G::TriangulationDAG) = G.graph

for op in (:out_deg, :in_deg, :deg,
    :dual_deg, :in_neighbors,
    :NE, :has, :delete!,
    :elist)
    @eval begin
        SimpleGraphs.$op(G::HistoryDAG) = SimpleGraphs.$op(graph(G))
        SimpleGraphs.$op(G::HistoryDAG, u) = SimpleGraphs.$op(graph(G), u)
        SimpleGraphs.$op(G::HistoryDAG, u, v) = SimpleGraphs.$op(graph(G), u, v)
    end
end

SimpleGraphs.out_neighbors(G::HistoryDAG, u) = graph(G).N[u] # The version in SimpleGraphs is allocating...

SimpleGraphs.add!(G::HistoryDAG, u) = SimpleGraphs.add!(graph(G), u)
"""
    SimpleGraphs.add!(G::HistoryDAG, u)

Add the triangle `u` to the graph `G`. Checks that cyclic permutations of `u`'s indices are not 
in `G` already, returning `false` otherwise.
"""
function SimpleGraphs.add!(G::HistoryDAG, u)
    for i = 1:3
        i₁, i₂, i₃ = i, (i % 3) + 1, ((i + 1) % 3) + 1
        has(G, TriangleType((u[i₁], u[i₂], u[i₃]))) && return false
    end
    return SimpleGraphs.add!(graph(G), u)
end
"""
    SimpleGraphs.add!(G::HistoryDAG, u, v)

Adds an edge between triangle `u ∈ G` and `v ∈ G` to the graph `G`. Checks 
are made for cyclic permutations of `u` and `v`'s indices.
"""
function SimpleGraphs.add!(G::HistoryDAG, u, v)
    for i in 1:3
        for j in 1:3
            i₁, i₂, i₃ = i, (i % 3) + 1, ((i + 1) % 3) + 1
            j₁, j₂, j₃ = j, (j % 3) + 1, ((j + 1) % 3) + 1
            T = TriangleType((u[i₁], u[i₂], u[i₃]))
            V = TriangleType((v[j₁], v[j₂], v[j₃]))
            if has(G, T) && has(G, V)
                SimpleGraphs.add!(graph(G), T, V)
                return true
            end
        end
    end
    return false
end

"""
    find_root(G::HistoryDAG; method=:brute)

Finds the root of the graph `G`, assuming only one root exists and the only node with in-degree zero is this root. 
There are two methods:

    - `method=:brute`: In this case, all vertices are searched until one is found with in-degree zero.
    - `method=:rng`: In this case, a random integer is selected, defniing the node to start at. 
    We then go upwards until the root is found; note that this only works as the graph is acyclic.
"""
function find_root(G::HistoryDAG; method=:brute)
    if method == :brute
        return _find_root_brute(G)
    elseif method == :rng
        return _find_root_rng(G)
    end
end
function _find_root_brute(G::HistoryDAG)
    for (k, v) in graph(G).NN
        length(v) == 0 && return k
    end
end
function _find_root_rng(G::HistoryDAG)
    verts = vlist(graph(G))
    num_verts = length(verts)
    starting_node = verts[rand(1:num_verts)]
    for _ in 1:num_verts
        num_choices = in_deg(G, starting_node)
        num_choices == 0 && return starting_node
        starting_node = in_neighbors(G, starting_node)[rand(1:num_choices)]
    end
end