###################################################
#/
#/
#/ Operations 
#/
#/
###################################################

#
# This script contains functions for performing certain 
# operations on triangulations. The contents are as 
# follows:
#   1. AddTriangle: Method for adding a triangle into a triangulation.
#   2. DeleteTriangle: Method for deleting a triangle from a triangulation.
#   3. FlipEdge: Method for flipping an edge in a triangulation. 
#   4. SplitEdge: Method for splitting a triangle into two along an edge.
#   5. SplitTriangle: Method for splitting a triangle into three about a point in its interior.
#   6. LegaliseEdge: Given an illegal edge, will make it legal by flipping edges in the triangulation.
#

###################################################
#/
#/
#/ AddTriangle
#/
#/
###################################################
function add_triangle!(i, j, k, T::Ts,
    adj::Adjacent{I,E}, adj2v::Adjacent2Vertex{I,Es,E},
    DG::DelaunayGraph{I}; update_ghost_edges=false) where {I,E,Es,Ts}
    V = triangle_type(Ts)
    Tᵢⱼₖ = construct_triangle(V, i, j, k)
    add_triangle!(T, Tᵢⱼₖ)
    ij_bnd = is_boundary_edge(i, j, adj)
    jk_bnd = is_boundary_edge(j, k, adj)
    ki_bnd = is_boundary_edge(k, i, adj)
    num_bnd_edges = count((ij_bnd, jk_bnd, ki_bnd))
    add_edge!(adj, i, j, k)
    add_edge!(adj, j, k, i)
    add_edge!(adj, k, i, j)
    add_edge!(adj2v, i, j, k)
    add_edge!(adj2v, j, k, i)
    add_edge!(adj2v, k, i, j)
    add_neighbour!(DG, k, i, j)
    add_neighbour!(DG, i, j)
    if num_bnd_edges == 1
        add_boundary_edges_single!(i, j, k, ij_bnd, jk_bnd, ki_bnd, T, adj, adj2v, DG; update_ghost_edges)
    elseif num_bnd_edges == 2
        add_boundary_edges_double!(i, j, k, ij_bnd, jk_bnd, ki_bnd, T, adj, adj2v, DG; update_ghost_edges)
    elseif length(T) == 1 # If this is true, then we just have a single triangle, and thus num_bnd_edges should be 3 
        add_boundary_edges_triple!(i, j, k, T, adj, adj2v, DG; update_ghost_edges)
    end
    return nothing
end
function add_triangle!(i, j, k, T, adj, adj2v, DG, HG::HistoryGraph{V}; update_ghost_edges=false) where {V}
    add_triangle!(i, j, k, T, adj, adj2v, DG; update_ghost_edges)
    Tᵢⱼₖ = construct_triangle(V, i, j, k)
    add_triangle!(HG, Tᵢⱼₖ)
    return nothing
end
function add_boundary_edges_single!(i, j, k, ij_bnd, jk_bnd, ki_bnd,
    T::Ts, adj::Adjacent{I,E}, adj2v, DG; update_ghost_edges=false) where {Ts,I,E}
    u, v, w = choose_uvw(ij_bnd, jk_bnd, ki_bnd, i, j, k)
    add_edge!(adj, u, w, I(BoundaryIndex))
    add_edge!(adj, w, v, I(BoundaryIndex))
    add_edge!(adj2v, I(BoundaryIndex), u, w)
    add_edge!(adj2v, I(BoundaryIndex), w, v)
    delete_edge!(adj2v, I(BoundaryIndex), u, v)
    if update_ghost_edges
        #= 
        For the Bowyer-Watson algorithm, we need to handle the ghost handles appropriately. In this case,  
        we are adding two ghost triangles uw∂ and wv∂, coming from the two new boundary edges uw and wv. 
        In addition, we need to make sure we delete the original ghost triangle uv∂. In doing this 
        deletion, though, we need to make sure the edges v∂ and ∂u do not get accidentally deleted.
        Also note that some of these updates are done above.
        =#
        # Ghost edges
        add_edge!(adj, w, I(BoundaryIndex), u)
        add_edge!(adj, I(BoundaryIndex), u, w)
        add_edge!(adj, v, I(BoundaryIndex), w)
        add_edge!(adj, I(BoundaryIndex), w, v)
        add_edge!(adj2v, u, w, I(BoundaryIndex))
        add_edge!(adj2v, w, I(BoundaryIndex), u)
        add_edge!(adj2v, w, v, I(BoundaryIndex))
        add_edge!(adj2v, v, I(BoundaryIndex), w)
        delete_edge!(adj2v, u, v, I(BoundaryIndex))
        delete_edge!(adj2v, v, I(BoundaryIndex), u)
        add_neighbour!(DG, I(BoundaryIndex), w) # Note that u and v are already in DG[∂]
        # Ghost triangles
        V = triangle_type(Ts)
        T1 = construct_triangle(V, u, w, I(BoundaryIndex))
        T2 = construct_triangle(V, w, v, I(BoundaryIndex))
        T3 = construct_triangle(V, u, v, I(BoundaryIndex))
        add_triangle!(T, T1, T2)
        delete_triangle!(T, T3)
    end
    return nothing
end
function add_boundary_edges_double!(i, j, k, ij_bnd, jk_bnd, ki_bnd,
    T::Ts, adj::Adjacent{I,E}, adj2v, DG; update_ghost_edges=false) where {Ts,I,E}
    u, v, w = choose_uvw(!ij_bnd, !jk_bnd, !ki_bnd, i, j, k)
    add_edge!(adj, v, u, I(BoundaryIndex))
    add_edge!(adj2v, I(BoundaryIndex), v, u)
    delete_edge!(adj2v, I(BoundaryIndex), v, w)
    delete_edge!(adj2v, I(BoundaryIndex), w, u)
    if update_ghost_edges
        #= 
        For the Bowyer-Watson algorithm, we need to handle the ghost handles appropriately. In this case,  
        we are adding a single ghost triangle vu∂, coming from the new boundary edge vu. 
        In addition, we need to make sure we delete the two original ghost triangles vw∂ and wu∂. In doing this 
        deletion, though, we need to make sure the edges ∂v and u∂ do not get accidentally deleted. 
        =#
        # Ghost edges
        add_edge!(adj, u, I(BoundaryIndex), v)
        add_edge!(adj, I(BoundaryIndex), v, u)
        delete_edge!(adj, w, I(BoundaryIndex))
        delete_edge!(adj, I(BoundaryIndex), w)
        add_edge!(adj2v, v, u, I(BoundaryIndex))
        add_edge!(adj2v, u, I(BoundaryIndex), v)
        delete_edge!(adj2v, w, u, I(BoundaryIndex))
        delete_edge!(adj2v, u, I(BoundaryIndex), w)
        delete_edge!(adj2v, v, w, I(BoundaryIndex))
        delete_edge!(adj2v, w, I(BoundaryIndex), v)
        delete_neighbour!(DG, I(BoundaryIndex), w) # u and v are still in DG[∂]
        # Ghost triangles
        V = triangle_type(Ts)
        T1 = construct_triangle(V, v, u, I(BoundaryIndex))
        T2 = construct_triangle(V, v, w, I(BoundaryIndex))
        T3 = construct_triangle(V, w, u, I(BoundaryIndex))
        add_triangle!(T, T1)
        delete_triangle!(T, T2, T3)
    end
    return nothing
end
function add_boundary_edges_triple!(i, j, k, T::Ts, adj::Adjacent{I,E},
    adj2v, DG; update_ghost_edges=false) where {Ts,I,E}
    add_edge!(adj, j, i, I(BoundaryIndex))
    add_edge!(adj, i, k, I(BoundaryIndex))
    add_edge!(adj, k, j, I(BoundaryIndex))
    add_edge!(adj2v, I(BoundaryIndex), j, i)
    add_edge!(adj2v, I(BoundaryIndex), i, k)
    add_edge!(adj2v, I(BoundaryIndex), k, j)
    if update_ghost_edges
        #= 
        For the Bowyer-Watson algorithm, we need to handle the ghost handles appropriately. In this case,  
        we are adding three ghost triangles ik∂, kj∂, and ji∂. Since the triangulation is initially empty, 
        we do not need to delete anything here.
        =#
        # Ghost edges
        add_edge!(adj, i, I(BoundaryIndex), j)
        add_edge!(adj, I(BoundaryIndex), j, i)
        add_edge!(adj, j, I(BoundaryIndex), k)
        add_edge!(adj, I(BoundaryIndex), k, j)
        add_edge!(adj, k, I(BoundaryIndex), i)
        add_edge!(adj, I(BoundaryIndex), i, k)
        add_edge!(adj2v, j, i, I(BoundaryIndex))
        add_edge!(adj2v, i, I(BoundaryIndex), j)
        add_edge!(adj2v, k, j, I(BoundaryIndex))
        add_edge!(adj2v, j, I(BoundaryIndex), k)
        add_edge!(adj2v, i, k, I(BoundaryIndex))
        add_edge!(adj2v, k, I(BoundaryIndex), i)
        add_neighbour!(DG, I(BoundaryIndex), i, j, k)
        # Ghost triangles
        V = triangle_type(Ts)
        T1 = construct_triangle(V, j, i, I(BoundaryIndex))
        T2 = construct_triangle(V, i, k, I(BoundaryIndex))
        T3 = construct_triangle(V, k, j, I(BoundaryIndex))
        add_triangle!(T, T1, T2, T3)
    end
    return nothing
end

###################################################
#/
#/
#/ DeleteTriangle
#/
#/
###################################################
function delete_triangle!(i, j, k, T::Ts,
    adj::Adjacent{I,E}, adj2v::Adjacent2Vertex{I,Es,E},
    DG::DelaunayGraph{I}; protect_boundary=false) where {I,E,Es,Ts}
    V = triangle_type(Ts)
    Tᵢⱼₖ = construct_triangle(V, i, j, k)
    delete_triangle!(T, Tᵢⱼₖ)
    delete_edge!(adj, i, j)
    delete_edge!(adj, j, k)
    delete_edge!(adj, k, i)
    delete_edge!(adj2v, i, j, k)
    delete_edge!(adj2v, j, k, i)
    delete_edge!(adj2v, k, i, j)
    ji_bnd = is_boundary_edge(j, i, adj)
    ik_bnd = is_boundary_edge(i, k, adj)
    kj_bnd = is_boundary_edge(k, j, adj)
    num_bnd_edges = !protect_boundary ? count((ji_bnd, ik_bnd, kj_bnd)) : 0 # If we are protecting the boundary, setting this count to zero will do that 
    ji_exists = edge_exists(j, i, adj) # Could also write this as ji_is_bnd || edge_exists(j, i, adj),
    ik_exists = edge_exists(i, k, adj) #    but this way makes the code below simpler.
    kj_exists = edge_exists(k, j, adj)
    protect_ji_neighbour = ji_exists && !ji_bnd
    protect_ik_neighbour = ik_exists && !ik_bnd
    protect_kj_neighbour = kj_exists && !kj_bnd
    !protect_ji_neighbour && delete_neighbour!(DG, i, j)
    !protect_ik_neighbour && delete_neighbour!(DG, k, i)
    !protect_kj_neighbour && delete_neighbour!(DG, j, k)
    if num_bnd_edges == 1
        delete_boundary_edges_single!(i, j, k, ji_bnd, ik_bnd, kj_bnd, adj, adj2v)
    elseif num_bnd_edges == 2
        delete_boundary_edges_double!(i, j, k, ji_bnd, ik_bnd, kj_bnd, adj, adj2v)
    elseif num_bnd_edges == 3 # length(T) == 0
        delete_boundary_edges_triple!(i, j, k, adj, adj2v)
    end
    return nothing
end
function delete_boundary_edges_single!(i, j, k, ji_bnd, ik_bnd, kj_bnd,
    adj::Adjacent{I,E}, adj2v::Adjacent2Vertex) where {I,E}
    u, v, w = choose_uvw(ji_bnd, kj_bnd, ik_bnd, i, j, k)
    delete_edge!(adj, v, u)
    delete_edge!(adj2v, I(BoundaryIndex), v, u)
    add_edge!(adj, v, w, I(BoundaryIndex))
    add_edge!(adj, w, u, I(BoundaryIndex))
    add_edge!(adj2v, I(BoundaryIndex), v, w)
    add_edge!(adj2v, I(BoundaryIndex), w, u)
    return nothing
end
function delete_boundary_edges_double!(i, j, k, ji_bnd, ik_bnd, kj_bnd,
    adj::Adjacent{I,E}, adj2v::Adjacent2Vertex) where {I,E}
    u, v, w = choose_uvw(!ji_bnd, !kj_bnd, !ik_bnd, i, j, k)
    delete_edge!(adj, u, w)
    delete_edge!(adj, w, v)
    delete_edge!(adj2v, I(BoundaryIndex), u, w)
    delete_edge!(adj2v, I(BoundaryIndex), w, v)
    add_edge!(adj, u, v, I(BoundaryIndex))
    add_edge!(adj2v, I(BoundaryIndex), u, v)
    return nothing
end
function delete_boundary_edges_triple!(i, j, k, adj::Adjacent{I,E}, adj2v) where {I,E}
    delete_edge!(adj, k, j)
    delete_edge!(adj, j, i)
    delete_edge!(adj, i, k)
    delete_edge!(adj2v, I(BoundaryIndex), k, j)
    delete_edge!(adj2v, I(BoundaryIndex), j, i)
    delete_edge!(adj2v, I(BoundaryIndex), i, k)
    return nothing
end

###################################################
#/
#/
#/ FlipEdge
#/
#/
###################################################
function flip_edge!(i, j, k, ℓ, T, adj, adj2v, DG)
    delete_triangle!(i, k, j, T, adj, adj2v, DG; protect_boundary=true)
    delete_triangle!(i, j, ℓ, T, adj, adj2v, DG; protect_boundary=true)
    add_triangle!(ℓ, k, j, T, adj, adj2v, DG)
    add_triangle!(ℓ, i, k, T, adj, adj2v, DG)
    return nothing
end
function flip_edge!(i, j, T, adj, adj2v, DG)
    ℓ = get_edge(adj, i, j)
    k = get_edge(adj, j, i)
    flip_edge!(i, j, k, ℓ, T, adj, adj2v, DG)
    return nothing
end
function flip_edge!(i, j, k, ℓ, T, adj, adj2v, DG, HG::HistoryGraph)
    flip_edge!(i, j, k, ℓ, T, adj, adj2v, DG)
    flip_edge!(i, j, k, ℓ, HG)
    return nothing
end
function flip_edge!(i, j, T, adj, adj2v, DG, HG::HistoryGraph)
    ℓ = get_edge(adj, i, j)
    k = get_edge(adj, j, i)
    flip_edge!(i, j, k, ℓ, T, adj, adj2v, DG, HG)
    return nothing
end
function flip_edge!(i, j, k, ℓ, HG::HistoryGraph{V}) where {V}
    Tᵢₖⱼ = construct_triangle(V, i, k, j)
    Tᵢⱼₗ = construct_triangle(V, i, j, ℓ)
    Tₗₖⱼ = construct_triangle(V, ℓ, k, j)
    Tₗᵢₖ = construct_triangle(V, ℓ, i, k)
    add_triangle!(HG, Tₗₖⱼ, Tₗᵢₖ)
    add_edge!(HG, Tᵢₖⱼ, Tₗₖⱼ, Tₗᵢₖ)
    add_edge!(HG, Tᵢⱼₗ, Tₗₖⱼ, Tₗᵢₖ)
    return nothing
end

###################################################
#/
#/
#/ SplitEdge
#/
#/
###################################################
function split_edge!(i, j, r, T, adj, adj2v, DG)
    k = get_edge(adj, i, j)
    delete_triangle!(i, j, k, T, adj, adj2v, DG; protect_boundary=!is_boundary_edge(j, i, adj))
    add_triangle!(i, r, k, T, adj, adj2v, DG)
    add_triangle!(r, j, k, T, adj, adj2v, DG)
    return nothing
end
function split_edge!(i, j, r, T, adj, adj2v, DG, HG::HistoryGraph{V}) where {V}
    split_edge!(i, j, r, T, adj, adj2v, DG)
    k = get_edge(adj, i, j)
    Tᵢⱼₖ = construct_triangle(V, i, j, k)
    Tᵢᵣₖ = construct_triangle(V, i, r, k)
    Tᵣⱼₖ = construct_triangle(V, r, j, k)
    add_triangle!(HG, Tᵢᵣₖ, Tᵣⱼₖ)
    add_edge!(HG, Tᵢⱼₖ, Tᵢᵣₖ, Tᵣⱼₖ)
    return nothing
end

###################################################
#/
#/
#/ SplitTriangle
#/
#/
###################################################
function split_triangle!(i, j, k, r, T, adj, adj2v, DG)
    delete_triangle!(i, j, k, T, adj, adj2v, DG; protect_boundary=true)
    add_triangle!(i, j, r, T, adj, adj2v, DG)
    add_triangle!(j, k, r, T, adj, adj2v, DG)
    add_triangle!(k, i, r, T, adj, adj2v, DG)
    return nothing
end
function split_triangle!(i, j, k, r, T, adj, adj2v, DG, HG::HistoryGraph)
    split_triangle!(i, j, k, r, T, adj, adj2v, DG)
    split_triangle!(i, j, k, r, HG)
    return nothing
end
function split_triangle!(i, j, k, r, HG::HistoryGraph{V}) where {V}
    Tᵢⱼᵣ = construct_triangle(V, i, j, r)
    Tⱼₖᵣ = construct_triangle(V, j, k, r)
    Tₖᵢᵣ = construct_triangle(V, k, i, r)
    Tᵢⱼₖ = construct_triangle(V, i, j, k)
    add_triangle!(HG, Tᵢⱼᵣ, Tⱼₖᵣ, Tₖᵢᵣ)
    add_edge!(HG, Tᵢⱼₖ, Tᵢⱼᵣ, Tⱼₖᵣ, Tₖᵢᵣ)
end

###################################################
#/
#/
#/ LegaliseEdge
#/
#/
###################################################
function legalise_edge!(i, j, r, T, HG, adj, adj2v, DG, pts)
    if !islegal(i, j, adj, pts)
        e = get_edge(adj, j, i)
        flip_edge!(i, j, e, r, T, adj, adj2v, DG, HG)
        legalise_edge!(i, e, r, T, HG, adj, adj2v, DG, pts)
        legalise_edge!(e, j, r, T, HG, adj, adj2v, DG, pts)
    end
    return nothing
end

