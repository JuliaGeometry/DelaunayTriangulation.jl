"""
    triangulate_rectangle(a, b, c, d, nx, ny; kwargs...) -> Triangulation

Triangulates the rectangle `[a, b] × [c, d]`.

# Arguments 
- `a`: The minimum `x`-coordinate. 
- `b`: The maximum `x`-coordinate.
- `c`: The minimum `y`-coordinate.
- `d`: The maximum `y`-coordinate.
- `nx`: The number of points in the `x`-direction.
- `ny`: The number of points in the `y`-direction.

# Keyword Arguments 
- `single_boundary=false`: If `true`, then the boundary nodes are stored as a contiguous section. Otherwise, the boundary is split into four sections, in the order 
   bottom, right, top, left.
- `delete_ghosts=false`: If `true`, then the ghost triangles are deleted. Otherwise, they are kept.
- `IntegerType::Type{I}=Int`: The type of the vertices. 
- `EdgeType::Type{E}=NTuple{2,IntegerType}`: The type of the edges.
- `TriangleType::Type{V}=NTuple{3,IntegerType}`: The type of the triangles.
- `EdgesType::Type{Es}=Set{EdgeType}`: The type of the edges container.
- `TrianglesType::Type{Ts}=Set{TriangleType}`: The type of the triangles container.

# Outputs 
- `tri`: The triangulation of the rectangle.
"""
@inline function triangulate_rectangle(
        a, b, c, d, nx, ny;
        single_boundary = false,
        delete_ghosts = false,
        IntegerType::Type{I} = Int,
        EdgeType::Type{E} = NTuple{2, IntegerType},
        TriangleType::Type{V} = NTuple{3, IntegerType},
        EdgesType::Type{Es} = Set{EdgeType},
        TrianglesType::Type{Ts} = Set{TriangleType},
    ) where {I, E, V, Es, Ts}
    return _triangulate_rectangle(a, b, c, d, nx, ny, I, E, V, Es, Ts, single_boundary, delete_ghosts)
end
@inline function _triangulate_rectangle(
        a, b, c, d, nx, ny,
        ::Type{I}, ::Type{E}, ::Type{V}, ::Type{Es}, ::Type{Ts},
        single_boundary, delete_ghosts,
    ) where {I, E, V, Es, Ts}
    T, sub2ind = get_lattice_triangles(nx, ny, Ts, V)
    points = get_lattice_points(a, b, c, d, nx, ny, sub2ind)
    boundary_nodes = get_lattice_boundary(nx, ny, sub2ind, Val(single_boundary), I)
    tri = Triangulation(points, T, boundary_nodes; IntegerType = I, EdgeType = E, TriangleType = V, EdgesType = Es, TrianglesType = Ts, delete_ghosts)
    compute_representative_points!(tri)
    return tri
end


"""
    get_lattice_triangles(nx, ny, Ts, V) 

Computes the triangles defining a lattice with `nx` and `ny` points in the `x`- and `y`-directions,
respectively.

See [`triangulate_rectangle`](@ref).

# Arguments 
- `nx`: The number of `x` points in the lattice.
- `ny`: The number of `y` points in the lattice.
- `Ts`: The type to use for representing a collection of triangles.
- `V`: The type to use for representing an individual triangle.

# Outputs 
- `T`: The collection of triangles.
- `sub2ind`: A map that takes cartesian indices `(i, j)` into the associated linear index along the lattice. See `LinearIndices`.
"""
@inline function get_lattice_triangles(nx, ny, ::Type{Ts}, ::Type{V}) where {Ts, V}
    T = Ts()
    sub2ind = LinearIndices((1:nx, 1:ny))
    for j in 1:(ny - 1)
        for i in 1:(nx - 1)
            u = sub2ind[CartesianIndex(i, j)]
            v = sub2ind[CartesianIndex(i + 1, j)]
            w = sub2ind[CartesianIndex(i, j + 1)]
            τ = construct_triangle(V, u, v, w)
            add_triangle!(T, τ)
            u = sub2ind[CartesianIndex(i, j + 1)]
            v = sub2ind[CartesianIndex(i + 1, j)]
            w = sub2ind[CartesianIndex(i + 1, j + 1)]
            τ = construct_triangle(V, u, v, w)
            add_triangle!(T, τ)
        end
    end
    return T, sub2ind
end


"""
    get_lattice_points(a, b, c, d, nx, ny, sub2ind)

Returns the points on a lattice `[a, b] × [c, d]`.

See [`triangulate_rectangle`](@ref).

# Arguments 
- `a`: The minimum `x`-coordinate. 
- `b`: The maximum `x`-coordinate.
- `c`: The minimum `y`-coordinate.
- `d`: The maximum `y`-coordinate.
- `nx`: The number of points in the `x`-direction.
- `ny`: The number of points in the `y`-direction.
- `sub2ind`: The map returned from [`get_lattice_triangles`](@ref).

# Outputs 
- `points`: The points on the lattice, where `points[sub2ind[CartesianIndex(i, j)]]` is the point at the `i`th `x` point and the `j`th `y` point,
"""
@inline function get_lattice_points(a, b, c, d, nx, ny, sub2ind)
    points = Vector{NTuple{2, Float64}}(undef, nx * ny)
    Δx = (b - a) / (nx - 1)
    Δy = (d - c) / (ny - 1)
    for j in 1:ny
        y = c + (j - 1) * Δy
        for i in 1:nx
            x = a + (i - 1) * Δx
            idx = sub2ind[CartesianIndex(i, j)]
            points[idx] = (x, y)
        end
    end
    return points
end


"""
    get_lattice_boundary(nx, ny, sub2ind, single_boundary, IntegerType)

Returns the boundary nodes defining a lattice. 

See [`triangulate_rectangle`](@ref).

# Arguments 
- `nx`: The number of `x` points in the lattice.
- `ny`: The number of `y` points in the lattice.
- `sub2ind`: The map returned from [`get_lattice_triangles`](@ref).
- `single_boundary`: If `true`, then the boundary nodes are stored as a contiguous section. Otherwise, the boundary is split into four sections, in the order 
   bottom, right, top, left.
- `IntegerType`: The type used for representing vertices. 

# Outputs 
- `boundary_nodes`: The boundary nodes, returned according to `single_boundary` as described above.
"""
@inline function get_lattice_boundary(nx, ny, sub2ind, single_boundary::Val{B}, ::Type{I}) where {B, I}
    b1 = Vector{I}(undef, nx)
    b2 = Vector{I}(undef, ny)
    b3 = Vector{I}(undef, nx)
    b4 = Vector{I}(undef, ny)
    for i in 1:nx
        b1[i] = sub2ind[CartesianIndex(i, 1)]
    end
    for j in 1:ny
        b2[j] = sub2ind[CartesianIndex(nx, j)]
    end
    for i in nx:-1:1
        b3[nx - i + 1] = sub2ind[CartesianIndex(i, ny)]
    end
    for j in ny:-1:1
        b4[ny - j + 1] = sub2ind[CartesianIndex(1, j)]
    end
    if is_true(single_boundary)
        popfirst!(b2)
        popfirst!(b3)
        popfirst!(b4)
        boundary_nodes = vcat(b1, b2, b3, b4)
    else
        boundary_nodes = [b1, b2, b3, b4]
    end
    return boundary_nodes
end