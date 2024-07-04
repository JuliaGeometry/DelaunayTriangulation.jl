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
function triangulate_rectangle(a, b, c, d, nx, ny;
    single_boundary=false,
    delete_ghosts=false,
    IntegerType::Type{I}=Int,
    EdgeType::Type{E}=NTuple{2,IntegerType},
    TriangleType::Type{V}=NTuple{3,IntegerType},
    EdgesType::Type{Es}=Set{EdgeType},
    TrianglesType::Type{Ts}=Set{TriangleType}) where {I,E,V,Es,Ts}

    ## Define the triangles
    T = Ts()
    sub2ind = LinearIndices((1:nx, 1:ny))
    idx = 1
    for j in 1:(ny-1)
        for i in 1:(nx-1)
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

    ## Define the points 
    points = Vector{NTuple{2,Float64}}(undef, nx * ny)
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

    ## Define the boundary nodes 
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
        b3[nx-i+1] = sub2ind[CartesianIndex(i, ny)]
    end
    for j in ny:-1:1
        b4[ny-j+1] = sub2ind[CartesianIndex(1, j)]
    end
    if single_boundary
        popfirst!(b2)
        popfirst!(b3)
        popfirst!(b4)
        boundary_nodes = vcat(b1, b2, b3, b4)
    else
        boundary_nodes = [b1, b2, b3, b4]
    end

    ## Now triangulate 
    tri = Triangulation(points, T, boundary_nodes; IntegerType, EdgeType, TriangleType,EdgesType, TrianglesType, delete_ghosts)
    compute_representative_points!(tri)
    return tri
end
