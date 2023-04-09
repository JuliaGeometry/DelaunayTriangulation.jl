"""
    triangulate_rectangle(a, b, c, d, nx, ny;
        single_boundary=false,
        add_ghost_triangles=true,
        IntegerType::Type{I}=Int64,
        EdgeType::Type{E}=NTuple{2,IntegerType},
        TriangleType::Type{V}=NTuple{3,IntegerType},
        EdgesType::Type{Es}=Set{EdgeType},
        TrianglesType::Type{Ts}=Set{TriangleType}) where {I,E,V,Es,Ts}

Computes a triangulation of the rectangular grid `[a, b] × [c, d]` with points `(xᵢ, yⱼ)`, 
`i = 1, …, nx`, `j = 1, …, ny`, where `xᵢ = a + (i-1)(b-a)/(nx-1)` and `yⱼ = b + (j-1)(d-c)/(ny-1)`. If 
the boundary of the rectangle should be considered as one single boundary, use `single_boundary = false`, and if 
you want the four sides of the boundary to be separated use `single_boundary = true`. 

Returns a [`Triangulation`](@ref) representing the triangulation.
"""
function triangulate_rectangle(a, b, c, d, nx, ny;
    single_boundary=false,
    add_ghost_triangles=true,
    IntegerType::Type{I}=Int64,
    EdgeType::Type{E}=NTuple{2,IntegerType},
    TriangleType::Type{V}=NTuple{3,IntegerType},
    EdgesType::Type{Es}=Set{EdgeType},
    TrianglesType::Type{Ts}=Set{TriangleType}) where {I,E,V,Es,Ts}

    ## Define the triangles
    T = initialise_triangles(Ts)
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
        boundary_nodes = [b1...,
            b2[(begin+1):end]...,
            b3[(begin+1):end]...,
            b4[(begin+1):end]...]
    else
        boundary_nodes = [b1, b2, b3, b4]
    end

    ## Now triangulate 
    tri = Triangulation(points, T, boundary_nodes; IntegerType, EdgeType, TriangleType,
        EdgesType, TrianglesType, add_ghost_triangles)
    compute_representative_points!(tri)
    return tri
end
