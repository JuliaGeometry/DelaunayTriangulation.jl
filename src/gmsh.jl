"""
    generate_mesh(x::Vector{Vector{Float64}}, y::Vector{Vector{Float64}}, ref; <keyword arguments>)
    generate_mesh(x::Vector{LinRange{Float64, Int64}}, y::Vector{LinRange{Float64, Int64}}, ref; <keyword arguments>)
    generate_mesh(x::Vector{Float64}, y::Vector{Float64}, ref; <keyword arguments>)

Generates the mesh of a domain when given its boundary, specified by `(x, y)` in counter-clockwise order, 
using `gmsh_path` to locate gmsh.

# Arguments 
- `x`: The `x`-coordinates of the boundary, given in counter-clockwise order.
- `y`: The `y`-coordinates of the boundary, given in counter-clockwise order.
- `ref`: The mesh refinement parameter.

If the boundary is to be regarded as a single type (e.g. a single boundary condition for the whole 
boundary), then the `x` and `y` can simply be vectors of points. Otherwise, if the boundary should be 
regarded as the union of multiple segment types, then the `x` and `y` should be vectors of vectors, 
with each entry being the vector of points for the corresponding boundary type. In this latter case, 
the points are still given in counter-clockwise order, with the additional requirement that the 
endpoints of each individual segment in `(x, y)` are connected, i.e. `x[i][end] = x[i+1][1]` and 
similarly for `y` (and, of course, `x[end][end] = x[1][1]` and similarly for `y`).

# Keyword Arguments 
- `mesh_algorithm = 6`: The meshing algorithm to use; see the gmsh documentation. Currently defaults to Frontal-Delaunay.
- `gmsh_path = "./gmsh-4.9.4-Windows64/gmsh.exe"`: A directory giving the location of the gmsh executable. 
- `verbosity = 0`: Verbosity level for the gmsh output; see the gmsh documentation.

# Outputs 
- `mesh`: The [`DelaunayTriangulation`](@ref) of the domain.
- `BN`: The indices of the boundary nodes.
"""
function generate_mesh end
@noinline function generate_mesh(x::Vector{Vector{Float64}}, y::Vector{Vector{Float64}}, ref;
    mesh_algorithm=6, verbosity=0, gmsh_path="./gmsh-4.9.4-Windows64/gmsh.exe",
    num_threads=Base.Threads.nthreads())
    ## Mesh generation
    all_x = reduce(vcat, x)
    all_y = reduce(vcat, y)
    Nbpts = length(all_x)
    keys = 1:length(x)
    shift_keys = keys
    open("meshgeometry.geo", "w") do fout
        # Define parameters
        #write(fout, "SetFactory(\"OpenCASCADE\");\n")
        write(fout, "r = $ref;\n")
        write(fout, "Mesh.Algorithm = $mesh_algorithm;\n")
        write(fout, "Mesh.Format = 1;\n")
        write(fout, "General.Verbosity = $verbosity;\n") #https://gitlab.onelab.info/gmsh/gmsh/-/issues/1133                                      
        # Start by defining the points 
        for k in 1:Nbpts
            write(fout, @sprintf("Point(%i) = {%2.15f,%2.15f,0,r};\n", k, all_x[k], all_y[k]))
        end
        left_starts = [1, cumsum(length.(x[1:end-1])) .+ 1...]
        right_starts = [cumsum(length.(x)) .- 1...]
        # Fit splines through the individual sections 
        for key in keys
            write(fout, "BSpline($key) = {")
            for k in left_starts[key]:(right_starts[key])
                write(fout, "$k,")
            end
            write(fout, "$(left_starts[(key%length(keys))+1])};\n")
        end
        # Now define the plane surfaces
        write(fout, "Line Loop(2881000) = {")
        for key in keys[1:end-1]
            write(fout, "$key, ")
        end
        write(fout, "$(keys[end])};\n")
        write(fout, "Plane Surface($(2length(keys)+1)) = {2881000};\n")
        # Now define the physical properties 
        for key in keys
            write(fout, "Physical Line($(key)) = {$key};\n")
        end
        write(fout, "Physical Surface($(5000 + length(keys))) = {$(2length(keys)+1)};\n")
        #write(fout, "Coherence;")
    end
    command = [gmsh_path; "meshgeometry.geo"; "-2"; "-format"; "msh2"; "-nt"; "$num_threads"]
    run(`$command`)
    ## Read the mesh
    elements = ElasticArray{Int64}(undef, 3, 0)
    boundary_nodes = [Vector{Int64}([]) for _ in 1:length(x)]
    nodes = ElasticArray{Float64}(undef, 2, 0)
    open("meshgeometry.msh", "r") do fid
        while true
            tline = readline(fid) # Start reading the file 
            if !(tline isa String) # Make sure the file isn't broken
                break
            end
            if tline == "\$MeshFormat" # Process the mesh format
                tline = readline(fid)
                if !(tline isa String)
                    break
                end
                MF = @scanf(tline, "%g%g%g", Float64, Int64, Int64)[2:end] # Read in the mesh format
                if MF ≠ (2.2, 0, 8)
                    error("Wrong mesh format selected.")
                end
                tline = readline(fid)
                if !(tline isa String) || !(tline == "\$EndMeshFormat")
                    break
                end
            elseif tline == "\$Nodes" # Process the nodes
                tline = readline(fid)
                no_nodes = @scanf(tline, "%g", Int64)[2] # Read the number of nodes
                append!(nodes, zeros(2, no_nodes))
                idx = 1
                while idx ≤ no_nodes
                    tline = readline(fid)
                    type = @scanf(tline, "%g%g%g%g", Int64, Float64, Float64, Float64)
                    nodes[:, idx] .= type[3:4]
                    idx += 1
                end
                tline = readline(fid)
                if !(tline isa String) || tline ≠ "\$EndNodes"
                    break
                end
            elseif tline == "\$Elements" # Process the elements 
                tline = readline(fid)
                if !(tline isa String)
                    break
                end
                no_elements = 0
                no_objects = @scanf(tline, "%g", Int64)[2]
                for _ = 1:no_objects
                    tline = readline(fid)
                    type = @scanf(tline, "%g%g%g%g%g%g%g%g", Int64, Int64, Int64, Int64, Int64, Int64, Int64, Int64)
                    if type[3] == 2
                        no_elements += 1
                        append!(elements, type[7:9])
                    end
                    if type[5] ∈ shift_keys
                        append!(boundary_nodes[type[5]], type[7:8])
                    end
                end
                tline = readline(fid)
                if !(tline isa String) || tline ≠ "\$EndElements"
                    break
                end
            else
                println("End of file or unknown type occurred.")
                break
            end
        end
    end
    T, adj, adj2v, DG = triangulate(elements, nodes, boundary_nodes)
    return T, adj, adj2v, DG, nodes, unique.(boundary_nodes)
end
@inline function generate_mesh(x::Vector{LinRange{Float64,Int64}}, y::Vector{LinRange{Float64,Int64}}, ref; mesh_algorithm=6, verbosity=0, gmsh_path="./gmsh-4.9.4-Windows64/gmsh.exe", num_threads=Base.Threads.nthreads())
    return generate_mesh(collect.(x), collect.(y), ref; verbosity, mesh_algorithm, gmsh_path, num_threads)
end
@inline function generate_mesh(x::Vector{Float64}, y::Vector{Float64}, ref; verbosity=0, mesh_algorithm=6, gmsh_path="./gmsh-4.9.4-Windows64/gmsh.exe", num_threads=Base.Threads.nthreads())
    return generate_mesh([x], [y], ref; verbosity, mesh_algorithm, gmsh_path, num_threads)
end

function triangulate(triangles,
    nodes::AbstractMatrix,
    boundary_nodes;
    IntegerType::Type{I}=Int64,
    EdgeType::Type{E}=NTuple{2,IntegerType},
    TriangleType::Type{V}=NTuple{3,IntegerType},
    EdgesType::Type{Es}=Set{EdgeType},
    TrianglesType::Type{Ts}=Set{TriangleType},
    sort_boundary=true) where {I,E,V,Es,Ts}
    boundary_nodes = deepcopy(boundary_nodes)
    boundary_nodes = reduce(vcat, boundary_nodes)
    sort_boundary && unique!(boundary_nodes)
    T = Ts()
    adj = Adjacent{I,E}()
    adj2v = Adjacent2Vertex{I,Es,E}()
    DG = DelaunayGraph{I}()
    itr = triangles isa AbstractMatrix ? eachcol(triangles) : triangles
    for (i, j, k) in itr
        add_triangle!(i, j, k, T, adj, adj2v, DG; protect_boundary=true)
    end
    for i in boundary_nodes
        add_neighbour!(DG, I(BoundaryIndex), i)
    end
    if sort_boundary
        cx, cy = sum(nodes; dims=2) / size(nodes, 2)
        θ = zeros(length(boundary_nodes))
        for (j, i) in pairs(boundary_nodes)
            p = @view nodes[:, i]
            x, y = p
            θ[j] = atan(y - cy, x - cx)
        end
        sort_idx = sortperm(θ)
        permute!(boundary_nodes, sort_idx)
        reverse!(boundary_nodes) # cw 
    end
    n = length(boundary_nodes)
    push!(boundary_nodes, boundary_nodes[begin])
    for i in 1:n
        u = boundary_nodes[i]
        v = boundary_nodes[i+1]
        add_edge!(adj2v, I(BoundaryIndex), u, v)
        add_edge!(adj, u, v, I(BoundaryIndex))
    end
    return T, adj, adj2v, DG
end
