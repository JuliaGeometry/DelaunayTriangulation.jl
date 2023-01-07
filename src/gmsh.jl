const LINE_LOOP_IDENTIFIER = 2881000
const PLANE_SURFACE_IDENTIFIER = 2LINE_LOOP_IDENTIFIER + 1
const PHYSICAL_SURFACE_IDENTIFIER = 4PLANE_SURFACE_IDENTIFIER + 38
const MESH_FORMAT = (2.2, 0, 8)

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
- `mesh`: The [`Triangulation`](@ref) of the domain.
- `BN`: The indices of the boundary nodes.
"""
function generate_mesh end
@noinline function generate_mesh(x::Vector{Vector{Float64}}, y::Vector{Vector{Float64}}, ref;
    mesh_algorithm=6, verbosity=0, gmsh_path="./gmsh-4.9.4-Windows64/gmsh.exe",
    num_threads=Base.Threads.nthreads())
    segment_numbers = write_gmsh(x, y, ref; mesh_algorithm, verbosity)
    run_gmsh(gmsh_path, num_threads)
    elements, nodes, boundary_nodes = read_gmsh(segment_numbers)
    T, adj, adj2v, DG = triangulate(elements, nodes, boundary_nodes)
    return Triangulation(T, adj, adj2v, DG, nodes), unique.(boundary_nodes)
end
@inline function generate_mesh(x::Vector{LinRange{Float64,Int64}}, y::Vector{LinRange{Float64,Int64}}, ref; mesh_algorithm=6, verbosity=0, gmsh_path="./gmsh-4.9.4-Windows64/gmsh.exe", num_threads=Base.Threads.nthreads())
    return generate_mesh(collect.(x), collect.(y), ref; verbosity, mesh_algorithm, gmsh_path, num_threads)
end
@inline function generate_mesh(x::Vector{Float64}, y::Vector{Float64}, ref; verbosity=0, mesh_algorithm=6, gmsh_path="./gmsh-4.9.4-Windows64/gmsh.exe", num_threads=Base.Threads.nthreads())
    return generate_mesh([x], [y], ref; verbosity, mesh_algorithm, gmsh_path, num_threads)
end

"""
    generate_mesh(a, b, c, d, ref;
        mesh_algorithm=6,
        verbosity=0,
        gmsh_path="./gmsh-4.9.4-Windows64/gmsh.exe",
        num_threads=Base.Threads.nthreads(),
        single_boundary=true)

Generates a mesh of the rectangle `[a, b] × [c, d]` using Gmsh.

# Arguments 
- `a, b, c, d`: These define the rectangle `[a, b] × [c, d]`.
- `ref`: The refinement parameter. 

# Keyword Arguments 
- `mesh_algorithm = 6`: The meshing algorithm to use; see the gmsh documentation. Currently defaults to Frontal-Delaunay.
- `gmsh_path = "./gmsh-4.9.4-Windows64/gmsh.exe"`: A directory giving the location of the gmsh executable. 
- `verbosity = 0`: Verbosity level for the gmsh output; see the gmsh documentation.
- `single_boundary=true`: If `true`, then the boundary is taken as a single line loop, whereas if `false` then boundary nodes are given for the individual segments, with the bottom, right, top, and left segments being the first, second, third, and fourth segments, respectively.

# Outputs 
- `mesh`: The [`Triangulation`](@ref) of the rectangle.
- `BN`: The indices of the boundary nodes.
"""
function generate_mesh(a, b, c, d, ref;
    mesh_algorithm=6,
    verbosity=0,
    gmsh_path="./gmsh-4.9.4-Windows64/gmsh.exe",
    num_threads=Base.Threads.nthreads(),
    single_boundary=true)
    write_gmsh(a, b, c, d, ref; mesh_algorithm,verbosity, single_boundary)
    run_gmsh(gmsh_path, num_threads)
    elements, nodes, boundary_nodes = read_gmsh(single_boundary ? [1] : 1:4)
    T, adj, adj2v, DG = triangulate(elements, nodes, boundary_nodes)
    return Triangulation(T, adj, adj2v, DG, nodes), unique.(boundary_nodes)
end

"""
    write_gmsh(x, y, ref; mesh_algorithm=6, verbosity=0)

Creates the `.geo` file for the mesh.
"""
function write_gmsh(x, y, ref; mesh_algorithm=6, verbosity=0)
    ## Mesh generation
    all_x = reduce(vcat, x)
    all_y = reduce(vcat, y)
    Nbpts = length(all_x)
    keys = 1:length(x)
    segment_numbers = keys
    open("meshgeometry.geo", "w") do fout
        # Define parameters
        _write_mesh_settings!(fout, ref, mesh_algorithm, verbosity)
        # Start by defining the points 
        for k in 1:Nbpts
            _write_point!(fout, k, all_x[k], all_y[k])
        end
        left_starts = [1, cumsum(length.(x[1:end-1])) .+ 1...] # left_starts[i] is the starting node of the ith segment
        right_starts = [cumsum(length.(x)) .- 1...]
        # Define the geometry
        for key in keys
            _write_spline!(fout, key, left_starts[key]:right_starts[key], left_starts[(key%length(keys))+1])
        end
        _write_line_loop!(fout, keys, LINE_LOOP_IDENTIFIER)
        _write_plane_surface!(fout, LINE_LOOP_IDENTIFIER, PLANE_SURFACE_IDENTIFIER)
        for key in keys
            _write_physical_line!(fout, key)
        end
        _write_physical_surface!(fout, PLANE_SURFACE_IDENTIFIER, PHYSICAL_SURFACE_IDENTIFIER)
    end
    return segment_numbers
end

"""
    write_gmsh(a, b, c, d, ref; mesh_algorithm=6, verbosity=0, single_boundary=true)

Creates the `.geo` file for the mesh of the rectangle `[a, b] × [c, d]`.
"""
function write_gmsh(a, b, c, d, ref; mesh_algorithm=6, verbosity=0, single_boundary=true)
    open("meshgeometry.geo", "w") do fout
        _write_mesh_settings!(fout, ref, mesh_algorithm, verbosity)
        _write_point!(fout, 1, a, c)
        _write_point!(fout, 2, b, c)
        _write_point!(fout, 3, b, d)
        _write_point!(fout, 4, a, d)
        _write_line!(fout, 1, 1, 2)
        _write_line!(fout, 2, 2, 3)
        _write_line!(fout, 3, 3, 4)
        _write_line!(fout, 4, 4, 1)
        _write_line_loop!(fout, 1:4, LINE_LOOP_IDENTIFIER)
        _write_plane_surface!(fout, LINE_LOOP_IDENTIFIER, PLANE_SURFACE_IDENTIFIER)
        if single_boundary
            _write_physical_line!(fout, 1:4, 1)
        else
            _write_physical_line!(fout, 1)
            _write_physical_line!(fout, 2)
            _write_physical_line!(fout, 3)
            _write_physical_line!(fout, 4)
        end
        _write_physical_surface!(fout, PLANE_SURFACE_IDENTIFIER, PHYSICAL_SURFACE_IDENTIFIER)
    end
    return nothing
end

function run_gmsh(gmsh_path, num_threads)
    command = [gmsh_path; "meshgeometry.geo"; "-2"; "-format"; "msh2"; "-nt"; "$num_threads"]
    run(`$command`)
    return nothing
end

"""
    read_gmsh(segment_numbers)

Reads the `.msh` file for the mesh.
"""
function read_gmsh(segment_numbers)
    elements = ElasticArray{Int64}(undef, 3, 0)
    boundary_nodes = [Vector{Int64}([]) for _ in segment_numbers]
    nodes = ElasticArray{Float64}(undef, 2, 0)
    open("meshgeometry.msh", "r") do fid
        while true
            _read_line!(fid, elements, nodes, boundary_nodes, segment_numbers) && break
        end
    end
    return elements, nodes, boundary_nodes
end

function _read_line!(fid, elements, nodes, boundary_nodes, segment_numbers)
    tline = readline(fid) # Start reading the file 
    _check_file(tline) && return true
    if tline == "\$MeshFormat" # Process the mesh format
        _read_mesh_format!(fid)
    elseif tline == "\$Nodes" # Process the nodes
        _read_nodes!(nodes, fid)
    elseif tline == "\$Elements" # Process the elements 
        _read_elements!(elements, boundary_nodes, fid, segment_numbers)
    else
        println("End of file or unknown type occurred.")
        return true
    end
    return false
end

function _write_mesh_settings!(fout, ref, mesh_algorithm, verbosity)
    write(fout, "r = $ref;\n")
    write(fout, "Mesh.Algorithm = $mesh_algorithm;\n")
    write(fout, "Mesh.Format = 1;\n")
    write(fout, "General.Verbosity = $verbosity;\n") #https://gitlab.onelab.info/gmsh/gmsh/-/issues/1133                                      
    return nothing
end

function _write_point!(fout, point_index, x, y)
    write(fout, @sprintf("Point(%i) = {%2.15f, %2.15f, 0, r};\n", point_index, x, y))
    return nothing
end

function _write_spline!(fout, spline_index, point_indices, last_index)
    write(fout, "BSpline($spline_index) = {")
    for k in point_indices
        write(fout, "$k,")
    end
    write(fout, "$last_index};\n") # need to connect with the next segment - this is also separate so that there's no comma
    return nothing
end

function _write_line!(fout, line_index, first_point, second_point)
    write(fout, "Line($line_index) = {$first_point, $second_point};\n")
    return nothing
end

function _write_line_loop!(fout, line_indices, identifier=LINE_LOOP_IDENTIFIER)
    write(fout, "Line Loop($identifier) = {")
    for line in @views line_indices[1:end-1]
        write(fout, "$line, ")
    end
    write(fout, "$(line_indices[end])};\n")
    return nothing
end

function _write_plane_surface!(fout, line_loop, identifier=PLANE_SURFACE_IDENTIFIER)
    write(fout, "Plane Surface($identifier) = {$line_loop};\n")
    return nothing
end

function _write_physical_line!(fout, line_index, identifier=line_index)
    write(fout, "Physical Line($identifier) = {$line_index};\n")
    return nothing
end
function _write_physical_line!(fout, line_indices::AbstractVector, identifier=PHYSICAL_LINE_IDENTIFIER)
    write(fout, "Physical Line($identifier) = {")
    for k in @views line_indices[begin:end-1]
        write(fout, "$k,")
    end
    write(fout, "$(line_indices[end])};\n")
    return nothing
end

function _write_physical_surface!(fout, plane_surface=PLANE_SURFACE_IDENTIFIER, identifier=PHYSICAL_SURFACE_IDENTIFIER)
    write(fout, "Physical Surface($identifier) = {$plane_surface};\n")
    return nothing
end

function _check_file(tline)
    return !(tline isa String)
end

function _read_mesh_format!(fid)
    tline = readline(fid)
    if !(tline isa String)
        return nothing
    end
    MF = @scanf(tline, "%g%g%g", Float64, Int64, Int64)[2:end] # Read in the mesh format
    if MF ≠ MESH_FORMAT
        error("Wrong mesh format, $MF, selected.")
    end
    tline = readline(fid)
    if !(tline isa String) || !(tline == "\$EndMeshFormat")
        return nothing
    end
    return nothing
end

function _read_nodes!(nodes, fid)
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
        return nothing
    end
    return nothing
end

function _read_elements!(elements, boundary_nodes, fid, segment_numbers)
    tline = readline(fid)
    if !(tline isa String)
        return nothing
    end
    no_elements = 0
    no_objects = @scanf(tline, "%g", Int64)[2]
    for _ in 1:no_objects
        tline = readline(fid)
        type = @scanf(tline, "%g%g%g%g%g%g%g%g", Int64, Int64, Int64, Int64, Int64, Int64, Int64, Int64)
        if type[3] == 2
            no_elements += 1
            @views append!(elements, type[7:9])
        end
        if type[5] ∈ segment_numbers
            @views append!(boundary_nodes[type[5]], type[7:8])
        end
    end
    tline = readline(fid)
    if !(tline isa String) || tline ≠ "\$EndElements"
        return nothing
    end
    return nothing
end