"""
    generate_mesh(x, y, ref;
        mesh_algorithm=6,
        gmsh_path="./gmsh-4.11.1-Windows64/gmsh.exe",
        verbosity=0,
        convert_result=true,
        add_ghost_triangles=false,
        check_args=true)

Using Gmsh, generates a mesh of the domain defined by `(x, y)`. 

# Arguments
- `x, y`

These are the coordinates defining the curves that define the boundaries of the domain. All 
curves are to be positively oriented, meaning the outermost boundary should be counter-clockwise 
while the interior boundaries should be clockwise.

There are three accepted forms for `x, y`:

-- `Vector{Vector{Vector{Float64}}}`: Multiple holes 

This form for `x` and `y` will separate the generated domain into separate curves, 
where `(x[1], y[1])` is the outer-most boundary and the remaining parts of `x`
and `y` define holes inside the domain. In general, `(x[m][n], y[m][n])` is the 
`n`th segment of the `m`th curve. It is assumed that `x[m][n][end] == x[m][n+1][begin]`,
and that `x[m][end][end] == x[m][begin][1]`, i.e. separate segments have shared endpoints.
Similar conditions hold for `y`. In this case, `(x[1], y[1])` should define a counter-clockwise 
curve while `(x[m], y[m])` should be a clockwise curve for `m > 1`.

-- `Vector{Vector{Float64}}`: One boundary, multiple segments 

This form for `x` and `y` uses just one boundary, but splits the boundary into multiple 
segments as described above. As in the case above, separate segments should have shared 
endpoints. In this case, the boundary should be provided counter-clockwise.

-- `Vector{Float64}`: One boundary, one segment 

This form for `x` and `y` uses just one boundary, and assumes the boundary is one 
continuous segment. It is assumed that `x[begin] == x[end]` and similarly for `y`. The boundary 
should be provided counter-clockwise.

- `ref`

This is the refinement parameter - smaller `ref` means more elements. 

# Keyword Arguments 
- `mesh_algorithm=6`

The algorithm to use for meshing. `6` means Frontal-Delaunay.
- `gmsh_path="./gmsh-4.11.1-Windows64/gmsh.exe"`

The location of the gmsh executable.
- `verbosity=0`

The verbosity level for Gmsh.
- `convert_result=true`

If `true`, the final result is converted into a [`Triangulation`](@ref) type. Otherwise, 
`(triangles, points, boundary_nodes)` is returned.

- `add_ghost_triangles=false`

If `convert_result`, then this declares whether or not ghost triangles should be added when 
converting the result into a [`Triangulation`](@ref) type. See also [`Triangulation`](@ref).

- `check_arguments=true`

Whether to verify the arguments have the correct orientation. 

# Outputs 
If `convert_result`, then the final result is a [`Triangulation`](@ref) type. Otherwise, 
the following values are returned:

- `elements`

The triangular elements for the mesh. 
- `nodes`

The nodes in the mesh. 
- `boundary_nodes`

The boundary nodes in the mesh. All boundaries are positively oriented relative to the interior, 
meaning the outermost boundary is counter-clockwise while the interior boundaries are clockwise.

The form of `boundary_nodes` matches the form of `x` and `y`, i.e. (also 
including matching endpoints):

-- `Vector{Vector{Vector{Float64}}}`

In this case, `boundary_nodes[m][n]` give the indices in `nodes` 
corresponding to the boundary specified by `(x[m][n], y[m][n])`. The nodes 
`boundary_nodes[1]` define a counter-clockwise curve, while `boundary_nodes[m]`
is clockwise for `m > 1`.

-- `Vector{Vector{Float64}}`

In this case, `boundary_nodes[n]` gives the indices in `nodes` 
corresponding to the boundary specified by `(x[n], y[n])`. The boundary nodes 
define a counter-clockwise curve in this case.

-- `Vector{Float64}`

In this case, `boundary_nodes` gives the indices in `nodes` 
corresponding to the boundary specified by `(x, y)`. The boundary nodes 
define a counter-clockwise curve in this case.

# Extended help 

The function proceeds in four steps:

1. Mesh generation

Here, we write a file "meshgeometry.geo" in the working directory. This file takes the form

    r = ref;
    Mesh.Algorithm = mesh_algorithm; 
    Mesh.Format = 1;
    General.Verbosity = 0;
    Point(<point index>) = {<x>, <y>, 0, r}; # For each point
    Line(<line index>) = {<initial point>, <final point>}; # For each line 
    Curve Loop(<boundary index>) = {<line 1>, <line 2>, ...}; # For each boundary 
    Plane Surface(1) = {<curve 1>, <curve 2>, ...}; # <curve 1> = 1 and is the outermost boundary, while <curve i> = i, i > 1, are boundaries of interior holes 
    Physical Curve(<last line index + i>) = {<line 1>, <line 2>, ...}; # For i ranging over the number of segments, and the lines represent that segment 
    Physical Surface(1) = {1};

Most importantly, every edge input into the function `generate_mesh` will be included in the mesh. An 
older version of this function previously used cubic splines for defining boundary curves, but this 
has the consequence that (1) not every edge put into the function is included, and (2) the boundary 
is not exactly represented. 

The function that handles this generation is `write_gmsh`.

2. Mesh writing 

The "meshgeometry.geo" file is then used to mesh the domain, running the terminal command 

    gmsh_path "meshgeometry.geo" -2 -format msh2 

This creates a file "meshgeometry.msh" in the same working directory. 

The function that handles this writing is `run_gmsh`.

3. Mesh reading 

Once "meshgeometry.msh" is created, we need to read it. The format used (\$MeshFormat) is 
2.2, but note that as of writing (13/01/2013), the most modern format is 4.1.

The "meshgeometry.msh" file is split into groups:

3a. MeshFormat 

This just reads off the format of the file used. This part of the file 
is read using `read_mesh_format!`.

3b. Nodes 

This lists the node indices and all the coordinates of the nodes, with the 
first line giving the number of nodes. A single line in this section, 
after the first, takes the form 

    <node index> <x> <y> 0

and we read this using `read_node_line`. This entire part of the file 
is read using `read_nodes!`.

3c. Elements 

The first line in this part of the file is the number of elements, though 
here elements refer to both the lines and the triangles. The lines (edges) 
are listed first, with each line taking the form 

    <line index> 1 2 <boundary index> <> <left node> <right node>

and will be in counter-clockwise order. After the lines are listed, all the
triangles follow, with each line in this part taking the form

    <triangle index> 2 2 1 1 <node 1> <node 2> <node 3>

with each triangle positively oriented. These lines are read using 
`read_element_line`. The entire part of the file is read using 
`read_elements!`.

4. Conversion to [`Triangulation`](@ref)

Once the file "meshgeometry.geo" has been read, we have a list of 
triangular elements, nodes, and boundary nodes. These need to all be 
converted into a [`Triangulation`](@ref) type, and a constructor of 
[`Triangulation`](@ref) is used to accomplish this.
"""
function generate_mesh(x::AAA, y::AAA, ref;
    mesh_algorithm=6,
    gmsh_path="./gmsh-4.11.1-Windows64/gmsh.exe",
    verbosity=0,
    convert_result=true,
    add_ghost_triangles=false,
    check_arguments=true) where {F<:
    Number,A<:AbstractVector{F},AA<:AbstractVector{A},AAA<:AbstractVector{AA}}
    identifiers, identifier_dict = write_gmsh(x, y, ref; mesh_algorithm, verbosity)
    run_gmsh(gmsh_path)
    elements, nodes, boundary_nodes = read_gmsh(identifiers, identifier_dict)
    check_arguments && check_args(nodes, boundary_nodes)
    if !convert_result
        return elements, nodes, boundary_nodes
    else
        tri = Triangulation(nodes, elements, boundary_nodes; add_ghost_triangles)
        compute_representative_points!(tri)
        return tri
    end
end

function generate_mesh(x::AA, y::AA, ref;
    mesh_algorithm=6,
    gmsh_path="./gmsh-4.11.1-Windows64/gmsh.exe",
    verbosity=0,
    convert_result=true,
    add_ghost_triangles=false,
    check_arguments=true) where {F<:
    Number,
    A<:
    AbstractVector{F},
    AA<:
    AbstractVector{A}}
    elements, nodes, boundary_nodes = generate_mesh([x], [y], ref; mesh_algorithm,
        gmsh_path, verbosity,
        convert_result=false,
        check_arguments)
    check_arguments && check_args(nodes, boundary_nodes[1])
    if !convert_result
        return elements, nodes, boundary_nodes[1]
    else
        tri = Triangulation(nodes, elements, boundary_nodes[1]; add_ghost_triangles)
        compute_representative_points!(tri)
        return tri
    end
end

function generate_mesh(x::A, y::A, ref;
    mesh_algorithm=6,
    gmsh_path="./gmsh-4.11.1-Windows64/gmsh.exe",
    verbosity=0,
    convert_result=true,
    add_ghost_triangles=false,
    check_arguments=true) where {F<:
    Number,
    A<:
    AbstractVector{F}}
    elements, nodes, boundary_nodes = generate_mesh([x], [y], ref; mesh_algorithm,
        gmsh_path, verbosity,
        convert_result=false,
        check_arguments)
    check_arguments && check_args(nodes, boundary_nodes[1])
    if !convert_result
        return elements, nodes, boundary_nodes[1]
    else
        tri = Triangulation(nodes, elements, boundary_nodes[1]; add_ghost_triangles)
        compute_representative_points!(tri)
        return tri
    end
end

"""
    generate_mesh(a, b, c, d, ref; 
        mesh_algorithm=6,
        gmsh_path="./gmsh-4.11.1-Windows64/gmsh.exe",
        verbosity=0,
        single_boundary=true,
        convert_result=true,
        add_ghost_triangles=false)

Generates a mesh of a rectangle `[a, b] × [c, d]`. Use `single_boundary=true` if 
each side of the rectangle should be treated the same, and `single_boundary=false` if 
you want boundary nodes for each side of the rectangle
     
See the main function [`generate_mesh`](@ref generate_mesh(::AAA, ::AAA, ::Any) where {F<:Number,A<:AbstractVector{F},AA<:AbstractVector{A},AAA<:AbstractVector{AA}}) for a description of the other 
arguments.
"""
function generate_mesh(a, b, c, d, ref;
    mesh_algorithm=6,
    gmsh_path="./gmsh-4.11.1-Windows64/gmsh.exe",
    verbosity=0,
    single_boundary=true,
    convert_result=true,
    add_ghost_triangles=false)
    if single_boundary
        x = [a, b, b, a, a]
        y = [c, c, d, d, c]
    else
        x = [[a, (a + b) / 2, b], [b, b, b], [b, (b + a) / 2, a], [a, a, a]]
        y = [[c, c, c], [c, (d + c) / 2, d], [d, d, d], [d, (d + c) / 2, c]]
    end
    tri = generate_mesh(x, y, ref; mesh_algorithm, gmsh_path, verbosity, convert_result,
        add_ghost_triangles)
    return tri
end

function write_gmsh(x, y, ref; mesh_algorithm, verbosity)
    num_lines = open("meshgeometry.geo", "w") do fout
        write_mesh_settings!(fout, ref, mesh_algorithm, verbosity)
        write_points!(fout, x, y)                                               # Defines all the points used (including duplicates)
        curve_ranges, segment_ranges, num_lines = write_curves!(fout, x, y)     # Defines all the lines, taking special care for connecting points 
        write_curve_loops!(fout, curve_ranges)                                  # This carefully connects all the segments, appropriately defining the separate curves as being connected
        write_plane_surface!(fout, length(curve_ranges))                        # Need to define a surface that connects all the regions, with the first boundary the outer-most boundary and all others holes
        write_physical_curves!(fout, segment_ranges, num_lines + 1)             # This attaches each segment with a specific physical curve, allowing us to identify it later and define the boundary nodes correctly
        write_physical_surface!(fout, 1, 1)                                     # Since we use physical properties above, Gmsh will only export features that are attached to something physical (unless Mesh.SaveAll = 1). Therefore, in order to get all the elements in the output, we assign a physical property to the whole domain
        return num_lines
    end
    identifiers = get_segment_identifiers(x) # Accessing the arrays
    num_segments = last(last(identifiers))
    shifted_identifiers = (1:num_segments) .+ num_lines # These are the identifiers in the Physical Curves above
    identifier_dict = get_segment_identifiers_dict(identifiers, shifted_identifiers) # Maps a Physical Curve to the indices (m, n) giving the place to put the boundary nodes 
    return identifiers, identifier_dict
end

function run_gmsh(gmsh_path)
    command = [gmsh_path; "meshgeometry.geo"; "-2"; "-format"; "msh2"]
    run(`$command`)
    return nothing
end

function read_gmsh(identifiers, identifier_dict)
    boundary_nodes = get_boundary_node_init(identifiers)
    elements = NTuple{3,Int64}[]
    nodes = NTuple{2,Float64}[]
    open("meshgeometry.msh", "r") do fid
        return read_mesh!(fid, elements, nodes, boundary_nodes, identifier_dict)
    end
    elements = Matrix(reinterpret(reshape, Int64, elements))
    nodes = Matrix(reinterpret(reshape, Float64, nodes))
    return elements, nodes, boundary_nodes
end

function get_segment_identifiers(x)
    num_segments = length.(x)
    key = 1
    keys = [zeros(Int64, i) for i in num_segments]
    for m in eachindex(x)
        for n in eachindex(x[m])
            keys[m][n] = key
            key += 1
        end
    end
    return keys
end

function get_segment_identifiers_dict(identifiers, shifted_identifiers)
    dict = Dict{Int64,NTuple{2,Int64}}()
    current_idx = firstindex(shifted_identifiers)
    for m in eachindex(identifiers)
        for n in eachindex(identifiers[m])
            dict[shifted_identifiers[current_idx]] = (m, n)
            current_idx += 1
        end
    end
    return dict
end

function get_boundary_node_init(identifiers)
    bn = [[Int64[] for _ in eachindex(identifiers[i])] for i in eachindex(identifiers)]
    for m in eachindex(identifiers)
        for n in eachindex(identifiers[m])
            sizehint!(bn[m][n], length(identifiers[m][n]))
        end
    end
    return bn
end

function write_mesh_settings!(fout, ref, mesh_algorithm, verbosity)
    write(fout, "r = $ref;\n")
    write(fout, "Mesh.Algorithm = $mesh_algorithm;\n")
    write(fout, "Mesh.Format = 1;\n")
    write(fout, "General.Verbosity = $verbosity;\n") #https://gitlab.onelab.info/gmsh/gmsh/-/issues/1133 
    return nothing
end

function write_point!(fout, idx, x, y)
    write(fout, "Point($idx) = {$x, $y, 0, r};\n")
    return nothing
end

function write_points!(fout, x::Vector{T}, y::Vector{T}, init=1) where {T<:Number}
    for i in eachindex(x, y)
        write_point!(fout, init, x[i], y[i])
        init += 1
    end
    return init
end
function write_points!(fout, x, y, init=1)
    for i in eachindex(x, y)
        init = write_points!(fout, x[i], y[i], init)
    end
    return init
end

function write_line!(fout, idx, i, j)
    write(fout, "Line($idx) = {$i, $j};\n")
    return nothing
end

function write_lines!(fout, x, y, init_point=1, init_line=1)
    for _ in firstindex(x):(lastindex(x)-1)
        write_line!(fout, init_line, init_point, init_point + 1)
        init_line += 1
        init_point += 1
    end
    return init_point, init_line
end

function write_segments!(fout, x, y, init_point=1, init_line=1)
    segment_line_idx = Vector{UnitRange{Int64}}(undef, length(x)) # This will store the range of line indices corresponding to each segment
    if length(x) == 1 # Only one segment
        orig_init_line = init_line
        orig_init_point = init_point
        seg_x = x[begin]
        seg_y = y[begin]
        # Only do begin:end-1, because it is assumed that seg_x[begin] == seg_x[end]
        init_point, init_line = @views write_lines!(fout, seg_x[begin:(end-1)],
            seg_y[begin:(end-1)], init_point,
            init_line)
        write_line!(fout, init_line, init_point, orig_init_point)
        init_point += 1
        init_line += 1
        segment_line_idx[1] = orig_init_line:(init_line-1)
        return init_point, init_line, segment_line_idx
    else
        orig_init_point_curve = init_point
        # Loop over each segment
        for (j, i) in enumerate(eachindex(x))
            seg_x = x[i]
            seg_y = y[i]
            orig_init_line = init_line
            if i == firstindex(x)
                # If we're at the start of the curve, we just go over each edge of the segment
                init_point, init_line = write_lines!(fout, seg_x, seg_y, init_point,
                    init_line)
                init_point += 1 # Need to move on to the next segment - skip the same initial point from the next segment and go on
            elseif firstindex(x) < i < lastindex(x)
                # We are on an interior segment. In this case, the initial point needs to come from 
                # the terinal point of the previous segment, but we still go all the way to the end.
                write_line!(fout, init_line, init_point - 1, init_point + 1)
                init_line += 1
                init_point += 1
                # Now having added the first point, go over the remainder of the segment 
                init_point, init_line = @views write_lines!(fout, seg_x[(begin+1):end],
                    seg_y[(begin+1):end],
                    init_point, init_line)
                init_point += 1
            elseif i == lastindex(x)
                # If we are on the very last segment (note that this is not the only segment, since we handle 
                # length(x) == 1 above), then we need to (1) handle the initial point as in the 
                # firstindex(x) < i < lastindex(x) case above, but also (2) the last point needs to be 
                # the first point on the curve
                # First, add the first point correctly
                write_line!(fout, init_line, init_point - 1, init_point + 1)
                init_line += 1
                init_point += 1
                # Add the points between the initial and terminal points
                init_point, init_line = @views write_lines!(fout,
                    seg_x[(begin+1):(end-1)],
                    seg_y[(begin+1):(end-1)],
                    init_point, init_line)
                # Now connect the last point with the initial point on the curve
                write_line!(fout, init_line, init_point, orig_init_point_curve)
                init_line += 1
                init_point += 1
            end
            segment_line_idx[j] = orig_init_line:(init_line-1) # -1 since we increment by 1 even at the end in the above functions
        end
    end
    return init_point, init_line, segment_line_idx
end

function write_curves!(fout, x, y, init_point=1, init_line=1)
    curve_line_idx = Vector{UnitRange{Int64}}(undef, length(x))
    segment_line_idx = Vector{Vector{UnitRange{Int64}}}(undef, length(x))
    for (j, i) in enumerate(eachindex(x))
        orig_init_line = init_line
        init_point, init_line, segment_line_idx[j] = write_segments!(fout, x[i], y[i],
            init_point, init_line)
        init_point += 1 # Need to +1 so that we skip the next point, which is the point that connects back with the initial point of the current segment i.e. not the next curve
        curve_line_idx[j] = orig_init_line:(init_line-1) # Need to -1 since write_segments! adds an extra +1 at the end
    end
    num_lines = init_line - 1
    return curve_line_idx, segment_line_idx, num_lines
end

function write_curve_loop!(fout, curve_idx, line_idxs)
    write(fout, "Curve Loop($curve_idx) = {")
    for i in firstindex(line_idxs):(lastindex(line_idxs)-1)
        write(fout, "$(line_idxs[i]), ")
    end
    write(fout, "$(line_idxs[end])};\n")
    return nothing
end

function write_curve_loops!(fout, curve_ranges)
    for i in eachindex(curve_ranges)
        write_curve_loop!(fout, i, curve_ranges[i])
    end
    return nothing
end

function write_plane_surface!(fout, num_curves)
    write(fout, "Plane Surface(1) = {")
    for i in 1:(num_curves-1)
        write(fout, "$i, ")
    end
    write(fout, "$num_curves};\n")
    return nothing
end

function write_physical_curve!(fout, id, segment_range)
    write(fout, "Physical Curve($id) = {")
    for i in firstindex(segment_range):(lastindex(segment_range)-1)
        write(fout, "$(segment_range[i]), ")
    end
    write(fout, "$(segment_range[end])};\n")
    return nothing
end

function write_physical_curves!(fout, segment_ranges::Vector{UnitRange{Int64}}, init=1)
    for segment_range in segment_ranges
        write_physical_curve!(fout, init, segment_range)
        init += 1
    end
    return init
end
function write_physical_curves!(fout, segment_ranges::Vector{Vector{UnitRange{Int64}}},
    init=1)
    for ranges in segment_ranges
        init = write_physical_curves!(fout, ranges, init)
    end
    return nothing
end

function write_physical_surface!(fout, surface_idx, identifier)
    write(fout, "Physical Surface($identifier) = {$surface_idx};\n")
    return nothing
end

function read_mesh!(fid, elements, nodes, boundary_nodes, identifier_dict)
    while true
        tline = readline(fid) # Start reading the file 
        is_not_string(tline) && continue
        if tline == "\$MeshFormat" # Process the mesh format
            read_mesh_format!(fid)
            continue
        elseif tline == "\$Nodes" # Process the nodes
            nodes = read_nodes!(nodes, fid)
            continue
        elseif tline == "\$Elements" # Process the elements 
            elements = read_elements!(elements, boundary_nodes, fid, identifier_dict)
            continue
        else
            println("End of file or unknown type occurred.")
            break
        end
        break
    end
    return nothing
end

is_not_string(tline::AbstractString) = false
is_not_string(tline) = true

function read_mesh_format(tline)
    all_vals = split(tline, ' ')
    v1 = parse(Float64, all_vals[1])
    v2 = parse(Int64, all_vals[2])
    v3 = parse(Int64, all_vals[3])
    return v1, v2, v3
end

function read_mesh_format!(fid)
    tline = readline(fid)
    is_not_string(tline) && return nothing
    MF = read_mesh_format(tline)
    if MF ≠ (2.2, 0, 8)
        error("Wrong mesh format, $MF, selected.")
    end
    tline = readline(fid)
    if is_not_string(tline) || !(tline == "\$EndMeshFormat")
        return nothing
    end
    return nothing
end

function read_node_line(tline) # A node line looks like "<node idx> <x> <y> <z> <ref>"
    all_vals = split(tline, ' ')
    x = parse(Float64, all_vals[2])
    y = parse(Float64, all_vals[3])
    return x, y
end

function read_num_objects(tline)
    return parse(Int64, tline)
end

function read_nodes!(nodes, fid)
    tline = readline(fid)
    is_not_string(tline) && return nothing
    no_nodes = read_num_objects(tline)
    for _ in 1:no_nodes
        tline = readline(fid)
        x, y = read_node_line(tline)
        push!(nodes, (x, y))
    end
    tline = readline(fid)
    if is_not_string(tline) || tline ≠ "\$EndNodes"
        return nothing
    end
    return nothing
end

function read_element_line(tline) # An element line (for a triangle, basic modification for a point, just see the code) looks like <elm idx> <elm type> <num tags> <physical> <node1> <node2> <node3>
    all_vals = split(tline, ' ')
    elm_type = parse(Int64, all_vals[2])
    physical_id = parse(Int64, all_vals[4])
    if elm_type == 1
        u = parse(Int64, all_vals[5])
        v = parse(Int64, all_vals[6])
        w = parse(Int64, all_vals[7])
    elseif elm_type == 2
        u = parse(Int64, all_vals[6])
        v = parse(Int64, all_vals[7])
        w = parse(Int64, all_vals[8])
    else
        throw("Invalid element type, $elm_type.")
    end
    return elm_type, physical_id, u, v, w
end

function read_elements!(elements, boundary_nodes, fid, identifier_dict)
    tline = readline(fid)
    is_not_string(tline) && return nothing
    no_objects = read_num_objects(tline)
    no_elements = 0
    for _ in 1:no_objects
        tline = readline(fid)
        elm_type, physical_id, u, v, w = read_element_line(tline)
        if elm_type == 2 # 2 means triangle
            no_elements += 1
            push!(elements, (u, v, w))
        end
        if physical_id ∈ Base.keys(identifier_dict)
            m, n = identifier_dict[physical_id]
            bn = boundary_nodes[m][n]
            if !isempty(bn) # In this case, "v" has already been added, so just add w
                append!(bn, w)
            else
                append!(bn, v, w)
            end
        end
    end
    tline = readline(fid)
    if is_not_string(tline) || tline ≠ "\$EndElements"
        return nothing
    end
    return nothing
end
