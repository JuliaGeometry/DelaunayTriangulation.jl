using DelaunayTriangulation
using Documenter
using Literate
using Test

DocMeta.setdocmeta!(DelaunayTriangulation, :DocTestSetup, :(using DelaunayTriangulation, Test);
    recursive=true)

const IS_LIVESERVER = false
const CLEANUP_FIGURES = true
const IS_GITHUB_ACTIONS = get(ENV, "GITHUB_ACTIONS", "false") == "true"
const IS_CI = get(ENV, "CI", "false") == "true"
function safe_include(filename)
    mod = @eval module $(gensym()) end
    return Base.include(mod, filename)
end

function clean_dir()
    for folder in ("tutorials", "applications")
        dir = joinpath(@__DIR__, "src", folder)
        files = readdir(dir)
        filter!(file -> endswith(file, ".md"), files)
        filter!(file -> file ∉ ("overview.md", "gmsh.md"), files)
        for file in files
            file_path = joinpath(dir, file)
            rm(file_path)
        end
        if !IS_LIVESERVER # draft documentation doesn't generate these files 
            generated_script_path = joinpath(dir, "generated")
            temp_script_path = joinpath(dir, "temp")
            for file in readdir(generated_script_path)
                rm(joinpath(generated_script_path, file))
            end
            for file in readdir(temp_script_path)
                rm(joinpath(temp_script_path, file))
            end
        end
        if CLEANUP_FIGURES
            fig_path = joinpath(@__DIR__, "src", folder, "figures")
            for file in readdir(fig_path)
                rm(joinpath(fig_path, file))
            end
        end
    end
end
clean_dir()

# When running docs locally, the EditURL is incorrect. For example, we might get 
#   ```@meta
#   EditURL = "<unknown>/docs/src/tutorials/constrained.jl"
#   ```
# We need to replace this EditURL if we are running the docs locally.
function update_edit_url(content)
    content = replace(content, "<unknown>" => "https://github.com/DanielVandH/DelaunayTriangulation.jl/tree/new-docs")
    content = replace(content, "temp/" => "") # as of Literate 2.14.1
    return content
end

# We can add the code to the end of each file in its uncommented form programatically 
function add_just_the_code_section(dir, file)
    file_path = joinpath(dir, file)
    new_file_path = joinpath(dir, "temp", file)
    cp(file_path, new_file_path)
    folder = splitpath(dir)[end] # tutorials or applications
    open(new_file_path, "a") do io
        write(io, "\n")
        write(io, "# ## Just the code\n")
        write(io, "# An uncommented version of this tutorial is given below.\n")
        write(io, "# You can view the source code for this [here](<unknown>/docs/src/$folder/@__NAME__.jl).\n")
        write(io, "\n")
        write(io, "# ```julia\n")
        write(io, "# @__CODE__\n")
        write(io, "# ```\n")
    end
    return new_file_path
end

# Now process all the literate files
for folder in ("tutorials", "applications")
    dir = joinpath(@__DIR__, "src", folder)
    outputdir = dir
    files = readdir(dir)
    filter!(file -> endswith(file, ".jl"), files)
    for file in files
        # See also https://github.com/Ferrite-FEM/Ferrite.jl/blob/d474caf357c696cdb80d7c5e1edcbc7b4c91af6b/docs/generate.jl for some of this
        file_path = joinpath(dir, file)
        file_path = add_just_the_code_section(dir, file)
        if !IS_LIVESERVER
            @testset "$(file)" begin
                safe_include(file_path)
            end
            script = Literate.script(file_path, joinpath(outputdir, "generated"))
            code = strip(read(script, String))
        else
            code = "Script outputs are not produced when producing draft documentation."
        end
        line_ending_symbol = occursin(code, "\r\n") ? "\r\n" : "\n"
        code_clean = join(filter(x -> !endswith(x, "#hide"), split(code, r"\n|\r\n")), line_ending_symbol)
        code_clean = replace(code_clean, r"^# This file was generated .*$"m => "")
        code_clean = strip(code_clean)
        post_strip(content) = replace(content, "@__CODE__" => code_clean)
        Literate.markdown(
            file_path,
            outputdir;
            documenter=true,
            postprocess=update_edit_url ∘ post_strip,
            credit=true
        )
    end
end

# All the pages to be included
const _PAGES = [
    "Introduction" => "index.md",
    "Tutorials" => [
        "Section Overview" => "tutorials/overview.md", # Introduction and installation
        "Unconstrained Triangulations" => "tutorials/unconstrained.md",
        "Constrained Triangulations" => [
            "Constrained Edges" => "tutorials/constrained_edges.md",
            "Outer Boundary" => "tutorials/constrained_outer_boundary.md",
            "Segmented Outer Boundary" => "tutorials/constrained_outer_boundary_segmented.md",
            "Domain with Interior Holes" => "tutorials/constrained_multiply_connected.md",
            "Domain with Interior Holes inside Interior Holes" => "tutorials/constrained_interior_within_interiors.md",
            "Disjoint Domains" => "tutorials/constrained_multipolygon.md",
        ],
        "Triangulation Operations" => [
            "Vertex Insertion and Deletion" => "tutorials/operations_vertex_insertion_deletion.md",
            "Edge Insertion" => "tutorials/operations_segment_insertion.md",
            "Adding or Clearing Ghost Triangles" => "tutorials/operations_ghost_triangles.md",
            "Edge Flipping" => "tutorials/operations_flip_edge.md",
            "Legalising an Edge" => "tutorials/operations_legalise_edge.md",
            "Edge Splitting" => "tutorials/operations_split_edge.md",
            "Triangle Splitting" => "tutorials/operations_split_triangle.md",
            "Locking and Unlocking the Convex Hull" => "tutorials/operations_convex_hull_locking.md",
            "Clearing Empty Features" => "tutorials/operations_feature_clearing.md",
        ],
        "Mesh Refinement" => "tutorials/refinement.md",
        "Triangulating Rectangular Regions" => "tutorials/lattice.md",
        "Gmsh Integration" => "tutorials/gmsh.md",
        "Triangulating Convex Polygons" => "tutorials/convex.md",
        "Voronoi Tessellations" => "tutorials/voronoi.md",
        "Clipped Voronoi Tessellations" => "tutorials/clipped.md",
        "Centroidal Voronoi Tessellations" => "tutorials/centroidal.md",
        "Point Location" => "tutorials/point_location.md",
        "Nearest Neighbour Queries" => "tutorials/nearest.md",
        "Convex Hulls" => "tutorials/convex_hull.md",
        "Pole of Inaccessibility" => "tutorials/pole_of_inaccessibility.md",
        "Using Custom Structs for Primitives and Boundaries" => "tutorials/custom_primitive.md"
    ],
    "Manual" => [
        "Section Overview" => "manual/overview.md",
        "Representing Primitives" => "manual/primitives.md",
        "Representing Boundaries" => "manual/boundaries.md",
        "Ghost Triangles" => "manual/ghost_triangles.md",
        "Geometrical Predicates" => "manual/predicates.md",
        "Triangulation Output" => "manual/triangulation_output.md",
        "Voronoi Tessellation Output" => "manual/voronoi_output.md",
    ],
    "API Reference" => [
        "Section Overview" => "api/overview.md",
        "List of Public Functions" => "api/public.md",
        "Primitive Interfaces" => "api/primitive_interfaces.md",
        "Triangulations" => "api/triangulations.md",
        "Triangulation Operations" => "api/operations.md",
        "Mesh Refinement" => "api/refinement.md",
        "Voronoi Tessellations" => "api/voronoi.md",
    ],
    "Extended Reference" => [
        "Section Overview" => "extended/overview.md",
        "All Data Structures" => "extended/data_structures.md",
        "Utility Functions" => "extended/utils.md",
    ],
    "Mathematical and Implementation Details" => [
        "Section Overview" => "math/overview.md",
        "Delaunay Triangulations" => "math/delaunay.md",
        "Constrained Delaunay Triangulations" => "math/constrained.md",
        "Triangulating Convex Polygons" => "math/convex.md",
        "Mesh Refinement" => "math/refinement.md",
        "Voronoi Tessellations" => "math/voronoi.md",
        "Clipped Voronoi Tessellations" => "math/clipped.md",
        "Centroidal Voronoi Tessellations" => "math/centroidal.md",
        "Point Location" => "math/point_location.md",
        "Nearest Neighbour Queries" => "math/nearest.md",
        "Convex Hulls" => "math/convex_hull.md",
        "Pole of Inaccessibility" => "math/pole_of_inaccessibility.md",
        "Triangulation Operations" => "math/operations.md",
    ],
    "Example Applications" => [
        "Section Overview" => "applications/overview.md",
        "Path-finding with Constrained Delaunay Triangulations" => "applications/pathfinding.md",
        "Interpolation" => "applications/interpolation.md", # also naturalneighbours.jl 
        "Cellular Biology" => "applications/cell_simulations.md",
        "PDE Discretisation" => "applications/pde_discretisation.md", # also finitevolumemethod.jl 
        "Image Compression" => "applications/image_compression.md", # see also https://d1wqtxts1xzle7.cloudfront.net/31255259/dfg99sirv-libre.pdf?1392197228=&response-content-disposition=inline%3B+filename%3DCentroidal_Voronoi_Tessellations_Applica.pdf&Expires=1690783507&Signature=d8s8javyhR743LoatXwziK84hklGFr77DE4Ns4DYcfm0ar19ZWZYlqRdZrUxzocNYZOa4oT4mrhh8WZ571BCa6-WDWQM4pNG0Zk0A9oZl4vuAzXBbKHLMt2cTXVms25Y7-bVBPYyQ8-YFNdTGg~5YibXW2kOxeoWcZo1JaBWYrOFezeg7DqZIY9smT0HtecVTHW1PjLUoJsnXbnTOF3My9NqXfY2ByXFWHcGb6U-KWvGntcHgnE8sxBdhAj9xPgehlbkygfIPY8mAmCbh7DIxcZ8HWKYaJfVvqTJOemFVx39dwi~Cwf-59eGBvFpnB2jUOVDsegPR40gz~Rqt3HCnA__&Key-Pair-Id=APKAJLOHF5GGSLRBV4ZA
        "Root Finding" => "applications/root_finding.md", # see also https://github.com/PioKow/GRPF and RootsAndPoles.jl
        "Counting Function Calls" => "applications/counting.md",
        "Edge Flipping Algorithm" => "applications/edge_flipping.md",
        "de Berg's Algorithm" => "applications/de_berg.md",
    ]
]

# Make sure we haven't forgotten any files
set = Set{String}()
for page in _PAGES
    if page[2] isa String
        push!(set, normpath(page[2]))
    else
        for _page in page[2]
            if _page[2] isa String
                push!(set, normpath(_page[2]))
            else
                for __page in _page[2]
                    push!(set, normpath(__page[2]))
                end
            end
        end
    end
end
missing_set = String[]
doc_dir = joinpath(@__DIR__, "src", "")
for (root, dir, files) in walkdir(doc_dir)
    for file in files
        filename = normpath(replace(joinpath(root, file), doc_dir => ""))
        if endswith(filename, ".md") && filename ∉ set
            push!(missing_set, filename)
        end
    end
end
!isempty(missing_set) && error("Missing files: $missing_set")

# Make and deploy
makedocs(;
    modules=[DelaunayTriangulation],
    authors="Daniel VandenHeuvel <danj.vandenheuvel@gmail.com>",
    repo="https://github.com/DanielVandH/DelaunayTriangulation.jl/blob/{commit}{path}#{line}",
    sitename="DelaunayTriangulation.jl",
    format=Documenter.HTML(;
        prettyurls=IS_CI,
        canonical="https://DanielVandH.github.io/DelaunayTriangulation.jl",
        edit_link="main",
        collapselevel=1,
        assets=String[]),
    linkcheck=false,
    strict=false,
    draft=IS_LIVESERVER,
    pages=_PAGES
)

deploydocs(;
    repo="github.com/DanielVandH/DelaunayTriangulation.jl",
    devbranch="main")

# Now that we are done, delete the literate files 
clean_dir()