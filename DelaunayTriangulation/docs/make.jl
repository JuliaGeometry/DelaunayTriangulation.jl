using DelaunayTriangulation
using Documenter
using Literate
using Test
using Dates

DocMeta.setdocmeta!(
    DelaunayTriangulation, :DocTestSetup, :(using DelaunayTriangulation, Test);
    recursive = true,
)

const IS_LIVESERVER = get(ENV, "LIVESERVER_ACTIVE", "false") == "true"
if IS_LIVESERVER
    using Revise
    Revise.revise()
end
const IS_GITHUB_ACTIONS = get(ENV, "GITHUB_ACTIONS", "false") == "true"
const IS_CI = get(ENV, "CI", "false") == "true"
function safe_include(filename)
    mod = @eval module $(gensym()) end
    return Base.include(mod, filename)
end
const session_tmp = mktempdir()

# When running docs locally, the EditURL is incorrect. For example, we might get 
#   ```@meta
#   EditURL = "<unknown>/docs/src/literate_tutorials/constrained.jl"
#   ```
# We need to replace this EditURL if we are running the docs locally. The last case is more complicated because, 
# after changing to use temporary directories, it can now look like...
#   ```@meta
#   EditURL = "../../../../../../../AppData/Local/Temp/jl_8nsMGu/cs1_just_the_code.jl"
#   ```
function update_edit_url(content, file, folder)
    content = replace(content, "<unknown>" => "https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main")
    content = replace(content, "temp/" => "") # as of Literate 2.14.1
    content = replace(content, r"EditURL\s*=\s*\"[^\"]*\"" => "EditURL = \"https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_$(folder)/$file\"")
    return content
end

# We can add the code to the end of each file in its uncommented form programatically.
function add_just_the_code_section(dir, file)
    file_name, file_ext = splitext(file)
    file_path = joinpath(dir, file)
    new_file_path = joinpath(session_tmp, file_name * "_just_the_code" * file_ext)
    cp(file_path, new_file_path, force = true)
    folder = splitpath(dir)[end] # literate_tutorials or literate_applications
    open(new_file_path, "a") do io
        write(io, "\n")
        write(io, "# ## Just the code\n")
        write(io, "# An uncommented version of this example is given below.\n")
        write(io, "# You can view the source code for this file [here](<unknown>/docs/src/$folder/@__NAME__.jl).\n")
        write(io, "\n")
        write(io, "# ```julia\n")
        write(io, "# @__CODE__\n")
        write(io, "# ```\n")
    end
    return new_file_path
end

# Now process all the literate files
ct() = Dates.format(now(), "HH:MM:SS")
for folder in ("tutorials", "applications")
    dir = joinpath(@__DIR__, "src", "literate_" * folder)
    outputdir = joinpath(@__DIR__, "src", folder)
    !isdir(outputdir) && mkpath(outputdir)
    files = readdir(dir)
    filter!(file -> endswith(file, ".jl") && !occursin("just_the_code", file), files)
    for file in files
        # See also https://github.com/Ferrite-FEM/Ferrite.jl/blob/d474caf357c696cdb80d7c5e1edcbc7b4c91af6b/docs/generate.jl for some of this
        file_path = joinpath(dir, file)
        #if !IS_LIVESERVER
        #    @info "[$(ct())] Processing $file: Testing"
        #    @testset "$(file)" begin
        #        safe_include(file_path)
        #    end
        #end
        new_file_path = add_just_the_code_section(dir, file)
        script = Literate.script(file_path, session_tmp, name = splitext(file)[1] * "_just_the_code_cleaned")
        code = strip(read(script, String))
        @info "[$(ct())] Processing $file: Converting markdown script"
        line_ending_symbol = occursin(code, "\r\n") ? "\r\n" : "\n"
        code_clean = join(filter(x -> !endswith(x, "#hide"), split(code, r"\n|\r\n")), line_ending_symbol)
        code_clean = replace(code_clean, r"^# This file was generated .*$"m => "")
        code_clean = strip(code_clean)
        post_strip = content -> replace(content, "@__CODE__" => code_clean)
        editurl_update = content -> update_edit_url(content, file, folder)
        Literate.markdown(
            new_file_path,
            outputdir;
            documenter = true,
            postprocess = editurl_update ∘ post_strip,
            credit = true,
            name = splitext(file)[1],
            execute = false,
            flavor = Literate.DocumenterFlavor(),
        )
    end
end

# All the pages to be included
const _PAGES = [
    "Introduction" => "index.md",
    "Tutorials" => [
        "Overview" => "tutorials/overview.md",
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
        ],
        "Mesh Refinement" => "tutorials/refinement.md",
        "Triangulating Rectangular Regions" => "tutorials/lattice.md",
        "Triangulating Convex Polygons" => "tutorials/convex.md",
        "Triangulating Curve-Bounded Domains" => "tutorials/curve_bounded.md",
        "Weighted Triangulations" => "tutorials/weighted.md",
        "Voronoi Tessellations" => "tutorials/voronoi.md",
        "Clipped Voronoi Tessellations" => [
            "Clipping to the Convex Hull" => "tutorials/clipped.md",
            "Clipping to a Rectangle" => "tutorials/clipped_rectangle.md",
            "Clipping to a Generic Convex Polygon" => "tutorials/clipped_polygon.md",
        ],
        "Centroidal Voronoi Tessellations" => "tutorials/centroidal.md",
        "Power Diagrams" => "tutorials/power.md",
        "Point Location" => "tutorials/point_location.md",
        "Nearest Neighbour Queries" => "tutorials/nearest.md",
        "Convex Hulls" => "tutorials/convex_hull.md",
        "Pole of Inaccessibility" => "tutorials/pole_of_inaccessibility.md",
        "Point-in-Polygon Testing" => "tutorials/point_in_polygon.md",
        "Using Custom Structs for Primitives and Boundaries" => "tutorials/custom_primitive.md",
    ],
    "Manual" => [
        "Overview" => "manual/overview.md",
        "Representing Primitives" => "manual/primitives.md",
        "Representing Boundaries" => "manual/boundaries.md",
        "Ghost (Negative) Vertices" => "manual/ghost_triangles.md",
        "Defining Curve-Bounded Domains" => "manual/curve_bounded.md",
        "Geometrical Predicates" => "manual/predicates.md",
        "Triangulation Output" => "manual/triangulation_output.md",
        "Voronoi Tessellation Output" => "manual/voronoi_output.md",
        "Predicate Kernels" => "manual/predicate_kernels.md",
    ],
    "API Reference" => [
        "Section Overview" => "api/overview.md",
        "Data Structures" => "api/data_structures.md",
        "Triangulations" => "api/triangulation.md",
        "Triangulation Operations" => "api/operations.md",
        "Voronoi Tessellations" => "api/voronoi.md",
        "Convex Hull" => "api/convex_hull.md",
        "Curves" => "api/curves.md",
        "Iterators" => "api/iterators.md",
        "Point Location" => "api/point_location.md",
        "Predicates" => "api/predicates.md",
        "Triangulation Statistics" => "api/statistics.md",
        "Primitive Interfaces" => "api/primitives.md",
        "Other" => "api/other.md",
    ],
    "Extended Reference" => [
        "Overview" => "extended/overview.md",
        "Data Structures" => "extended/data_structures.md",
        "Algorithm Internals" => "extended/algorithms.md",
        "Utility Functions" => "extended/utils.md",
    ],
    "Mathematical Details" => [
        "Overview" => "math/overview.md",
        "Delaunay Triangulations" => "math/delaunay.md",
        "Constrained Delaunay Triangulations" => "math/constrained.md",
        "Triangulating Convex Polygons" => "math/convex.md",
        "Mesh Refinement" => "math/refinement.md",
        "Curves" => "math/curves.md",
        "Triangulating Curve-Bounded Domains" => "math/curve_bounded.md",
        "Weighted Delaunay Triangulations" => "math/weighted.md",
        "Voronoi Tessellations" => "math/voronoi.md",
        "Clipped Voronoi Tessellations" => "math/clipped.md",
        "Centroidal Voronoi Tessellations" => "math/centroidal.md",
        "Power Diagrams" => "math/power.md",
    ],
    "Example Applications" => [
        "Overview" => "applications/overview.md",
        "Interpolation" => "applications/interpolation.md", # also naturalneighbours.jl 
        "Cellular Biology" => "applications/cell_simulations.md",
        "Solving PDEs" => "applications/pde_discretisation.md", # also finitevolumemethod.jl 
    ],
    "Terminology" => "terminology.md",
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
    modules = [DelaunayTriangulation],
    authors = "Daniel VandenHeuvel <danj.vandenheuvel@gmail.com>",
    sitename = "DelaunayTriangulation.jl",
    format = Documenter.HTML(;
        prettyurls = IS_CI,
        canonical = "https://JuliaGeometry.github.io/DelaunayTriangulation.jl",
        edit_link = "main",
        size_threshold = 8000 * 2^10,
        size_threshold_warn = 1000 * 2^10,
        size_threshold_ignore = ["api/triangulation.md", "extended/data_structures.md"],
        collapselevel = 1,
        assets = String[],
        mathengine = MathJax3(
            Dict(
                :loader => Dict("load" => ["[tex]/physics"]),
                :tex => Dict(
                    "inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
                    "tags" => "ams",
                    "packages" => ["base", "ams", "autoload", "physics"],
                ),
            ),
        ),
    ),
    linkcheck = false,
    warnonly = true,
    draft = IS_LIVESERVER,
    pages = _PAGES,
    pagesonly = true,
)

deploydocs(;
    repo = "github.com/JuliaGeometry/DelaunayTriangulation.jl",
    devbranch = "main",
    push_preview = true,
)
