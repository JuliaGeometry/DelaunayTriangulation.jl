using DelaunayTriangulation
using Documenter
using Literate
using Test

DocMeta.setdocmeta!(DelaunayTriangulation, :DocTestSetup, :(using DelaunayTriangulation, Test);
    recursive=true)

const IS_GITHUB_ACTIONS = get(ENV, "GITHUB_ACTIONS", "false") == "true"
const IS_CI = get(ENV, "CI", "false") == "true"
function safe_include(filename)
    mod = @eval module $(gensym()) end
    return Base.include(mod, filename)
end

# When running docs locally, the EditURL is incorrect. For example, we might get 
#   ```@meta
#   EditURL = "<unknown>/docs/src/tutorials/constrained.jl"
#   ```
# We need to replace this EditURL if we are running the docs locally.
function update_edit_url(content)
    content = replace(content, "<unknown>" => "https://github.com/DanielVandH/DelaunayTriangulation.jl/tree/new-docs")
    return content
end

for folder in ("tutorials", "applications")
    dir = joinpath(@__DIR__, "src", folder)
    outputdir = joinpath(dir, "generated")
    files = readdir(dir)
    filter!(file -> endswith(file, ".jl"), files)
    for file in files
        file_path = joinpath(dir, file)
        @testset "$(file)" begin
            safe_include(file_path)
        end
        Literate.markdown(
            file_path,
            outputdir;
            documenter=true,
            postprocess=update_edit_url,
            credit=true
        )
    end
end

function relink_literate_files(name_path)
    section_name, path = name_path
    folder, filename = splitpath(path)
    if filename ∉ ("gmsh.md", "overview.md")
        return section_name => joinpath(folder, "generated", filename)
    end
    return name_path
end

const _PAGES = [
    "Introduction" => "index.md",
    "Tutorials" => relink_literate_files.([
        "Section Overview" => "tutorials/overview.md", # Introduction and installation
        "Unconstrained Triangulations" => "tutorials/unconstrained.md",
        "Constrained Triangulations" => "tutorials/constrained.md",
        "Dynamic Triangulations" => "tutorials/dynamic.md",
        "Mesh Refinement" => "tutorials/refinement.md",
        "Triangulating Rectangular Regions" => "tutorials/lattice.md",
        "Gmsh Integration" => "tutorials/gmsh.md",
        "Triangulating Convex Polygons" => "tutorials/convex.md",
        "Voronoi Tessellations" => "tutorials/voronoi.md",
        "Clipped Voronoi Tessellations" => "tutorials/clipped.md",
        "Centroidal Voronoi Tessellations" => "tutorials/centroidal.md",
        "Point Location" => "tutorials/point_location.md",
        "Nearest Neighbour Queries" => "tutorials/nearest.md",
        "Computing Convex Hulls" => "tutorials/convex_hull.md",
        "Computing the Pole of Inaccessibility" => "tutorials/pole_of_inaccessibility.md"
    ]),
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
    "Example Applications" => relink_literate_files.([
        "Section Overview" => "applications/overview.md",
        "Path-finding with Constrained Delaunay Triangulations" => "applications/pathfinding.md",
        "Interpolation" => "applications/interpolation.md", # also naturalneighbours.jl 
        "Cellular Biology" => "applications/cell_simulations.md",
        "PDE Discretisation" => "applications/pde_discretisation.md", # also finitevolumemethod.jl 
        "Image Compression" => "applications/image_compression.md", # see also https://d1wqtxts1xzle7.cloudfront.net/31255259/dfg99sirv-libre.pdf?1392197228=&response-content-disposition=inline%3B+filename%3DCentroidal_Voronoi_Tessellations_Applica.pdf&Expires=1690783507&Signature=d8s8javyhR743LoatXwziK84hklGFr77DE4Ns4DYcfm0ar19ZWZYlqRdZrUxzocNYZOa4oT4mrhh8WZ571BCa6-WDWQM4pNG0Zk0A9oZl4vuAzXBbKHLMt2cTXVms25Y7-bVBPYyQ8-YFNdTGg~5YibXW2kOxeoWcZo1JaBWYrOFezeg7DqZIY9smT0HtecVTHW1PjLUoJsnXbnTOF3My9NqXfY2ByXFWHcGb6U-KWvGntcHgnE8sxBdhAj9xPgehlbkygfIPY8mAmCbh7DIxcZ8HWKYaJfVvqTJOemFVx39dwi~Cwf-59eGBvFpnB2jUOVDsegPR40gz~Rqt3HCnA__&Key-Pair-Id=APKAJLOHF5GGSLRBV4ZA
        "Root Finding" => "applications/root_finding.md", # see also https://github.com/PioKow/GRPF and RootsAndPoles.jl
        "Counting Function Calls" => "applications/counting.md",
    ])
]

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
    pages=_PAGES
)

deploydocs(;
    repo="github.com/DanielVandH/DelaunayTriangulation.jl",
    devbranch="main")
