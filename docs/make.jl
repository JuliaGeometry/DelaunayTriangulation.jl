using DelaunayTriangulation
using Documenter

DocMeta.setdocmeta!(DelaunayTriangulation, :DocTestSetup, :(using DelaunayTriangulation);
    recursive=true)

makedocs(;
    modules=[DelaunayTriangulation],
    authors="Daniel VandenHeuvel <danj.vandenheuvel@gmail.com>",
    repo="https://github.com/DanielVandH/DelaunayTriangulation.jl/blob/{commit}{path}#{line}",
    sitename="DelaunayTriangulation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://DanielVandH.github.io/DelaunayTriangulation.jl",
        edit_link="main",
        assets=String[]),
    pages=["Home" => "index.md",
        "Triangulations" => [
            "Unconstrained Triangulations" => "triangulations/unconstrained.md",
            "Lattice" => "triangulations/lattice.md",
            "Gmsh" => "triangulations/gmsh.md",
            "Plotting" => "triangulations/plotting.md"
        ],
        "Triangulation Algorithms" => [
            "Bowyer-Watson Algorithm" => "tri_algs/bowyer.md"
        ],
        "Ghost Triangles and Boundary Handling" => "boundary_handling.md",
        "Other Features" => [
            "Point Location" => "other_features/point_location.md",
            "Pole of Inaccessibility and Polygons" => "other_features/pole_of_inaccessibility.md",
            "Convex Hull" => "other_features/convex_hull.md"
        ],
        "Data Structures" => [
            "Adjacent" => "data_structures/adjacent.md",
            "Adjacent2Vertex" => "data_structures/adjacent2vertex.md",
            "Graph" => "data_structures/graph.md",
            "Representative Coordinates" => "data_structures/representative.md",
            "Convex Hull" => "data_structures/convex_hull.md",
            "Triangulation" => "data_structures/triangulation.md"
        ],
        "Primitive Interfaces" => [
            "General and Defaults" => "interface/interface.md",
            "Triangles" => "interface/triangles.md",
            "Edges" => "interface/edges.md",
            "Points" => "interface/points.md",
            "Boundary Nodes" => "interface/boundary_nodes.md",
            "Example" => "interface/example.md"],
        "Predicates" => "predicates.md",
        "Operations" => "operations.md",
        "Other Utilities" => "utils.md"
    ])

deploydocs(;
    repo="github.com/DanielVandH/DelaunayTriangulation.jl",
    devbranch="main")
