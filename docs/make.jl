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
            "Constrained Triangulations" => "triangulations/constrained.md",
            "Mesh Refinement" => "triangulations/refinement.md",
            "Lattice" => "triangulations/lattice.md",
            "Convex Polygons" => "triangulations/convex.md"
        ],
        "Voronoi Tessellations" => [
            "Voronoi Tessellations" => "tessellations/voronoi.md",
            "Clipped Voronoi Tessellations" => "tessellations/clipped.md",
            "Centroidal Voronoi Tessellation" => "tessellations/lloyd.md",
        ],
        "Boundary Handling" => "boundary_handling.md",
        "Data Structures" => [
            "Adjacent" => "data_structures/adjacent.md",
            "Adjacent2Vertex" => "data_structures/adjacent2vertex.md",
            "Graph" => "data_structures/graph.md",
            "Convex Hull" => "data_structures/convex_hull.md",
            "Triangulation" => "data_structures/triangulation.md",
            "Statistics" => "data_structures/statistics.md",
            "Voronoi Tessellation" => "data_structures/voronoi.md"
        ],
        "Operations" => "operations.md",
        "Other Features" => [
            "Point Location" => "other_features/point_location.md",
            "Pole of Inaccessibility and Polygons" => "other_features/pole_of_inaccessibility.md",
            "Convex Hull" => "other_features/convex_hull.md"
        ],
        "Primitive Interfaces" => [
            "General and Defaults" => "interface/interface.md",
            "Triangles" => "interface/triangles.md",
            "Edges" => "interface/edges.md",
            "Points" => "interface/points.md",
            "Boundary Nodes" => "interface/boundary_nodes.md",
            "Example" => "interface/example.md",
            "Application: Counting Function Calls" => "interface/counting.md"],
        "Predicates" => "predicates.md",
        "Other Utilities" => "utils.md",
        "Triangulation Algorithms" => [
            "Bowyer-Watson Algorithm" => "tri_algs/bowyer.md",
            "Chew's Algorithm for Triangulating Convex Polygons" => "tri_algs/convex.md",
            "Constrained Triangulations" => "tri_algs/constrained.md",
        ]
    ])

deploydocs(;
    repo="github.com/DanielVandH/DelaunayTriangulation.jl",
    devbranch="main",
    push_preview=true)
