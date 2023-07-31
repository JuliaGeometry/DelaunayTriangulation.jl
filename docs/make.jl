using DelaunayTriangulation
using Documenter

DocMeta.setdocmeta!(DelaunayTriangulation, :DocTestSetup, :(using DelaunayTriangulation, Test);
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
    pages=[
        "Introduction" => [
            "Introduction" => "introduction.md",
            "Installation" => "installation.md",
            "Referencing" => "referencing.md",
            "Similar Packages" => "similar.md",
        ],
        "Tutorials" => [
            "Installation and Overview" => "tutorials/installation.md", # Introduction and installation
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
        "Extended Manual" => [
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
        ]
    ])

#=
pages=["Home" => "index.md",
    "Triangulations" => [
        "Unconstrained Triangulations" => "triangulations/unconstrained.md",
        "Constrained Triangulations" => "triangulations/constrained.md",
        "Mesh Refinement" => "triangulations/refinement.md",
        "Lattice" => "triangulations/lattice.md",
        "Gmsh" => "triangulations/gmsh.md",
        "Plotting" => "triangulations/plotting.md",
        "Convex Polygons" => "triangulations/convex.md"
    ],
    "Voronoi Tessellations" => [
        "Voronoi Tessellations" => "tessellations/voronoi.md",
        "Clipped Voronoi Tessellations" => "tessellations/clipped.md",
        "Centroidal Voronoi Tessellation" => "tessellations/lloyd.md",
        "Plotting" => "tessellations/plotting.md"
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
=#

deploydocs(;
    repo="github.com/DanielVandH/DelaunayTriangulation.jl",
    devbranch="main")
