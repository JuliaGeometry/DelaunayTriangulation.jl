using DelaunayTriangulation
using Test
using SafeTestsets

include("helper_functions.jl")

@testset verbose = true "DelaunayTriangulation.jl" begin
    @testset verbose = true "Triangulation" begin
        @info "Testing Rectangular Triangulation"
        @safetestset "Rectangular Triangulation" begin
            include("triangulation/rectangle.jl")
        end
        @info "Testing Bowyer Watson"
        @safetestset "Bowyer Watson" begin
            include("triangulation/bowyer_watson.jl")
        end
        @info "Testing Triangulate"
        @safetestset "Triangulate" begin
            include("triangulation/triangulate.jl")
        end
        @info "Testing Convex Triangulation"
        @safetestset "Convex Triangulation" begin
            include("triangulation/convex_triangulation.jl")
        end
        @info "Testing Constrained Triangulation"
        @safetestset "Constrained Triangulation" begin
            include("triangulation/constrained.jl")
        end
        #  @info "Testing Weighted Triangulation"
        # Needs to be setup properly before we do more testing of it. Please see the 
        # comments at the end of the weighted.jl file (the one included in the comment below).
        #  @safetestset "Weighted Triangulation" begin
        #      include("triangulation/weighted.jl")
        #  end
        @info "Testing check_args" 
        @safetestset "check_args" begin
            include("triangulation/check_args.jl")
        end
    end

    @testset verbose = true "Interfaces" begin
        @info "Testing Triangles"
        @safetestset "Triangles" begin
            include("interfaces/triangles.jl")
        end
        @info "Testing Edges"
        @safetestset "Edges" begin
            include("interfaces/edges.jl")
        end
        @info "Testing Points"
        @safetestset "Points" begin
            include("interfaces/points.jl")
        end
        @info "Testing Boundary Nodes"
        @safetestset "Boundary Nodes" begin
            include("interfaces/boundary_nodes.jl")
        end
    end

    @testset verbose = true "Data Structures" begin
        @info "Testing Adjacent"
        @safetestset "Adjacent" begin
            include("data_structures/adjacent.jl")
        end
        @info "Testing Adjacent2Vertex"
        @safetestset "Adjacent2Vertex" begin
            include("data_structures/adjacent2vertex.jl")
        end
        @info "Testing Graph"
        @safetestset "Graph" begin
            include("data_structures/graph.jl")
        end
        @info "Testing ConvexHull"
        @safetestset "ConvexHull" begin
            include("data_structures/convex_hull.jl")
        end
        @info "Testing Triangulation"
        @safetestset "Triangulation" begin
            include("data_structures/triangulation.jl")
        end
        @info "Testing Representative Points"
        @safetestset "Representative Points" begin
            include("data_structures/representative.jl")
        end
        @info "Testing Statistics"
        @safetestset "Statistics" begin
            include("data_structures/statistics.jl")
        end
        @info "Testing ShuffledPolygonLinkedList"
        @safetestset "ShuffledPolygonLinkedList" begin
            include("data_structures/shuffled_polygon_linked_list.jl")
        end
        @info "Testing TriangulationCache"
        @safetestset "TriangulationCache" begin
            include("data_structures/triangulation_cache.jl")
        end
        @info "Testing MaxPriorityQueue"
        @safetestset "MaxPriorityQueue" begin
            include("data_structures/max_priority_queue.jl")
        end
        @info "Testing Queue"
        @safetestset "Queue" begin
            include("data_structures/queue.jl")
        end
        @info "Testing Curves"
        @safetestset "Curves" begin
            include("data_structures/curves.jl")
        end
        @info "Testing RTree"
        @safetestset "RTree" begin
            include("data_structures/rtree.jl")
        end
        @info "Testing BST"
        @safetestset "BST" begin
            include("data_structures/bst.jl")
        end
        @info "Testing PolygonHierarchy"
        @safetestset "PolygonHierarchy" begin
            include("data_structures/polygon_hierarchy.jl")
        end
    end

    @testset verbose = true "Predicates" begin
        @info "Testing Certificate"
        @safetestset "Certificate" begin
            include("predicates/certificate.jl")
        end
        @info "Testing Boundaries and Ghosts"
        @safetestset "Boundaries and Ghosts" begin
            include("predicates/boundaries_and_ghosts.jl")
        end
        @info "Testing General"
        @safetestset "General" begin
            include("predicates/general.jl")
        end
        @info "Testing Index and Ghost Handling"
        @safetestset "Index and Ghost Handling" begin
            include("predicates/index_and_ghost_handling.jl")
        end
    end

    @testset verbose = true "Utilities" begin
        @info "Testing Utilities"
        @safetestset "Utilities" begin
            include("utils.jl")
        end
        @info "Testing Geometric Utilities"
        @safetestset "Geometric Utilities" begin
            include("geo_utils.jl")
        end
        @info "Testing Triangulation Validation"
        @safetestset "Triangulation Validation" begin
            include("helper_function_tests.jl")
        end
    end

    @testset verbose = true "Point Location" begin
        @info "Testing Brute Force"
        @safetestset "Brute Force" begin
            include("point_location/brute_force.jl")
        end
        @info "Testing Initial Point Selection"
        @safetestset "Initial Point Selection" begin
            include("point_location/select_initial_point.jl")
        end
        @info "Testing Initial Triangle Selection"
        @safetestset "Initial Triangle Selection" begin
            include("point_location/select_initial_triangle_interior_node.jl")
        end
        @info "Testing Interior Edge Intersections"
        @safetestset "Interior Edge Intersections" begin
            include("point_location/interior_edge_intersections.jl")
        end
        @info "Testing Ghost Search"
        @safetestset "Ghost Search" begin
            include("point_location/ghost_search.jl")
        end
        @info "Testing Jump and March"
        @safetestset "Jump and March" begin
            include("point_location/jump_and_march.jl")
        end
        @info "Testing Polygon Location"
        @safetestset "Polygon Location" begin
            include("point_location/find_polygon.jl")
        end
    end

    @testset verbose = true "Operations" begin
        @info "Testing add_triangle!"
        @safetestset "add_triangle!" begin
            include("operations/add_triangle.jl")
        end
        @info "Testing delete_triangle!"
        @safetestset "delete_triangle!" begin
            include("operations/delete_triangle.jl")
        end
        @info "Testing add_ghost_triangles!"
        @safetestset "add_ghost_triangles!" begin
            include("operations/add_ghost_triangles.jl")
        end
        @info "Testing delete_ghost_triangles!"
        @safetestset "delete_ghost_triangles!" begin
            include("operations/delete_ghost_triangles.jl")
        end
        @info "Testing add_point!"
        @safetestset "add_point!" begin
            include("operations/add_point.jl")
        end
        @info "Testing flip_edge!"
        @safetestset "flip_edge!" begin
            include("operations/flip_edge.jl")
        end
        @info "Testing split_triangle!"
        @safetestset "split_triangle!" begin
            include("operations/split_triangle.jl")
        end
        @info "Testing split_edge!"
        @safetestset "split_edge!" begin
            include("operations/split_edge.jl")
        end
        @info "Testing legalise_edge!"
        @safetestset "legalise_edge!" begin
            include("operations/legalise_edge.jl")
        end
        @info "Testing delete_point!"
        @safetestset "delete_point!" begin
            include("operations/delete_point.jl")
        end
        @info "Testing (un)lock_convex_hull!.jl"
        @safetestset "(un)lock_convex_hull!.jl" begin
            include("operations/lock_convex_hull.jl")
        end
        @info "Testing delete_holes!"
        @safetestset "delete_holes!" begin
            include("operations/delete_holes.jl")
        end
    end

    @testset verbose = true "Constrained Triangulation" begin
        @info "Testing Segment Location"
        @safetestset "Segment Location" begin
            include("constrained_triangulation/segment_location.jl")
        end
        @info "Testing Segment Insertion"
        @safetestset "Segment Insertion" begin
            include("constrained_triangulation/segment_insertion.jl")
        end
    end

    @testset verbose = true "Refinement" begin
        @info "Testing Refinement"
        @safetestset "Refinement" begin 
            include("refinement/refine.jl")
        end
        @info "Testing Curve-Bounded Refinement"
        @safetestset "Curve-Bounded Refinement" begin
            include("refinement/curve_bounded.jl")
        end
    end

    @info "Testing Voronoi"
    @safetestset "Voronoi" begin
        include("voronoi/voronoi.jl")
    end
end