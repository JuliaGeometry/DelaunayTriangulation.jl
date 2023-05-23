using Test
using LinearAlgebra
using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using DataStructures

include("../helper_functions.jl")

@testset "Working with identifiers" begin
    @testset "Computing segment identifiers" begin
        x = [[zeros(120), zeros(127)], [zeros(50)], [zeros(51), zeros(50)], [zeros(1000)]]
        ids = DT.get_segment_identifiers(x)
        @test ids == [[1, 2], [3], [4, 5], [6]]
    end

    @testset "Computing Dict of segment identifiers" begin
        dict = DT.get_segment_identifiers_dict([[1, 2], [3], [4, 5, 6], [7, 8]], 13:20)
        @test dict == Dict(13 => (1, 1),
            14 => (1, 2),
            15 => (2, 1),
            16 => (3, 1),
            17 => (3, 2),
            18 => (3, 3),
            19 => (4, 1),
            20 => (4, 2))
    end
end

@testset "Initialising the boundary node vectors" begin
    bn = DT.get_boundary_node_init([[1, 2], [3], [4, 5, 6]])
    @test bn == [[Int[], Int[]],
        [Int[]],
        [Int[], Int[], Int[]]]
end

@testset "Defining mesh settings" begin
    fname = tempname()
    open(fname, "w") do fout
        return DT.write_mesh_settings!(fout, 0.1, 6, 0)
    end
@test read(fname, String) == "r = 0.1;
Mesh.Algorithm = 6;
Mesh.Format = 1;
General.Verbosity = 0;
"
end

@testset "Writing points" begin
    @testset "Writing single points" begin
        fname = tempname()
        open(fname, "w") do fout
            DT.write_point!(fout, 1, 0.5, 3.2)
            DT.write_point!(fout, 2, -0.3, 0)
            return DT.write_point!(fout, 5, 2.3333292, 0.222)
        end
@test read(fname, String) == "Point(1) = {0.5, 3.2, 0, r};
Point(2) = {-0.3, 0, 0, r};
Point(5) = {2.3333292, 0.222, 0, r};
"
    end

    @testset "Writing a vector of points" begin
        x = [0.3, 0.22, -0.2, -0.5, 1.3]
        y = [0.0, 5.3, 2.6, 4.3, 0.222]
        fname = tempname()
        open(fname, "w") do fout
            return DT.write_points!(fout, x, y)
        end
@test read(fname, String) == "Point(1) = {0.3, 0.0, 0, r};
Point(2) = {0.22, 5.3, 0, r};
Point(3) = {-0.2, 2.6, 0, r};
Point(4) = {-0.5, 4.3, 0, r};
Point(5) = {1.3, 0.222, 0, r};
"
        xx = [2.3, 2.22, -0.22, -3.5, 0.3]
        yy = [5.0, 1.2, 7.3, 4.377, -0.3102]
        open(fname, "w") do fout
            DT.write_points!(fout, x, y)
            return DT.write_points!(fout, xx, yy, 71)
        end
@test read(fname, String) == "Point(1) = {0.3, 0.0, 0, r};
Point(2) = {0.22, 5.3, 0, r};
Point(3) = {-0.2, 2.6, 0, r};
Point(4) = {-0.5, 4.3, 0, r};
Point(5) = {1.3, 0.222, 0, r};
Point(71) = {2.3, 5.0, 0, r};
Point(72) = {2.22, 1.2, 0, r};
Point(73) = {-0.22, 7.3, 0, r};
Point(74) = {-3.5, 4.377, 0, r};
Point(75) = {0.3, -0.3102, 0, r};
"
    end

    @testset "Writing a vector of vector of points" begin
        x = [[0.3, -0.5, 1.0], [1.2, 5.3], [1.0], [2.0, 3.3, 5.3, 10.0]]
        y = [[4.2, 2.2, 0.5], [-3.3, 0.0], [2.3], [5.3, 13.3, 20.9, 2.9]]
        fname = tempname()
        open(fname, "w") do fout
            return DT.write_points!(fout, x, y)
        end
@test read(fname, String) == "Point(1) = {0.3, 4.2, 0, r};
Point(2) = {-0.5, 2.2, 0, r};
Point(3) = {1.0, 0.5, 0, r};
Point(4) = {1.2, -3.3, 0, r};
Point(5) = {5.3, 0.0, 0, r};
Point(6) = {1.0, 2.3, 0, r};
Point(7) = {2.0, 5.3, 0, r};
Point(8) = {3.3, 13.3, 0, r};
Point(9) = {5.3, 20.9, 0, r};
Point(10) = {10.0, 2.9, 0, r};
"
        x = [[2.3, 5.0], [1.0], [2.3], [2.5, 13.3, 20.9, 20.0, 21.1, -3.2], [2.0, 3.0, 4.0]]
        y = [[3.5, 5.1], [1.2], [12.0], [7.7, 0.1, -0.3, 2.2, 1.3, 2.5], [5.5, 2.3, -19.9]]
        fname = tempname()
        open(fname, "w") do fout
            return DT.write_points!(fout, x, y, 17)
        end
@test read(fname, String) == "Point(17) = {2.3, 3.5, 0, r};
Point(18) = {5.0, 5.1, 0, r};
Point(19) = {1.0, 1.2, 0, r};
Point(20) = {2.3, 12.0, 0, r};
Point(21) = {2.5, 7.7, 0, r};
Point(22) = {13.3, 0.1, 0, r};
Point(23) = {20.9, -0.3, 0, r};
Point(24) = {20.0, 2.2, 0, r};
Point(25) = {21.1, 1.3, 0, r};
Point(26) = {-3.2, 2.5, 0, r};
Point(27) = {2.0, 5.5, 0, r};
Point(28) = {3.0, 2.3, 0, r};
Point(29) = {4.0, -19.9, 0, r};
"
    end

    ## Writing a vector of a vector of a vector of points 
    @testset "Writing a vector of vector of vector of points" begin
        x1 = [[0.3, -0.5, 1.0], [1.2, 5.3], [1.0], [2.0, 3.3, 5.3, 10.0]]
        y1 = [[4.2, 2.2, 0.5], [-3.3, 0.0], [2.3], [5.3, 13.3, 20.9, 2.9]]
        x2 = [[2.3, 5.0], [1.0], [2.3], [2.5, 13.3, 20.9, 20.0, 21.1, -3.2], [2.0, 3.0, 4.0]]
        y2 = [[3.5, 5.1], [1.2], [12.0], [7.7, 0.1, -0.3, 2.2, 1.3, 2.5], [5.5, 2.3, -19.9]]
        x = [x1, x2]
        y = [y1, y2]
        fname = tempname()
        open(fname, "w") do fout
            return DT.write_points!(fout, x, y)
        end
@test read(fname, String) == "Point(1) = {0.3, 4.2, 0, r};
Point(2) = {-0.5, 2.2, 0, r};
Point(3) = {1.0, 0.5, 0, r};
Point(4) = {1.2, -3.3, 0, r};
Point(5) = {5.3, 0.0, 0, r};
Point(6) = {1.0, 2.3, 0, r};
Point(7) = {2.0, 5.3, 0, r};
Point(8) = {3.3, 13.3, 0, r};
Point(9) = {5.3, 20.9, 0, r};
Point(10) = {10.0, 2.9, 0, r};
Point(11) = {2.3, 3.5, 0, r};
Point(12) = {5.0, 5.1, 0, r};
Point(13) = {1.0, 1.2, 0, r};
Point(14) = {2.3, 12.0, 0, r};
Point(15) = {2.5, 7.7, 0, r};
Point(16) = {13.3, 0.1, 0, r};
Point(17) = {20.9, -0.3, 0, r};
Point(18) = {20.0, 2.2, 0, r};
Point(19) = {21.1, 1.3, 0, r};
Point(20) = {-3.2, 2.5, 0, r};
Point(21) = {2.0, 5.5, 0, r};
Point(22) = {3.0, 2.3, 0, r};
Point(23) = {4.0, -19.9, 0, r};
"
        fname = tempname()
        open(fname, "w") do fout
            return DT.write_points!(fout, x, y, 23)
        end
@test read(fname, String) == "Point(23) = {0.3, 4.2, 0, r};
Point(24) = {-0.5, 2.2, 0, r};
Point(25) = {1.0, 0.5, 0, r};
Point(26) = {1.2, -3.3, 0, r};
Point(27) = {5.3, 0.0, 0, r};
Point(28) = {1.0, 2.3, 0, r};
Point(29) = {2.0, 5.3, 0, r};
Point(30) = {3.3, 13.3, 0, r};
Point(31) = {5.3, 20.9, 0, r};
Point(32) = {10.0, 2.9, 0, r};
Point(33) = {2.3, 3.5, 0, r};
Point(34) = {5.0, 5.1, 0, r};
Point(35) = {1.0, 1.2, 0, r};
Point(36) = {2.3, 12.0, 0, r};
Point(37) = {2.5, 7.7, 0, r};
Point(38) = {13.3, 0.1, 0, r};
Point(39) = {20.9, -0.3, 0, r};
Point(40) = {20.0, 2.2, 0, r};
Point(41) = {21.1, 1.3, 0, r};
Point(42) = {-3.2, 2.5, 0, r};
Point(43) = {2.0, 5.5, 0, r};
Point(44) = {3.0, 2.3, 0, r};
Point(45) = {4.0, -19.9, 0, r};
"
    end
end

@testset "Writing lines" begin
    @testset "Writing a line" begin
        fname = tempname()
        open(fname, "w") do fout
            DT.write_line!(fout, 7, 2, 3)
            return DT.write_line!(fout, 1, 1, 2)
        end
@test read(fname, String) == "Line(7) = {2, 3};
Line(1) = {1, 2};
"
    end

    @testset "Writing a set of lines from a vector of numbers" begin
        x = rand(5)
        y = rand(5)
        fname = tempname()
        open(fname, "w") do fout
            return DT.write_lines!(fout, x, y)
        end
@test read(fname, String) == "Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
"
        fname = tempname()
        open(fname, "w") do fout
            return DT.write_lines!(fout, x, y, 2, 7)
        end
@test read(fname, String) == "Line(7) = {2, 3};
Line(8) = {3, 4};
Line(9) = {4, 5};
Line(10) = {5, 6};
"
    end

    @testset "Writing a set of lines from a vector of vector of numbers" begin
        @testset "Single segment" begin
            x = [rand(5)]
            y = [rand(5)]
            fname = tempname()
            open(fname, "w") do fout
                global segment_ranges
                return _, _, segment_ranges = DT.write_segments!(fout, x, y)
            end
@test read(fname, String) == "Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
"
            @test segment_ranges == [1:4]

            x = [rand(7)]
            y = [rand(7)]
            fname = tempname()
            open(fname, "w") do fout
                global segment_ranges
                return _, _, segment_ranges = DT.write_segments!(fout, x, y, 5, 10)
            end
@test read(fname, String) == "Line(10) = {5, 6};
Line(11) = {6, 7};
Line(12) = {7, 8};
Line(13) = {8, 9};
Line(14) = {9, 10};
Line(15) = {10, 5};
"
            @test segment_ranges == [10:15]
        end

        @testset "Two segments" begin
            x = [rand(5), rand(5)]
            y = [rand(5), rand(5)]
            fname = tempname()
            open(fname, "w") do fout
                global segment_ranges
                return _, _, segment_ranges = DT.write_segments!(fout, x, y)
            end
@test read(fname, String) == "Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 7};
Line(6) = {7, 8};
Line(7) = {8, 9};
Line(8) = {9, 1};
"
            @test segment_ranges == [1:4, 5:8]

            x = [rand(3), rand(5)]
            y = [rand(3), rand(5)]
            fname = tempname()
            open(fname, "w") do fout
                global segment_ranges
                return _, _, segment_ranges = DT.write_segments!(fout, x, y)
            end
@test read(fname, String) == "Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 5};
Line(4) = {5, 6};
Line(5) = {6, 7};
Line(6) = {7, 1};
"
            @test segment_ranges == [1:2, 3:6]

            x = [rand(5), rand(4)]
            y = [rand(5), rand(4)]
            fname = tempname()
            open(fname, "w") do fout
                global segment_ranges
                return _, _, segment_ranges = DT.write_segments!(fout, x, y, 10, 5)
            end
@test read(fname, String) == "Line(5) = {10, 11};
Line(6) = {11, 12};
Line(7) = {12, 13};
Line(8) = {13, 14};
Line(9) = {14, 16};
Line(10) = {16, 17};
Line(11) = {17, 10};
"
            @test segment_ranges == [5:8, 9:11]
        end

        @testset "Three segments" begin
            x = [rand(5), rand(5), rand(5)]
            y = [rand(5), rand(5), rand(5)]
            fname = tempname()
            open(fname, "w") do fout
                global segment_ranges
                return _, _, segment_ranges = DT.write_segments!(fout, x, y)
            end
@test read(fname, String) == "Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 7};
Line(6) = {7, 8};
Line(7) = {8, 9};
Line(8) = {9, 10};
Line(9) = {10, 12};
Line(10) = {12, 13};
Line(11) = {13, 14};
Line(12) = {14, 1};
"
            @test segment_ranges == [1:4, 5:8, 9:12]

            x = [rand(6), rand(4), rand(3)]
            y = [rand(6), rand(4), rand(3)]
            fname = tempname()
            open(fname, "w") do fout
                global segment_ranges
                return _, _, segment_ranges = DT.write_segments!(fout, x, y)
            end
@test read(fname, String) == "Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 8};
Line(7) = {8, 9};
Line(8) = {9, 10};
Line(9) = {10, 12};
Line(10) = {12, 1};
"
            @test segment_ranges == [1:5, 6:8, 9:10]
            fname = tempname()
            open(fname, "w") do fout
                global segment_ranges
                return _, _, segment_ranges = DT.write_segments!(fout, x, y, 2, 7)
            end
@test read(fname, String) == "Line(7) = {2, 3};
Line(8) = {3, 4};
Line(9) = {4, 5};
Line(10) = {5, 6};
Line(11) = {6, 7};
Line(12) = {7, 9};
Line(13) = {9, 10};
Line(14) = {10, 11};
Line(15) = {11, 13};
Line(16) = {13, 2};
"
            @test segment_ranges == [7:11, 12:14, 15:16]
        end
    end

    @testset "Writing a set of lines from vector of vector of vector of numbers" begin
        @testset "One curve, one segment" begin
            x = [[rand(5)]]
            y = [[rand(5)]]
            fname = tempname()
            open(fname, "w") do fout
                global curve_range, segment_range, nc
                return curve_range, segment_range, nc = DT.write_curves!(fout, x, y)
            end
@test read(fname, String) == "Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
"
            @test curve_range == [1:4]
            @test segment_range == [[1:4]]
            @test nc == 4
        end

        @testset "One curve, two segments" begin
            x = [[rand(5), rand(4)]]
            y = [[rand(5), rand(4)]]
            fname = tempname()
            open(fname, "w") do fout
                global curve_range, segment_range, nc
                return curve_range, segment_range, nc = DT.write_curves!(fout, x, y, 10, 5)
            end
@test read(fname, String) == "Line(5) = {10, 11};
Line(6) = {11, 12};
Line(7) = {12, 13};
Line(8) = {13, 14};
Line(9) = {14, 16};
Line(10) = {16, 17};
Line(11) = {17, 10};
"
            @test segment_range == [[5:8, 9:11]]
            @test curve_range == [5:11]
            @test nc == 11
        end

        @testset "Two curves, one segment each" begin
            x = [[rand(5)], [rand(4)]] # [[1, 2, 3, 4, 1], [6, 7, 8, 6]]
            y = [[rand(5)], [rand(4)]]
            fname = tempname()
            open(fname, "w") do fout
                global curve_range, segment_range, nc
                return curve_range, segment_range, nc = DT.write_curves!(fout, x, y)
            end
@test read(fname, String) == "Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {6, 7};
Line(6) = {7, 8};
Line(7) = {8, 6};
"
            @test curve_range == [1:4, 5:7]
            @test segment_range == [[1:4], [5:7]]
            @test nc == 7
        end

        @testset "Two curves, differing segments" begin
            x = [[rand(5)], [rand(4), rand(10)]] ## [[1, 2, 3, 4, 1], [6, 7, 8, 9, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]]
            y = [[rand(5)], [rand(4), rand(10)]]
            fname = tempname()
            open(fname, "w") do fout
                global curve_range, segment_range, nc
                return curve_range, segment_range, nc = DT.write_curves!(fout, x, y)
            end
@test read(fname, String) == "Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {6, 7};
Line(6) = {7, 8};
Line(7) = {8, 9};
Line(8) = {9, 11};
Line(9) = {11, 12};
Line(10) = {12, 13};
Line(11) = {13, 14};
Line(12) = {14, 15};
Line(13) = {15, 16};
Line(14) = {16, 17};
Line(15) = {17, 18};
Line(16) = {18, 6};
"
            @test curve_range == [1:4, 5:16]
            @test segment_range == [[1:4], [5:7, 8:16]]
            @test nc == 16

            x = [[rand(5)], [rand(4), rand(10), rand(6)]]
            y = [[rand(5)], [rand(4), rand(10), rand(6)]]
            fname = tempname()
            open(fname, "w") do fout
                global curve_range, segment_range, nc
                return curve_range, segment_range, nc = DT.write_curves!(fout, x, y)
            end
@test read(fname, String) == "Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {6, 7};
Line(6) = {7, 8};
Line(7) = {8, 9};
Line(8) = {9, 11};
Line(9) = {11, 12};
Line(10) = {12, 13};
Line(11) = {13, 14};
Line(12) = {14, 15};
Line(13) = {15, 16};
Line(14) = {16, 17};
Line(15) = {17, 18};
Line(16) = {18, 19};
Line(17) = {19, 21};
Line(18) = {21, 22};
Line(19) = {22, 23};
Line(20) = {23, 24};
Line(21) = {24, 6};
"
            @test curve_range == [1:4, 5:21]
            @test segment_range == [[1:4], [5:7, 8:16, 17:21]]
            @test nc == 21
        end

        @testset "Another two curve example" begin
            x = [[rand(5), rand(5)], [rand(5)]]
            y = [[rand(5), rand(5)], [rand(5)]]
            fname = tempname()
            open(fname, "w") do fout
                return DT.write_curves!(fout, x, y)
            end
@test read(fname, String) == "Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 7};
Line(6) = {7, 8};
Line(7) = {8, 9};
Line(8) = {9, 1};
Line(9) = {11, 12};
Line(10) = {12, 13};
Line(11) = {13, 14};
Line(12) = {14, 11};
"
        end

        @testset "Three curves" begin
            x1 = [rand(5)]
            x2 = [rand(3), rand(5)]
            x3 = [rand(4), rand(5), rand(10)]
            y1 = [rand(5)]
            y2 = [rand(3), rand(5)]
            y3 = [rand(4), rand(5), rand(10)]
            x = [x1, x2, x3]
            y = [y1, y2, y3]
            fname = tempname()
            open(fname, "w") do fout
                global curve_range, segment_range, nc
                return curve_range, segment_range, nc = DT.write_curves!(fout, x, y)
            end
@test read(fname, String) == "Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {6, 7};
Line(6) = {7, 8};
Line(7) = {8, 10};
Line(8) = {10, 11};
Line(9) = {11, 12};
Line(10) = {12, 6};
Line(11) = {14, 15};
Line(12) = {15, 16};
Line(13) = {16, 17};
Line(14) = {17, 19};
Line(15) = {19, 20};
Line(16) = {20, 21};
Line(17) = {21, 22};
Line(18) = {22, 24};
Line(19) = {24, 25};
Line(20) = {25, 26};
Line(21) = {26, 27};
Line(22) = {27, 28};
Line(23) = {28, 29};
Line(24) = {29, 30};
Line(25) = {30, 31};
Line(26) = {31, 14};
"
            @test curve_range == [1:4, 5:10, 11:26]
            @test segment_range == [[1:4], [5:6, 7:10], [11:13, 14:17, 18:26]]
            @test nc == 26
        end

        @testset "Many curves" begin
            x1 = [rand(5)]
            x2 = [rand(3), rand(5)]
            x3 = [rand(4), rand(5), rand(10)]
            x4 = [rand(8), rand(6)]
            x5 = [rand(6)]
            y1 = [rand(5)]
            y2 = [rand(3), rand(5)]
            y3 = [rand(4), rand(5), rand(10)]
            y4 = [rand(8), rand(6)]
            y5 = [rand(6)]
            x = [x1, x2, x3, x4, x5]
            y = [y1, y2, y3, y4, y5]
            fname = tempname()
            open(fname, "w") do fout
                global segment_ranges, curve_ranges, nc
                return curve_ranges, segment_ranges, nc = DT.write_curves!(fout, x, y)
            end
            fname2 = tempname()
            open(fname2, "w") do fout
                global seg1, seg2, seg3, seg4, seg5
                init_point, init_line, seg1 = DT.write_segments!(fout, x1, y1, 1, 1)
                init_point, init_line, seg2 = DT.write_segments!(fout, x2, y2, init_point + 1,
                    init_line)
                init_point, init_line, seg3 = DT.write_segments!(fout, x3, y3, init_point + 1,
                    init_line)
                init_point, init_line, seg4 = DT.write_segments!(fout, x4, y4, init_point + 1,
                    init_line)
                return init_point, init_line, seg5 = DT.write_segments!(fout, x5, y5, init_point + 1,
                    init_line)
            end
            @test read(fname, String) == read(fname2, String)
            @test curve_ranges == [1:4, 5:10, 11:26, 27:38, 39:43]
            @test segment_ranges == [seg1, seg2, seg3, seg4, seg5]
            @test nc == 43
        end
    end
end

@testset "Writing curves" begin
    @testset "Writing a curve loop" begin
        fname = tempname()
        open(fname, "w") do fout
            return DT.write_curve_loop!(fout, 1, 1:7)
        end
        @test read(fname, String) == "Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7};\n"

        fname = tempname()
        open(fname, "w") do fout
            return DT.write_curve_loop!(fout, 4, 13:21)
        end
        @test read(fname, String) == "Curve Loop(4) = {13, 14, 15, 16, 17, 18, 19, 20, 21};\n"
    end

    @testset "Writing multiple curve loops" begin
        curve_range = [1:4, 5:10, 11:15]
        fname = tempname()
        open(fname, "w") do fout
            return DT.write_curve_loops!(fout, curve_range)
        end
@test read(fname, String) == "Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8, 9, 10};
Curve Loop(3) = {11, 12, 13, 14, 15};
"

        curve_range = [1:4]
        open(fname, "w") do fout
            return DT.write_curve_loops!(fout, curve_range)
        end
        @test read(fname, String) == "Curve Loop(1) = {1, 2, 3, 4};\n"
    end
end

@testset "Writing a plane surface" begin
    fname = tempname()
    open(fname, "w") do fout
        return DT.write_plane_surface!(fout, 1)
    end
    @test read(fname, String) == "Plane Surface(1) = {1};\n"

    fname = tempname()
    open(fname, "w") do fout
        return DT.write_plane_surface!(fout, 4)
    end
    @test read(fname, String) == "Plane Surface(1) = {1, 2, 3, 4};\n"
end

@testset "Defining physical features" begin
    @testset "Writing a physical curve" begin
        fname = tempname()
        open(fname, "w") do fout
            return DT.write_physical_curve!(fout, 2, 5:7)
        end
        @test read(fname, String) == "Physical Curve(2) = {5, 6, 7};\n"

        fname = tempname()
        open(fname, "w") do fout
            return DT.write_physical_curve!(fout, 1, 5)
        end
        @test read(fname, String) == "Physical Curve(1) = {5};\n"
    end

    @testset "Writing multiple physical curves" begin
        fname = tempname()
        rng = [1:4, 5:9, 10:11]
        open(fname, "w") do fout
            return DT.write_physical_curves!(fout, rng, 1)
        end
@test read(fname, String) == "Physical Curve(1) = {1, 2, 3, 4};
Physical Curve(2) = {5, 6, 7, 8, 9};
Physical Curve(3) = {10, 11};
"

        fname = tempname()
        rng = [1:4, 5:10, 11:20]
        open(fname, "w") do fout
            return DT.write_physical_curves!(fout, rng, 4)
        end
@test read(fname, String) == "Physical Curve(4) = {1, 2, 3, 4};
Physical Curve(5) = {5, 6, 7, 8, 9, 10};
Physical Curve(6) = {11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
"

        fname = tempname()
        rng = [[1:4], [5:10]]
        open(fname, "w") do fout
            return DT.write_physical_curves!(fout, rng, 11)
        end
@test read(fname, String) == "Physical Curve(11) = {1, 2, 3, 4};
Physical Curve(12) = {5, 6, 7, 8, 9, 10};
"
    end

    @testset "Writing a physical surface" begin
        fname = tempname()
        open(fname, "w") do fout
            return DT.write_physical_surface!(fout, 1, 1)
        end
        @test read(fname, String) == "Physical Surface(1) = {1};\n"
    end
end

@testset "Reading lines" begin
    @testset "Reading the node lines" begin
        tline = "2 0.7071067811865476 0.7071067811865475 0 r"
        xx, yy = DT.read_node_line(tline)
        @test xx == 0.7071067811865476
        @test yy == 0.7071067811865475
    end

    @testset "Reading the element lines" begin
        tline = "8 1 2 13 4 16 5"
        elm_type, physical_id, u, v, w = DT.read_element_line(tline)
        @test elm_type == 1
        @test physical_id == 13
        @test u == 4
        @test v == 16
        @test w == 5

        tline = "51 2 2 1 1 32 12 36"
        elm_type, physical_id, u, v, w = DT.read_element_line(tline)
        @test elm_type == 2
        @test physical_id == 1
        @test u == 32
        @test v == 12
        @test w == 36
    end

    @testset "Reading the number of objects" begin
        tline = "72"
        @test DT.read_num_objects(tline) == 72
    end

    @testset "Reading the mesh format line" begin
        tline = "2.2 0 8"
        @test DT.read_mesh_format(tline) == (2.2, 0, 8)
    end
end

if !(get(ENV, "CI", "false") == "true")
    @testset "Meshing" begin
        @testset "Writing an example mesh" begin
            θ1 = LinRange(0, π, 5)
            θ2 = LinRange(π, 2π, 5)
            ϕ = LinRange(2π, 0, 5)
            R = 1.0
            r = 0.5
            x1 = @. R * cos(θ1)
            y1 = @. R * sin(θ1)
            x2 = @. R * cos(θ2)
            y2 = @. R * sin(θ2)
            x3 = @. r * cos(ϕ)
            y3 = @. r * sin(ϕ)
            x = [[x1, x2], [x3]]
            y = [[y1, y2], [y3]]

            id, id_dict = DT.write_gmsh(x, y, 0.5; mesh_algorithm=6, verbosity=0)

            @test id == [[1, 2], [3]]
            @test id_dict == Dict(13 => (1, 1), 14 => (1, 2), 15 => (2, 1))
@test read("meshgeometry.geo", String) == "r = 0.5;
Mesh.Algorithm = 6;
Mesh.Format = 1;
General.Verbosity = 0;
Point(1) = {1.0, 0.0, 0, r};
Point(2) = {0.7071067811865476, 0.7071067811865475, 0, r};
Point(3) = {6.123233995736766e-17, 1.0, 0, r};
Point(4) = {-0.7071067811865475, 0.7071067811865476, 0, r};
Point(5) = {-1.0, 1.2246467991473532e-16, 0, r};
Point(6) = {-1.0, 1.2246467991473532e-16, 0, r};
Point(7) = {-0.7071067811865477, -0.7071067811865475, 0, r};
Point(8) = {-1.8369701987210297e-16, -1.0, 0, r};
Point(9) = {0.7071067811865474, -0.7071067811865477, 0, r};
Point(10) = {1.0, -2.4492935982947064e-16, 0, r};
Point(11) = {0.5, -1.2246467991473532e-16, 0, r};
Point(12) = {-9.184850993605148e-17, -0.5, 0, r};
Point(13) = {-0.5, 6.123233995736766e-17, 0, r};
Point(14) = {3.061616997868383e-17, 0.5, 0, r};
Point(15) = {0.5, 0.0, 0, r};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 7};
Line(6) = {7, 8};
Line(7) = {8, 9};
Line(8) = {9, 1};
Line(9) = {11, 12};
Line(10) = {12, 13};
Line(11) = {13, 14};
Line(12) = {14, 11};
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Curve Loop(2) = {9, 10, 11, 12};
Plane Surface(1) = {1, 2};
Physical Curve(13) = {1, 2, 3, 4};
Physical Curve(14) = {5, 6, 7, 8};
Physical Curve(15) = {9, 10, 11, 12};
Physical Surface(1) = {1};
"

            DT.run_gmsh("./gmsh-4.11.1-Windows64/gmsh.exe")

@test read("meshgeometry.msh", String) == "\$MeshFormat\r
2.2 0 8\r
\$EndMeshFormat\r
\$Nodes\r
36\r
1 1 0 0\r
2 0.7071067811865476 0.7071067811865475 0\r
3 6.123233995736766e-17 1 0\r
4 -0.7071067811865475 0.7071067811865476 0\r
5 -1 1.224646799147353e-16 0\r
6 -0.7071067811865477 -0.7071067811865475 0\r
7 -1.83697019872103e-16 -1 0\r
8 0.7071067811865474 -0.7071067811865477 0\r
9 0.5 -1.224646799147353e-16 0\r
10 -9.184850993605148e-17 -0.5 0\r
11 -0.5 6.123233995736766e-17 0\r
12 3.061616997868383e-17 0.5 0\r
13 0.853553390593673 0.35355339059231 0\r
14 0.3535533905941049 0.8535533905929296 0\r
15 -0.3535533905924035 0.8535533905936342 0\r
16 -0.8535533905929296 0.3535533905941049 0\r
17 -0.8535533905936216 -0.3535533905924342 0\r
18 -0.3535533905942161 -0.8535533905928835 0\r
19 0.3535533905924341 -0.8535533905936216 0\r
20 0.8535533905928643 -0.3535533905942622 0\r
21 0.2500000000005643 -0.2499999999994356 0\r
22 -0.2499999999994354 -0.2500000000005645 0\r
23 -0.2500000000008388 0.2499999999991612 0\r
24 0.2499999999994376 0.2500000000005623 0\r
25 -0.3585358256888765 0.5492598171197338 0\r
26 0.5492598171197998 0.3585358256888155 0\r
27 0.3585358256888558 -0.5492598171197742 0\r
28 -0.5492598171197743 -0.3585358256888564 0\r
29 -0.5732233047030799 0.2803300858900932 0\r
30 0.5732233047030577 -0.2803300858900962 0\r
31 -0.2803300858900968 -0.5732233047030627 0\r
32 0.2803300858901325 0.5732233047030406 0\r
33 -0.7355093454939969 -0.05737017598737389 0\r
34 0.7355093454939962 0.05737017598731059 0\r
35 -0.05737017598736074 0.7355093454939913 0\r
36 0.05737017598736395 -0.7355093454939916 0\r
\$EndNodes\r
\$Elements\r
72\r
1 1 2 13 1 1 13\r
2 1 2 13 1 13 2\r
3 1 2 13 2 2 14\r
4 1 2 13 2 14 3\r
5 1 2 13 3 3 15\r
6 1 2 13 3 15 4\r
7 1 2 13 4 4 16\r
8 1 2 13 4 16 5\r
9 1 2 14 5 5 17\r
10 1 2 14 5 17 6\r
11 1 2 14 6 6 18\r
12 1 2 14 6 18 7\r
13 1 2 14 7 7 19\r
14 1 2 14 7 19 8\r
15 1 2 14 8 8 20\r
16 1 2 14 8 20 1\r
17 1 2 15 9 9 21\r
18 1 2 15 9 21 10\r
19 1 2 15 10 10 22\r
20 1 2 15 10 22 11\r
21 1 2 15 11 11 23\r
22 1 2 15 11 23 12\r
23 1 2 15 12 12 24\r
24 1 2 15 12 24 9\r
25 2 2 1 1 30 20 34\r
26 2 2 1 1 29 16 33\r
27 2 2 1 1 31 18 36\r
28 2 2 1 1 32 14 35\r
29 2 2 1 1 20 1 34\r
30 2 2 1 1 18 7 36\r
31 2 2 1 1 16 5 33\r
32 2 2 1 1 14 3 35\r
33 2 2 1 1 23 12 25\r
34 2 2 1 1 22 11 28\r
35 2 2 1 1 21 10 27\r
36 2 2 1 1 24 9 26\r
37 2 2 1 1 25 4 29\r
38 2 2 1 1 27 8 30\r
39 2 2 1 1 28 6 31\r
40 2 2 1 1 26 2 32\r
41 2 2 1 1 28 11 33\r
42 2 2 1 1 25 12 35\r
43 2 2 1 1 27 10 36\r
44 2 2 1 1 26 9 34\r
45 2 2 1 1 8 20 30\r
46 2 2 1 1 6 18 31\r
47 2 2 1 1 4 16 29\r
48 2 2 1 1 2 14 32\r
49 2 2 1 1 11 29 33\r
50 2 2 1 1 12 32 35\r
51 2 2 1 1 10 31 36\r
52 2 2 1 1 9 30 34\r
53 2 2 1 1 13 2 26\r
54 2 2 1 1 15 4 25\r
55 2 2 1 1 19 8 27\r
56 2 2 1 1 17 6 28\r
57 2 2 1 1 15 25 35\r
58 2 2 1 1 13 26 34\r
59 2 2 1 1 17 28 33\r
60 2 2 1 1 19 27 36\r
61 2 2 1 1 11 23 29\r
62 2 2 1 1 9 21 30\r
63 2 2 1 1 10 22 31\r
64 2 2 1 1 12 24 32\r
65 2 2 1 1 5 17 33\r
66 2 2 1 1 7 19 36\r
67 2 2 1 1 3 15 35\r
68 2 2 1 1 1 13 34\r
69 2 2 1 1 23 25 29\r
70 2 2 1 1 21 27 30\r
71 2 2 1 1 22 28 31\r
72 2 2 1 1 24 26 32\r
\$EndElements\r
"

            elements, nodes, boundary_nodes = DT.read_gmsh(id, id_dict)

            @test nodes ≈
                  [1 0
                0.7071067811865476 0.7071067811865475
                6.123233995736766e-17 1
                -0.7071067811865475 0.7071067811865476
                -1 1.224646799147353e-16
                -0.7071067811865477 -0.7071067811865475
                -1.83697019872103e-16 -1
                0.7071067811865474 -0.7071067811865477
                0.5 -1.224646799147353e-16
                -9.184850993605148e-17 -0.5
                -0.5 6.123233995736766e-17
                3.061616997868383e-17 0.5
                0.853553390593673 0.35355339059231
                0.3535533905941049 0.8535533905929296
                -0.3535533905924035 0.8535533905936342
                -0.8535533905929296 0.3535533905941049
                -0.8535533905936216 -0.3535533905924342
                -0.3535533905942161 -0.8535533905928835
                0.3535533905924341 -0.8535533905936216
                0.8535533905928643 -0.3535533905942622
                0.2500000000005643 -0.2499999999994356
                -0.2499999999994354 -0.2500000000005645
                -0.2500000000008388 0.2499999999991612
                0.2499999999994376 0.2500000000005623
                -0.3585358256888765 0.5492598171197338
                0.5492598171197998 0.3585358256888155
                0.3585358256888558 -0.5492598171197742
                -0.5492598171197743 -0.3585358256888564
                -0.5732233047030799 0.2803300858900932
                0.5732233047030577 -0.2803300858900962
                -0.2803300858900968 -0.5732233047030627
                0.2803300858901325 0.5732233047030406
                -0.7355093454939969 -0.05737017598737389
                0.7355093454939962 0.05737017598731059
                -0.05737017598736074 0.7355093454939913
                0.05737017598736395 -0.7355093454939916]'

            @test elements ==
                  [30 20 34
                29 16 33
                31 18 36
                32 14 35
                20 1 34
                18 7 36
                16 5 33
                14 3 35
                23 12 25
                22 11 28
                21 10 27
                24 9 26
                25 4 29
                27 8 30
                28 6 31
                26 2 32
                28 11 33
                25 12 35
                27 10 36
                26 9 34
                8 20 30
                6 18 31
                4 16 29
                2 14 32
                11 29 33
                12 32 35
                10 31 36
                9 30 34
                13 2 26
                15 4 25
                19 8 27
                17 6 28
                15 25 35
                13 26 34
                17 28 33
                19 27 36
                11 23 29
                9 21 30
                10 22 31
                12 24 32
                5 17 33
                7 19 36
                3 15 35
                1 13 34
                23 25 29
                21 27 30
                22 28 31
                24 26 32]'
            bn1 = [1, 13, 2, 14, 3, 15, 4, 16, 5]
            bn2 = [5, 17, 6, 18, 7, 19, 8, 20, 1]
            bn3 = [9, 21, 10, 22, 11, 23, 12, 24, 9]
            @test boundary_nodes == [[bn1, bn2], [bn3]]
        end

        @testset "Meshing a single curve with multiple segments" begin
            θ1 = LinRange(0, π, 5)
            θ2 = LinRange(π, 2π, 5)
            ϕ = LinRange(2π, 0, 5)
            R = 1.0
            r = 0.5
            x1 = @. R * cos(θ1)
            y1 = @. R * sin(θ1)
            x2 = @. R * cos(θ2)
            y2 = @. R * sin(θ2)
            x3 = @. r * cos(ϕ)
            y3 = @. r * sin(ϕ)
            x = [[x1, x2], [x3]]
            y = [[y1, y2], [y3]]
            x = [x1, x2]
            y = [y1, y2]
            elements, nodes, boundary_nodes = generate_mesh(x, y, 0.5; convert_result=false)
            g1 = read("meshgeometry.geo", String)
            m1 = read("meshgeometry.msh", String)
            _elements, _nodes, _boundary_nodes = generate_mesh([x], [y], 0.5; convert_result=false)
            g2 = read("meshgeometry.geo", String)
            m2 = read("meshgeometry.msh", String)
            @test g1 == g2
            @test m1 == m2
            @test boundary_nodes == _boundary_nodes[1]
            @test elements == _elements
            @test nodes == _nodes
        end
        @testset "Meshing a single curve" begin
            θ1 = LinRange(0, π, 5)
            θ2 = LinRange(π, 2π, 5)
            ϕ = LinRange(2π, 0, 5)
            R = 1.0
            r = 0.5
            x1 = @. R * cos(θ1)
            y1 = @. R * sin(θ1)
            x2 = @. R * cos(θ2)
            y2 = @. R * sin(θ2)
            x3 = @. r * cos(ϕ)
            y3 = @. r * sin(ϕ)
            x = [[x1, x2], [x3]]
            y = [[y1, y2], [y3]]
            x = [x1..., x2...]
            y = [y1..., y2...]
            elements, nodes, boundary_nodes = generate_mesh(x, y, 0.5; convert_result=false, check_arguments=false)
            g1 = read("meshgeometry.geo", String)
            m1 = read("meshgeometry.msh", String)
            _elements, _nodes, _boundary_nodes = generate_mesh([[x]], [[y]], 0.5;
                convert_result=false, check_arguments=false)
            g2 = read("meshgeometry.geo", String)
            m2 = read("meshgeometry.msh", String)
            @test g1 == g2
            @test m1 == m2
            @test boundary_nodes == _boundary_nodes[1][1]
            @test elements == _elements
            @test nodes == _nodes
        end
        @testset "Meshing a square" begin
            x = [0.0, 2.0, 2.0, 0.0, 0.0]
            y = [0.0, 0.0, 2.0, 2.0, 0.0]
            elements1, nodes1, boundary_nodes1 = generate_mesh(x, y, 0.5; convert_result=false)
            fig = Figure()
            ax = Axis(fig[1, 1])
            mesh!(ax, nodes1, elements1')
            x = [[0.0, 1.0, 2.0], [2.0, 2.0, 2.0], [2.0, 1.0, 0.0], [0.0, 0.0, 0.0]]
            y = [[0.0, 0.0, 0.0], [0.0, 1.0, 2.0], [2.0, 2.0, 2.0], [2.0, 1.0, 0.0]]
            elements2, nodes2, boundary_nodes2 = generate_mesh(x, y, 0.5; convert_result=false)
            ax = Axis(fig[1, 2])
            mesh!(ax, nodes2, elements2')
            elements3, nodes3, boundary_nodes3 = generate_mesh(0.0, 2.0, 0.0, 2.0, 0.5;
                convert_result=false)
            ax = Axis(fig[2, 1])
            mesh!(ax, nodes3, elements3')
            elements4, nodes4, boundary_nodes4 = generate_mesh(0.0, 2.0, 0.0, 2.0, 0.5;
                single_boundary=false,
                convert_result=false)
            ax = Axis(fig[2, 2])
            mesh!(ax, nodes4, elements4')
            @test elements1 == elements3
            @test nodes1 == nodes3
            @test boundary_nodes1 == boundary_nodes1
            @test elements2 == elements4
            @test nodes2 == nodes4
            @test boundary_nodes2 == boundary_nodes4
        end
    end
end

if !(get(ENV, "CI", "false") == "true")
    @testset "Triangulating" begin
        @testset "A square: One boundary with no ghost triangles" begin
            tri = generate_mesh(0.0, 2.0, 0.0, 2.0, 0.5; add_ghost_triangles=false)
            elements, nodes, bn = generate_mesh(0.0, 2.0, 0.0, 2.0, 0.5; convert_result=false)
            @test tri.points == nodes
            @test tri.triangles == Set{NTuple{3,Int}}((T[1], T[2], T[3]) for T in eachcol(elements))
            adj = tri.adjacent.adjacent
            adj = DefaultDict(DT.DefaultAdjacentValue, adj)
            for T in tri.triangles
                i, j, k = T
                @test adj[(i, j)] == k
                @test adj[(j, k)] == i
                @test adj[(k, i)] == j
            end
            adj2v = tri.adjacent2vertex.adjacent2vertex
            for (i, S) in adj2v
                for (j, k) in S
                    if all(≥(DT.FirstPointIndex), (i, j, k))
                        @test adj[(j, k)] == i
                        @test adj[(k, i)] == j
                        @test adj[(i, j)] == k
                    end
                end
            end
            for (j, k) in adj2v[DT.BoundaryIndex]
                @test adj[(j, k)] == DT.BoundaryIndex
                @test j ∈ bn && k ∈ bn
            end
            @test tri.boundary_nodes == bn
            g = tri.graph
            @test sort(collect(g.graph[DT.BoundaryIndex])) == sort(unique(bn))
            DT.delete_vertex!(g, DT.BoundaryIndex)
            for (i, j) in g.graph.E
                k = adj[(i, j)]
                el = (i, j, k)
                e = (i, j)
                if k < DT.FirstPointIndex
                    k = adj[(j, i)]
                    el = (k, j, i)
                    e = (j, i)
                end
                τ, pred = DT.contains_triangle(el, tri.triangles)
                @test DT.compare_triangles(τ, el)
                @test pred
                @test e ∈ adj2v[k]
            end
            @test tri.constrained_edges == Set{NTuple{2,Int}}()
            @test tri.boundary_map == OrderedDict(DT.BoundaryIndex => bn)
            for i in 1:(lastindex(bn)-1)
                local u, v, w
                u = bn[i]
                v = bn[i+1]
                w = adj[(v, u)]
                @test adj[(v, u)] == DT.BoundaryIndex
                @test adj[(u, w)] == DT.DefaultAdjacentValue
                @test adj[(w, v)] == DT.DefaultAdjacentValue
                @test (v, u) ∈ adj2v[DT.BoundaryIndex]
            end
            for (boundary_index, segment_index) in tri.boundary_map
                for (j, k) in adj2v[boundary_index]
                    @test adj[(j, k)] == boundary_index
                    @test j ∈ DT.get_boundary_nodes(bn, segment_index) &&
                        k ∈ DT.get_boundary_nodes(bn, segment_index)
                end
            end
        end

        @testset "A square: One boundary with ghost triangles" begin
            tri = generate_mesh(0.0, 2.0, 0.0, 2.0, 0.5; add_ghost_triangles=true)
            elements, nodes, bn = generate_mesh(0.0, 2.0, 0.0, 2.0, 0.5; convert_result=false)
            @test tri.points == nodes
            elset = Set{NTuple{3,Int}}((T[1], T[2], T[3]) for T in eachcol(elements))
            ng = 0
            for T in tri.triangles
                i, j, k = indices(T)
                if all(≥(DT.FirstPointIndex), (i, j, k))
                    @test T ∈ elset
                elseif i < DT.FirstPointIndex
                    @test j ∈ bn && k ∈ bn
                    @test DefaultDict(DT.DefaultAdjacentValue, tri.adjacent.adjacent)[(j, k)] == i
                    @test j ∈ tri.graph.graph[DT.BoundaryIndex] && k ∈ tri.graph.graph[DT.BoundaryIndex]
                    @test (j, k) ∈ tri.adjacent2vertex.adjacent2vertex[DT.BoundaryIndex]
                    ng += 1
                elseif j < DT.FirstPointIndex
                    @test i ∈ bn && k ∈ bn
                    @test DefaultDict(DT.DefaultAdjacentValue, tri.adjacent.adjacent)[(k, i)] == j
                    @test i ∈ tri.graph.graph[DT.BoundaryIndex] && k ∈ tri.graph.graph[DT.BoundaryIndex]
                    @test (k, i) ∈ tri.adjacent2vertex.adjacent2vertex[DT.BoundaryIndex]
                    ng += 1
                elseif k < DT.FirstPointIndex
                    @test i ∈ bn && j ∈ bn
                    @test DefaultDict(DT.DefaultAdjacentValue, tri.adjacent.adjacent)[(i, j)] == k
                    @test i ∈ tri.graph.graph[DT.BoundaryIndex] && j ∈ tri.graph.graph[DT.BoundaryIndex]
                    @test (i, j) ∈ tri.adjacent2vertex.adjacent2vertex[DT.BoundaryIndex]
                    ng += 1
                end
            end
            @test ng == DT.num_boundary_edges(bn)
            adj = tri.adjacent.adjacent
            adj = DefaultDict(DT.DefaultAdjacentValue, adj)
            for T in tri.triangles
                i, j, k = T
                @test adj[(i, j)] == k
                @test adj[(j, k)] == i
                @test adj[(k, i)] == j
            end
            adj2v = tri.adjacent2vertex.adjacent2vertex
            for ((u, v), w) in adj
                @test (u, v) ∈ adj2v[w]
            end
            for (i, S) in adj2v
                for (j, k) in S
                    @test adj[(j, k)] == i
                    @test adj[(k, i)] == j
                    @test adj[(i, j)] == k
                end
            end
            for (j, k) in adj2v[DT.BoundaryIndex]
                @test adj[(j, k)] == DT.BoundaryIndex
                @test j ∈ bn && k ∈ bn
            end
            @test tri.boundary_nodes == bn
            g = tri.graph
            for (i, j) in g.graph.E
                k = adj[(i, j)]
                el = (i, j, k)
                τ, pred = DT.contains_triangle(el, tri.triangles)
                @test DT.compare_triangles(τ, el)
                @test (i, j) ∈ adj2v[k]
            end
            @test tri.constrained_edges == Set{NTuple{2,Int}}()
            @test tri.boundary_map == OrderedDict(DT.BoundaryIndex => bn)
            for i in 1:(lastindex(bn)-1)
                local u, v, w
                u = bn[i]
                v = bn[i+1]
                w = adj[(v, u)]
                @test adj[(v, u)] == DT.BoundaryIndex
                @test adj[(u, w)] == v
                @test adj[(w, v)] == u
                @test (v, u) ∈ adj2v[DT.BoundaryIndex]
                @test u ∈ tri.graph.graph[DT.BoundaryIndex] && v ∈ tri.graph.graph[DT.BoundaryIndex]
            end
            @test sort(collect(g.graph[DT.BoundaryIndex])) == sort(unique(bn))
            for (boundary_index, segment_index) in tri.boundary_map
                for (j, k) in adj2v[boundary_index]
                    @test adj[(j, k)] == boundary_index
                    @test j ∈ DT.get_boundary_nodes(bn, segment_index) &&
                        k ∈ DT.get_boundary_nodes(bn, segment_index)
                end
            end
            for i in length(tri.boundary_nodes):-1:2
                local u, v
                u = tri.boundary_nodes[i]
                v = tri.boundary_nodes[i-1]
                @test DT.get_adjacent(tri, u, v) == DT.BoundaryIndex
            end
        end

        @testset "A square: Multiple boundaries with no ghost triangles" begin
            tri = generate_mesh(0.0, 2.0, 0.0, 2.0, 0.5; add_ghost_triangles=false,
                single_boundary=false)
            elements, nodes, bn = generate_mesh(0.0, 2.0, 0.0, 2.0, 0.5; convert_result=false,
                single_boundary=false)
            @test tri.points == nodes
            @test tri.triangles == Set{NTuple{3,Int}}((T[1], T[2], T[3]) for T in eachcol(elements))
            @test tri.boundary_nodes == bn
            @test tri.boundary_map == OrderedDict(-(1:4) .=> 1:4)
            @test collect(Base.keys(tri.boundary_map)) == [-1, -2, -3, -4]
            adj = tri.adjacent.adjacent
            adj = DefaultDict(DT.DefaultAdjacentValue, adj)
            for T in tri.triangles
                i, j, k = T
                @test adj[(i, j)] == k
                @test adj[(j, k)] == i
                @test adj[(k, i)] == j
            end
            adj2v = tri.adjacent2vertex.adjacent2vertex
            for ((u, v), w) in adj
                @test (u, v) ∈ adj2v[w]
            end
            for (i, S) in adj2v
                for (j, k) in S
                    if all(≥(DT.FirstPointIndex), (i, j, k))
                        @test adj[(j, k)] == i
                        @test adj[(k, i)] == j
                        @test adj[(i, j)] == k
                    end
                end
            end
            for (boundary_index, segment_index) in tri.boundary_map
                for (j, k) in adj2v[boundary_index]
                    @test adj[(j, k)] == boundary_index
                    @test j ∈ bn[segment_index] && k ∈ bn[segment_index]
                end
            end
            @test tri.boundary_nodes == bn
            g = tri.graph
            @test sort(collect(g.graph[DT.BoundaryIndex])) == sort(unique(bn[1]))
            @test sort(collect(g.graph[DT.BoundaryIndex-1])) == sort(unique(bn[2]))
            @test sort(collect(g.graph[DT.BoundaryIndex-2])) == sort(unique(bn[3]))
            @test sort(collect(g.graph[DT.BoundaryIndex-3])) == sort(unique(bn[4]))
            for j in 1:4
                for i in 1:(lastindex(bn[j])-1)
                    local u, v, w
                    u = bn[j][i]
                    v = bn[j][i+1]
                    w = adj[(v, u)]
                    @test adj[(v, u)] == DT.BoundaryIndex - j + 1
                    @test adj[(u, w)] == DT.DefaultAdjacentValue
                    @test adj[(w, v)] == DT.DefaultAdjacentValue
                    @test (v, u) ∈ adj2v[DT.BoundaryIndex-j+1]
                    @test u ∈ tri.graph.graph[DT.BoundaryIndex-j+1] &&
                        v ∈ tri.graph.graph[DT.BoundaryIndex-j+1]
                end
            end
            DT.delete_boundary_vertices_from_graph!(g)
            for (i, j) in g.graph.E
                k = adj[(i, j)]
                el = (i, j, k)
                e = (i, j)
                if k < DT.FirstPointIndex
                    k = adj[(j, i)]
                    el = (k, j, i)
                    e = (j, i)
                end
                τ, pred = DT.contains_triangle(el, tri.triangles)
                @test DT.compare_triangles(τ, el)
                @test pred
                @test e ∈ adj2v[k]
            end
            @test tri.constrained_edges == Set{NTuple{2,Int}}()
        end

        @testset "A square: Multiple boundaries with ghost triangles" begin
            tri = generate_mesh(0.0, 2.0, 0.0, 2.0, 0.5; add_ghost_triangles=true,
                single_boundary=false)
            elements, nodes, bn = generate_mesh(0.0, 2.0, 0.0, 2.0, 0.5; convert_result=false,
                single_boundary=false)
            @test tri.points == nodes
            elset = Set{NTuple{3,Int}}((T[1], T[2], T[3]) for T in eachcol(elements))
            ng = 0
            for T in tri.triangles
                i, j, k = indices(T)
                if all(≥(DT.FirstPointIndex), (i, j, k))
                    @test T ∈ elset
                elseif i < DT.FirstPointIndex
                    n = tri.boundary_map[i]
                    @test j ∈ bn[n] && k ∈ bn[n]
                    @test DefaultDict(DT.DefaultAdjacentValue, tri.adjacent.adjacent)[(j, k)] == i
                    @test j ∈ tri.graph.graph[i] && k ∈ tri.graph.graph[i]
                    @test (j, k) ∈ tri.adjacent2vertex.adjacent2vertex[i]
                    ng += 1
                elseif j < DT.FirstPointIndex
                    n = tri.boundary_map[j]
                    @test k ∈ bn[n] && i ∈ bn[n]
                    @test DefaultDict(DT.DefaultAdjacentValue, tri.adjacent.adjacent)[(k, i)] == j
                    @test k ∈ tri.graph.graph[j] && i ∈ tri.graph.graph[j]
                    @test (k, i) ∈ tri.adjacent2vertex.adjacent2vertex[j]
                    ng += 1
                elseif k < DT.FirstPointIndex
                    n = tri.boundary_map[k]
                    @test i ∈ bn[n] && j ∈ bn[n]
                    @test DefaultDict(DT.DefaultAdjacentValue, tri.adjacent.adjacent)[(i, j)] == k
                    @test i ∈ tri.graph.graph[k] && j ∈ tri.graph.graph[k]
                    @test (i, j) ∈ tri.adjacent2vertex.adjacent2vertex[k]
                    ng += 1
                end
            end
            @test tri.boundary_nodes == bn
            @test tri.boundary_map == OrderedDict(-(1:4) .=> 1:4)
            @test collect(Base.keys(tri.boundary_map)) == [-1, -2, -3, -4]
            adj = tri.adjacent.adjacent
            adj = DefaultDict(DT.DefaultAdjacentValue, adj)
            for T in tri.triangles
                i, j, k = T
                @test adj[(i, j)] == k
                @test adj[(j, k)] == i
                @test adj[(k, i)] == j
            end
            adj2v = tri.adjacent2vertex.adjacent2vertex
            for (i, S) in adj2v
                for (j, k) in S
                    @test adj[(j, k)] == i
                    @test adj[(k, i)] == j
                    @test adj[(i, j)] == k
                end
            end
            for (boundary_index, segment_index) in tri.boundary_map
                for (j, k) in adj2v[boundary_index]
                    @test adj[(j, k)] == boundary_index
                    @test adj[(k, boundary_index)] == j
                    @test adj[(boundary_index, j)] == k
                    @test j ∈ bn[segment_index] && k ∈ bn[segment_index]
                end
            end
            @test tri.boundary_nodes == bn
            for j in 1:4
                for i in 1:(lastindex(bn[j])-1)
                    local u, v, w
                    u = bn[j][i]
                    v = bn[j][i+1]
                    w = adj[(v, u)]
                    @test adj[(v, u)] == DT.BoundaryIndex - j + 1
                    @test adj[(u, w)] == v
                    @test adj[(w, v)] == u
                    @test (v, u) ∈ adj2v[DT.BoundaryIndex-j+1]
                    @test u ∈ tri.graph.graph[DT.BoundaryIndex-j+1] &&
                        v ∈ tri.graph.graph[DT.BoundaryIndex-j+1]
                end
            end
            g = tri.graph
            for (i, j) in g.graph.E
                k = adj[(i, j)]
                if k == DT.DefaultAdjacentValue # since the graph is undirected, the ghost edges could be in the wrong direction
                    @test i < DT.FirstPointIndex || j < DT.FirstPointIndex
                    i, j = j, i
                    k = adj[(i, j)]
                end
                el = (i, j, k)
                τ, pred = DT.contains_triangle(el, tri.triangles)
                @test DT.compare_triangles(τ, el)
                @test (i, j) ∈ adj2v[k]
            end
            @test tri.constrained_edges == Set{NTuple{2,Int}}()
            @test DT.get_adjacent(tri, 1, 16) == DT.BoundaryIndex - 3
            @test DT.get_adjacent(tri, 16, 8) == DT.BoundaryIndex - 3
            @test DT.get_adjacent(tri, 8, 15) == DT.BoundaryIndex - 3
            @test DT.get_adjacent(tri, 15, 7) == DT.BoundaryIndex - 3
            @test DT.get_adjacent(tri, 7, 14) == DT.BoundaryIndex - 2
            @test DT.get_adjacent(tri, 14, 6) == DT.BoundaryIndex - 2
            @test DT.get_adjacent(tri, 6, 13) == DT.BoundaryIndex - 2
            @test DT.get_adjacent(tri, 13, 5) == DT.BoundaryIndex - 2
            @test DT.get_adjacent(tri, 5, 12) == DT.BoundaryIndex - 1
            @test DT.get_adjacent(tri, 12, 4) == DT.BoundaryIndex - 1
            @test DT.get_adjacent(tri, 4, 11) == DT.BoundaryIndex - 1
            @test DT.get_adjacent(tri, 11, 3) == DT.BoundaryIndex - 1
            @test DT.get_adjacent(tri, 3, 10) == DT.BoundaryIndex
            @test DT.get_adjacent(tri, 10, 2) == DT.BoundaryIndex
            @test DT.get_adjacent(tri, 2, 9) == DT.BoundaryIndex
            @test DT.get_adjacent(tri, 9, 1) == DT.BoundaryIndex
            @test sum(length(tri.adjacent2vertex.adjacent2vertex[DT.BoundaryIndex])) +
                sum(length(tri.adjacent2vertex.adjacent2vertex[DT.BoundaryIndex-1])) +
                sum(length(tri.adjacent2vertex.adjacent2vertex[DT.BoundaryIndex-2])) +
                sum(length(tri.adjacent2vertex.adjacent2vertex[DT.BoundaryIndex-3])) == 16
        end

        @testset "A complicated geometry" begin
            x, y = complicated_geometry()
            tri = generate_mesh(x, y, 0.1; convert_result=true, add_ghost_triangles=true)

            elements, nodes, bn = generate_mesh(x, y, 0.1; convert_result=false)
            @test tri.points == nodes
            elset = Set{NTuple{3,Int}}((T[1], T[2], T[3]) for T in eachcol(elements))
            ng = 0
            for T in tri.triangles
                i, j, k = indices(T)
                if all(≥(DT.FirstPointIndex), (i, j, k))
                    @test T ∈ elset
                elseif i < DT.FirstPointIndex
                    m, n = tri.boundary_map[i]
                    @test j ∈ bn[m][n] && k ∈ bn[m][n]
                    @test DefaultDict(DT.DefaultAdjacentValue, tri.adjacent.adjacent)[(j, k)] == i
                    @test j ∈ tri.graph.graph[i] && k ∈ tri.graph.graph[i]
                    @test (j, k) ∈ tri.adjacent2vertex.adjacent2vertex[i]
                    ng += 1
                elseif j < DT.FirstPointIndex
                    m, n = tri.boundary_map[j]
                    @test k ∈ bn[m][n] && i ∈ bn[m][n]
                    @test DefaultDict(DT.DefaultAdjacentValue, tri.adjacent.adjacent)[(k, i)] == j
                    @test k ∈ tri.graph.graph[j] && i ∈ tri.graph.graph[j]
                    @test (k, i) ∈ tri.adjacent2vertex.adjacent2vertex[j]
                    ng += 1
                elseif k < DT.FirstPointIndex
                    m, n = tri.boundary_map[k]
                    @test i ∈ bn[m][n] && j ∈ bn[m][n]
                    @test DefaultDict(DT.DefaultAdjacentValue, tri.adjacent.adjacent)[(i, j)] == k
                    @test i ∈ tri.graph.graph[k] && j ∈ tri.graph.graph[k]
                    @test (i, j) ∈ tri.adjacent2vertex.adjacent2vertex[k]
                    ng += 1
                end
            end
            @test tri.boundary_nodes == bn
            @test tri.boundary_map ==
                OrderedDict(-1 => (1, 1), -2 => (1, 2), -3 => (1, 3), -4 => (1, 4),
                -5 => (2, 1),
                -6 => (3, 1),
                -7 => (4, 1), -8 => (4, 2), -9 => (4, 3), -10 => (4, 4),
                -11 => (5, 1))
            @test collect(Base.keys(tri.boundary_map)) == [-1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11]
            adj = tri.adjacent.adjacent
            adj = DefaultDict(DT.DefaultAdjacentValue, adj)
            for T in tri.triangles
                i, j, k = T
                @test adj[(i, j)] == k
                @test adj[(j, k)] == i
                @test adj[(k, i)] == j
            end
            adj2v = tri.adjacent2vertex.adjacent2vertex
            for (i, S) in adj2v
                for (j, k) in S
                    @test adj[(j, k)] == i
                    @test adj[(k, i)] == j
                    @test adj[(i, j)] == k
                end
            end
            for (boundary_index, segment_index) in tri.boundary_map
                for (j, k) in adj2v[boundary_index]
                    @test adj[(j, k)] == boundary_index
                    @test adj[(k, boundary_index)] == j
                    @test adj[(boundary_index, j)] == k
                    @test j ∈ DT.get_boundary_nodes(bn, segment_index) &&
                        k ∈ DT.get_boundary_nodes(bn, segment_index)
                end
            end
            @test tri.boundary_nodes == bn
            current_index = DT.BoundaryIndex
            for m in 1:DT.num_curves(bn)
                bn_m = DT.get_boundary_nodes(bn, m)
                for n in 1:DT.num_segments(bn_m)
                    bn_n = DT.get_boundary_nodes(bn_m, n)
                    for ℓ in 1:(lastindex(bn_m)-1)
                        local u, v, w
                        u = bn_n[ℓ]
                        v = bn_n[ℓ+1]
                        w = adj[(v, u)]
                        @test adj[(v, u)] == current_index
                        @test adj[(u, w)] == v
                        @test adj[(w, v)] == u
                        @test (v, u) ∈ adj2v[current_index]
                        @test u ∈ tri.graph.graph[current_index] && v ∈ tri.graph.graph[current_index]
                    end
                    current_index -= 1
                end
            end
            g = tri.graph
            for (i, j) in g.graph.E
                k = adj[(i, j)]
                if k == DT.DefaultAdjacentValue # since the graph is undirected, the ghost edges could be in the wrong direction
                    @test i < DT.FirstPointIndex || j < DT.FirstPointIndex
                    i, j = j, i
                    k = adj[(i, j)]
                end
                el = (i, j, k)
                τ, pred = DT.contains_triangle(el, tri.triangles)
                @test DT.compare_triangles(τ, el)
                @test (i, j) ∈ adj2v[k]
            end
            @test tri.constrained_edges == Set{NTuple{2,Int}}()

            bn1 = tri.boundary_nodes[1]
            for i in eachindex(bn1)
                for (u, v) in zip(bn1[i][end:-1:2], bn1[i][(end-1):-1:1])
                    @test DT.get_adjacent(tri, u, v) == DT.BoundaryIndex - i + 1
                end
            end
            bn2, bn3, bn4, bn5 = tri.boundary_nodes[(begin+1):end]
            for i in eachindex(bn2)
                for (u, v) in zip(bn2[i][end:-1:2], bn2[i][(end-1):-1:1])
                    @test DT.get_adjacent(tri, u, v) == DT.BoundaryIndex - 5 + 1
                end
            end
            for i in eachindex(bn3)
                for (u, v) in zip(bn3[i][end:-1:2], bn3[i][(end-1):-1:1])
                    @test DT.get_adjacent(tri, u, v) == DT.BoundaryIndex - 6 + 1
                end
            end
            for i in eachindex(bn4)
                for (u, v) in zip(bn4[i][end:-1:2], bn4[i][(end-1):-1:1])
                    @test DT.get_adjacent(tri, u, v) == (-7, -8, -9, -10)[i]
                end
            end
            for i in eachindex(bn5)
                for (u, v) in zip(bn5[i][end:-1:2], bn5[i][(end-1):-1:1])
                    @test DT.get_adjacent(tri, u, v) == -11
                end
            end
            @test DT.get_adjacent(tri, 111, 337) == -10
            @test DT.get_adjacent(tri, 337, 122) == -10
            @test DT.get_adjacent(tri, 122, 336) == -10
            @test DT.get_adjacent(tri, 336, 121) == -10
            @test DT.get_adjacent(tri, 121, 335) == -10
            @test DT.get_adjacent(tri, 335, 120) == -10
            @test DT.get_adjacent(tri, 120, 334) == -9
            @test DT.get_adjacent(tri, 334, 333) == -9
            @test DT.get_adjacent(tri, 117, 307) == -8
            @test tri.convex_hull.points == tri.points
            @test tri.convex_hull.indices == convex_hull(tri.points).indices
        end

        @testset "Another specific example: Tricuspoid" begin
            a = 4 / 5
            t = LinRange(0, 2π, 100)
            x = @. a * (2cos(t) + cos(2t))
            y = @. a * (2sin(t) - sin(2t))
            tri = generate_mesh(x, y, 0.1)
            tri2 = generate_mesh(x, y, 1.0)
            fig = Figure()
            ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", width=300, height=300,
                title=L"(a):$ $ Dense mesh", titlealign=:left)
            triplot!(ax, tri)
            ax = Axis(fig[1, 2], xlabel=L"x", ylabel=L"y", width=300, height=300,
                title=L"(b):$ $  Coarse mesh", titlealign=:left)
            triplot!(ax, tri2)
            resize_to_layout!(fig)
        end

        @testset "Another specific example with multiple segments" begin
            # The first segment 
            t = LinRange(0, 1 / 4, 25)
            x1 = cos.(2π * t)
            y1 = sin.(2π * t)
            # The second segment 
            t = LinRange(0, -3, 25)
            x2 = collect(t)
            y2 = repeat([1.0], length(t))
            # The third segment 
            t = LinRange(1, 0, 25)
            x3 = -3.0 .+ (1 .- t) .* sin.(t)
            y3 = collect(t)
            # The fourth segment 
            t = LinRange(0, 1, 25)
            x4 = collect(-3.0(1 .- t))
            y4 = collect(0.98t)
            # The fifth segment 
            x5 = [0.073914, 0.0797, 0.1522, 0.1522, 0.2, 0.28128, 0.3659, 0.4127, 0.3922, 0.4068, 0.497, 0.631, 0.728, 0.804, 0.888, 1.0]
            y5 = [0.8815, 0.8056, 0.80268, 0.73258, 0.6, 0.598, 0.5777, 0.525, 0.4346, 0.3645, 0.3032, 0.2886, 0.2623, 0.1367, 0.08127, 0.0]
            # Now combine the vectors 
            x = [x1, x2, x3, x4, x5]
            y = [y1, y2, y3, y4, y5]
            # Mesh 
            tri = generate_mesh(x, y, 0.05)
            fig = Figure()
            ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", width=600, height=300)
            triplot!(ax, tri)
            colors = [:red, :blue, :orange, :purple, :darkgreen]
            bn_map = get_boundary_map(tri)
            for (i, segment_index) in enumerate(values(bn_map))
                bn_nodes = get_boundary_nodes(tri, segment_index)
                lines!(ax, get_points(tri)[:, bn_nodes], color=colors[i], linewidth=4)
            end
            resize_to_layout!(fig)
        end

        @testset "Complicated example for docs" begin
            x1 = [collect(LinRange(0, 2, 4)),
                collect(LinRange(2, 2, 4)),
                collect(LinRange(2, 0, 4)),
                collect(LinRange(0, 0, 4))]
            y1 = [collect(LinRange(0, 0, 4)),
                collect(LinRange(0, 6, 4)),
                collect(LinRange(6, 6, 4)),
                collect(LinRange(6, 0, 4))]
            r = 0.5
            h = k = 0.6
            θ = LinRange(2π, 0, 50)
            x2 = [h .+ r .* cos.(θ)]
            y2 = [k .+ r .* sin.(θ)]
            r = 0.2
            h = 1.5
            k = 0.5
            x3 = [h .+ r .* cos.(θ)]
            y3 = [k .+ r .* sin.(θ)]
            x4 = reverse(reverse.([collect(LinRange(1, 1.5, 4)),
                collect(LinRange(1.5, 1.5, 4)),
                collect(LinRange(1.5, 1, 4)),
                collect(LinRange(1, 1, 4))]))
            y4 = reverse(reverse.([collect(LinRange(2, 2, 4)),
                collect(LinRange(2, 5, 4)),
                collect(LinRange(5, 5, 4)),
                collect(LinRange(5, 2, 4))]))
            x5 = [reverse([0.2, 0.5, 0.75, 0.75, 0.2, 0.2])]
            y5 = [reverse([2.0, 2.0, 3.0, 4.0, 5.0, 2.0])]
            x = [x1, x2, x3, x4, x5]
            y = [y1, y2, y3, y4, y5]
            tri = generate_mesh(x, y, 0.2)
            fig, ax, sc = triplot(tri; recompute_centers=true, show_ghost_edges=true, convex_hull_linestyle=:solid, convex_hull_linewidth=4)
            xlims!(ax, -0.5, 2.5)
            ylims!(ax, -0.5, 6.5)
        end
    end
end