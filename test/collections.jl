############################################
##
## PRIMITIVE COLLECTIONS 
##
############################################
@testset "Triangles" begin
    T₁ = Triangle(1, 2, 3)
    T₂ = Triangle(5, 8, 2)
    T₃ = Triangle(17, 29, 15)
    T₄ = Triangle(-2, 0, 5)
    T₅ = Triangle(17, 30, 72)
    T₆ = Triangle(5, 50, 101)
    T_set = Set{Triangle{Int}}([T₁, T₂, T₃, T₄, T₅, T₆])
    T_set_struct = Triangles(T_set)
    @test triangles(T_set_struct) == T_set
    T_vec = [T₁, T₂, T₃, T₄, T₅, T₆]
    T_vec_struct = Triangles(T_vec)
    @test triangles(T_vec_struct) == T_set
    T_tup = [T₁, T₂, T₃, T₄, T₅, T₆]
    T_tup_struct = Triangles(T_tup)
    @test triangles(T_tup_struct) == T_set
    @test isconcretetype(typeof(T_tup_struct))
    T_set_param = Triangles{Int64,Triangle{Int64}}(T_set)
    @test T_set_param == T_set_struct
    @test Triangles(T₁, T₂, T₃, T₄, T₅, T₆).triangles == T_set

    T = Triangle(1, 2, 3)
    Ts = Triangles([T])
    DT.add_triangle!(Ts, Triangle(4, 5, 7))
    @test Ts.triangles == Set{Triangle{Int64}}([Triangle(1, 2, 3), Triangle(4, 5, 7)])
    DT.delete_triangle!(Ts, Triangle(1, 2, 3))
    @test Ts.triangles == Set{Triangle{Int64}}([Triangle(4, 5, 7)])
    DT.add_triangle!(Ts, Triangle(5, 3, 2), Triangle(10, 11, 12), Triangle(13, 15, 19))
    @test Ts.triangles == Set{Triangle{Int64}}([Triangle(4, 5, 7), Triangle(5, 3, 2), Triangle(10, 11, 12), Triangle(13, 15, 19)])
    DT.delete_triangle!(Ts, Triangle(2, 5, 3), Triangle(19, 13, 15))
    @test Ts.triangles == Set{Triangle{Int64}}([Triangle(4, 5, 7), Triangle(10, 11, 12)])
    DT.delete_triangle!(Ts, Triangle(4, 5, 7))
    @test Ts.triangles == Set{Triangle{Int64}}([Triangle(10, 11, 12)])

    T₁ = Triangle(1, 2, 3)
    T₂ = Triangle(5, 8, 2)
    T₃ = Triangle(17, 29, 15)
    T₄ = Triangle(-2, 0, 5)
    T₅ = Triangle(17, 30, 72)
    T₆ = Triangle(5, 50, 101)
    T_set = Set{Triangle{Int}}([T₁, T₂, T₃, T₄, T₅, T₆])
    T_tris = Triangles(T_set)
    T_testset = Set{Triangle{Int}}()
    for T in T_tris
        push!(T_testset, T)
    end
    @test T_set == T_testset
    @test length(T_tris) == 6
    @test eltype(T_set) == Triangle{Int64}
    @test eltype(Triangles{Int,Triangle{Int}}) == Triangle{Int64}
    T_testset2 = Set{Triangle{Int}}()
    T_testsetint2 = Int64[]
    for (i, T) in enumerate(T_tris)
        push!(T_testset2, T)
        push!(T_testsetint2, i)
    end
    @test T_set == T_testset2
    @test T_testsetint2 == [1, 2, 3, 4, 5, 6]
end

@testset "Points" begin
    p₁ = Point(2.0, 5.0)
    p₂ = Point(5.0, 1.7)
    p₃ = Point(2.2, 2.2)
    p₄ = Point(-17.0, 5.0)
    pts_vec = [p₁, p₂, p₃, p₄]
    pts_vec_struct = Points(pts_vec)
    @test size(pts_vec_struct) == (4,)
    @test eachindex(pts_vec_struct) == Base.OneTo(4)
    @test points(pts_vec_struct) == pts_vec
    pts_tup = (p₁, p₂, p₃, p₄)
    pts_tup_struct = Points(pts_tup)
    @test points(pts_tup_struct) == pts_vec
    @test isconcretetype(typeof(pts_tup_struct))
    @test DT.get_point(pts_tup_struct, 1) == p₁
    @test DT.get_point(pts_tup_struct, 2) == p₂
    @test DT.get_point(pts_tup_struct, 3) == p₃
    @test DT.get_point(pts_tup_struct, 4) == p₄
    @test length(pts_tup_struct) == 4
    add_point!(pts_tup_struct, Point(5.7, 13.3))
    @test points(pts_tup_struct) == [p₁, p₂, p₃, p₄, Point(5.7, 13.3)]
    @test length(pts_tup_struct) == 5

    p0 = Float64[5, 5]
    p1 = Float64[4.5, 2.5]
    p2 = Float64[2.5, 1.5]
    p3 = Float64[3, 3.5]
    p4 = Float64[0, 2]
    p5 = Float64[1, 5]
    p6 = Float64[1, 3]
    p7 = Float64[4, -1]
    p8 = Float64[-1, 4]
    pts = [p0, p1, p2, p3, p4, p5, p6, p7, p8]
    qts = Points(pts)
    @test getx.(points(qts)) == first.(pts)
    @test gety.(points(qts)) == last.(pts)
    rts = Points(p0, p1, p2, p3, p4, p5, p6, p7, p8)
    @test getx.(points(rts)) == first.(pts)
    @test gety.(points(rts)) == last.(pts)
    @test typeof(rts) == Points{Float64,Vector{Float64},Point{Float64,Vector{Float64}}}

    pt_test = Vector{Point{Float64,Vector{Float64}}}()
    for pt in rts
        push!(pt_test, pt)
    end
    @test pt_test == rts.points
    @test length(rts) == 9
    @test eltype(rts) == Point{Float64,Vector{Float64}}
    @test eltype(Points{Float64,Vector{Float64},Point{Float64,Vector{Float64}}}) == Point{Float64,Vector{Float64}}

    p₁ = Point(2.0, 5.0)
    p₂ = Point(5.0, 1.7)
    p₃ = Point(2.2, 2.2)
    p₄ = Point(-17.0, 5.0)
    pts_vec = [p₁, p₂, p₃, p₄]
    pts = Points(pts_vec)
    @test pts.points == pts_vec
    @test DT.xcentroid(pts) == (-17.0 + 5.0) / 2
    @test DT.ycentroid(pts) == (5.0 + 1.7) / 2
    @test DT.xmin(pts) == -17.0
    @test DT.xmax(pts) == 5.0
    @test DT.ymin(pts) == 1.7
    @test DT.ymax(pts) == 5.0
    @test DT.width(pts) == 22.0
    @test DT.height(pts) == 3.3
    @test DT.max_width_height(pts) == 22.0
    @test DT.get_point(pts, 1) == p₁
    @test DT.get_point(pts, 2) == p₂
    @test DT.get_point(pts, 3) == p₃
    @test DT.get_point(pts, 4) == p₄
    @test DT.get_point(pts, DT.LowerRightBoundingIndex) == Point(-6.0 + DT.BoundingTriangleShift * 22.0, 3.35 - 22.0)
    @test DT.get_point(pts, DT.LowerLeftBoundingIndex) == Point(-6.0 - DT.BoundingTriangleShift * 22.0, 3.35 - 22.0)
    @test DT.get_point(pts, DT.UpperBoundingIndex) == Point(-6.0, 3.35 + DT.BoundingTriangleShift * 22.0)
    @test DT.get_point(pts, DT.LowerRightBoundingIndex) == pts.lower_right_bounding_triangle_coords
    @test DT.get_point(pts, DT.LowerLeftBoundingIndex) == pts.lower_left_bounding_triangle_coords
    @test DT.get_point(pts, DT.UpperBoundingIndex) == pts.upper_bounding_triangle_coords
    @test_throws BoundsError DT.get_point(pts, 0)
    @test_throws BoundsError DT.get_point(pts, -5)
    @test_throws BoundsError DT.get_point(pts, 17)

    pts_copy = deepcopy(pts)
    Random.seed!(292991)
    shuffle!(pts_copy)
    pts_copy2 = deepcopy(pts_vec)
    Random.seed!(292991)
    shuffle!(pts_copy2)
    @test pts_copy.points == pts_copy2

    p₁ = Point(2.0, 5.0)
    p₂ = Point(5.0, 1.7)
    p₃ = Point(2.2, 2.2)
    p₄ = Point(-17.0, 5.0)
    pts_vec = [p₁, p₂, p₃, p₄]
    pts = Points(pts_vec)
    @test firstindex(pts)==1
    @test lastindex(pts)==4
end
