############################################
##
## PREDICATES
##
############################################
@testset "Orient with coordinates rather than tuples" begin
    for _ in 1:2381
        pᵢ = 62rand(2)
        pⱼ = 37rand(2)
        pₖ = 87rand(2)
        @test orient(pᵢ, pⱼ, pₖ) == orient(pᵢ..., pⱼ..., pₖ...)
        pᵢ = 13 .* Tuple(rand(2))
        pⱼ = 5 .* Tuple(rand(2))
        pₖ = 53 .* Tuple(rand(2))
        @test orient(pᵢ, pⱼ, pₖ) == orient(pᵢ..., pⱼ..., pₖ...)
    end
end

@testset "Orientation of a triangle's points" begin
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
    @test DT.isoriented(DT.construct_triangle(NTuple{3,Int64}, 4, 6, 7), pts) == 1
    @test DT.isoriented(DT.construct_triangle(NTuple{3,Int64}, 4, 7, 6), pts) == -1
    @test DT.isoriented(DT.construct_triangle(NTuple{3,Int64}, 4, 2, 3), pts) == -1
    @test DT.isoriented(DT.construct_triangle(NTuple{3,Int64}, 4, 7, 3), pts) == 1
    @test DT.isoriented(DT.construct_triangle(NTuple{3,Int64}, 5, 7, 9), pts) == 1
    @test DT.isoriented(DT.construct_triangle(NTuple{3,Int64}, 5, 9, 7), pts) == -1
    @test DT.isoriented(DT.construct_triangle(NTuple{3,Int64}, 3, 8, 5), pts) == -1
    @test DT.isoriented(DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3), [[1.0, 2.0], [1.0, 5.0], [1.0, 8.0]]) == 0
    pts = [[0, -3.0], [3.0, 0.0], [0.0, 3.0], [-3.0, 0.0]]
    @test DT.isoriented(DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3), pts) == DT.isoriented(DT.construct_triangle(NTuple{3,Int64}, 2, 3, 4), pts) == DT.isoriented(DT.construct_triangle(NTuple{3,Int64}, 4, 1, 2), pts) == 1
    p₁ = [5.7044025422189, 1.801603986463]
    p₂ = [8.3797127128527, 5.8924221838871]
    p₃ = [2.8875415689061, 6.2038339497809]
    @test orient(p₁, p₂, p₃) == 1
    @test orient(p₁, p₂, [10.0, 4.0]) == -1
    p₁ = [5.0, 1.0]
    p₂ = [5.0, 6.0]
    @test orient(p₁, p₂, [5.0, 5.0]) == 0
    @test orient(p₁, p₂, [5.0, 2.0]) == 0
    @test orient(p₂, p₁, [5.0, 2.0]) == 0
end

@testset "Orientation of ghost edges" begin
    p1 = @SVector[-3.32, 3.53]
    p2 = @SVector[-5.98, 2.17]
    p3 = @SVector[-6.36, -1.55]
    pts = [p1, p2, p3]
    DT.compute_centroid!(pts)
    T = Set{NTuple{3,Int64}}()
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    DG = DT.DelaunayGraph{Int64}()
    DT.add_triangle!(1, 2, 3, T, adj, adj2v, DG; update_ghost_edges=true)
    @test all(DT.isoriented(T, pts) == 1 for T in T)
    DT.add_triangle!(T, (2, 3, 0))
    @test !all(DT.isoriented(T, pts) == 1 for T in T)
end

@testset "Testing if points are in a given circle" begin
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
    @test DT.isincircle(pts, 5, 7, 6, 9) == 1
    @test DT.isincircle(pts, 5, 7, 6, 3) == -1
    @test DT.isincircle(pts, 5, 7, 6, 3) == -1
    @test DT.isincircle(pts, 5, 7, 6, 6) == 0
    @test DT.isincircle(pts, 3, 2, 1, 4) == 1
    @test DT.isincircle(pts, 3, 2, 1, 6) == 1
    @test DT.isincircle(pts, 3, 2, 1, 7) == 1
    @test DT.isincircle(pts, 3, 2, 1, 5) == -1
    @test DT.isincircle(pts, 3, 2, 1, 8) == -1
    @test incircle(p4, p6, p5, p8) == 1
end

@testset "Testing if a point is to the left of a line" begin
    A = [4.6, 3.2]
    B = [3.2, 2.2]
    C = [3.4, 3.2]
    @test DT.isleftofline([A, B, C], 3, 2, 1) == 1
    @test DT.isleftofline([A, B, C], 3, 1, 2) == -1
    C = [5.8, 3.6]
    @test DT.isleftofline([A, B, C], 3, 2, 1) == -1
    A = [1.0, 7.0]
    B = [1.0, 1.0]
    C = [1.0, 5.0]
    @test DT.isleftofline([A, B, C], 3, 2, 1) == 0
    @test DT.isleftofline([A, B, C], 3, 1, 2) == 0
    A = [2.123933267613, 7.1892809338214]
    B = [-1.5542939635314, 3.3935384556756]
    C = [2.8172732249214, 5.085758012496]
    @test DT.isleftofline([A, B, C], 3, 2, 1) == -1
    @test DT.isleftofline([A, B, C], 3, 1, 2) == 1
    C = [-2.8172732249214, 5.085758012496]
    @test DT.isleftofline([A, B, C], 3, 2, 1) == 1
    @test DT.isleftofline([A, B, C], 3, 1, 2) == -1
    @test DT.isleftofline([A, B, C], 3, 2, 1) == 1
    @test DT.isleftofline([A, B, C], 3, 1, 2) == -1
    D = [-6.3, -2.77]
    E = [-3.46, -2.07]
    F = [-5.22, -1.41]
    @test DT.isleftofline([D, E, F], 1, 2, 3) == 1
    @test DT.isleftofline([D, E, F], 2, 1, 3) == -1
    F = [-4.88, -4.57]
    @test DT.isleftofline([D, E, F], 1, 2, 3) == -1
    @test DT.isleftofline([D, E, F], 2, 1, 3) == 1
end

@testset "Testing if a point is to the left of a line defined by edges of the bounding triangle" begin
    p1 = [2.0, 3.0]
    p2 = [5.7, 2.3]
    p3 = [17.0, -2.0]
    p = [0.0, 0.0]
    pts = [p1, p2, p3, p]
    i, j = DT.LowerRightBoundingIndex, DT.UpperBoundingIndex
    @test DT.isleftofline(pts, i, j, 4) == 1
    i, j = DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex
    @test DT.isleftofline(pts, i, j, 4) == -1
    i, j = DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex
    @test DT.isleftofline(pts, i, j, 4) == 1
    i, j = DT.UpperBoundingIndex, DT.LowerRightBoundingIndex
    @test DT.isleftofline(pts, i, j, 4) == -1
    i, j = DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex
    @test DT.isleftofline(pts, i, j, 4) == 1
    i, j = DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex
    @test DT.isleftofline(pts, i, j, 4) == -1
end

@testset "Testing if we can find where a point lies to a segment connecting to the bounding triangle" begin
    pᵢ = [17.0, 2.0]
    p = [2.0, -2.5]
    pts = [pᵢ, p]
    @test DT.isleftofline(pts, 1, DT.LowerRightBoundingIndex, 2) == -1
    @test DT.isleftofline(pts, DT.LowerRightBoundingIndex, 1, 2) == 1
    pᵢ = [-8.7999865, 21.89396]
    pts = [pᵢ, p]
    @test DT.isleftofline(pts, 1, DT.LowerRightBoundingIndex, 2) == -1
    @test DT.isleftofline(pts, DT.LowerRightBoundingIndex, 1, 2) == 1
    pᵢ = [-12.49835, -10.738]
    pts = [pᵢ, p]
    @test DT.isleftofline(pts, 1, DT.LowerRightBoundingIndex, 2) == 1
    @test DT.isleftofline(pts, DT.LowerRightBoundingIndex, 1, 2) == -1
    pᵢ = [20.0, 20.0]
    p = [0.0, 20.0]
    pts = [pᵢ, p]
    @test DT.isleftofline(pts, 1, DT.LowerRightBoundingIndex, 2) == -1
    @test DT.isleftofline(pts, DT.LowerRightBoundingIndex, 1, 2) == 1
    pᵢ = [20.569333, 8.405]
    p = [8.743, -17.13]
    pts = [pᵢ, p]
    @test DT.isleftofline(pts, 1, DT.LowerLeftBoundingIndex, 2) == 1
    @test DT.isleftofline(pts, DT.LowerLeftBoundingIndex, 1, 2) == -1
    p = [-98.23, 26.9]
    pts = [pᵢ, p]
    @test DT.isleftofline(pts, 1, DT.LowerLeftBoundingIndex, 2) == -1
    @test DT.isleftofline(pts, DT.LowerLeftBoundingIndex, 1, 2) == 1
    @test DT.isleftofline(pts, 1, DT.UpperBoundingIndex, 2) == 1
    @test DT.isleftofline(pts, DT.UpperBoundingIndex, 1, 2) == -1
end

@testset "Testing Boolean intriangle values" begin
    e = [
        1 1 1
        1 1 0
        1 1 -1
        1 0 1
        1 0 0
        1 0 -1
        1 -1 1
        1 -1 0
        1 -1 -1
        0 1 1
        0 1 0
        0 1 -1
        0 0 1
        0 0 0
        0 0 -1
        0 -1 1
        0 -1 0
        0 -1 -1
        -1 1 1
        -1 1 0
        -1 1 -1
        -1 0 1
        -1 0 0
        -1 0 -1
        -1 -1 1
        -1 -1 0
        -1 -1 -1
    ]
    ev = [1, 0, -1, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, -1, 0, -1]
    @test all(DT.isintriangle(e[i, 1], e[i, 2], e[i, 3]) == ev[i] for i in eachindex(ev))
end

@testset "Testing that a point is in a triangle, outside, or on the edge" begin
    p1 = Float64[5, 5]
    p2 = Float64[4.5, 2.5]
    p3 = Float64[2.5, 1.5]
    p4 = Float64[3, 3.5]
    p5 = Float64[0, 2]
    p6 = Float64[1, 5]
    p7 = Float64[1, 3]
    p8 = Float64[4, -1]
    p9 = Float64[-1, 4]
    pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9]
    T1 = DT.construct_triangle(NTuple{3,Int64}, 4, 1, 6)
    T2 = DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1)
    T3 = DT.construct_triangle(NTuple{3,Int64}, 3, 2, 4)
    T4 = DT.construct_triangle(NTuple{3,Int64}, 8, 1, 2)
    T5 = DT.construct_triangle(NTuple{3,Int64}, 8, 2, 3)
    T6 = DT.construct_triangle(NTuple{3,Int64}, 8, 3, 5)
    T7 = DT.construct_triangle(NTuple{3,Int64}, 5, 3, 7)
    T8 = DT.construct_triangle(NTuple{3,Int64}, 3, 4, 7)
    T9 = DT.construct_triangle(NTuple{3,Int64}, 5, 7, 9)
    T10 = DT.construct_triangle(NTuple{3,Int64}, 7, 6, 9)
    T11 = DT.construct_triangle(NTuple{3,Int64}, 7, 4, 6)
    T = [T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11]
    pt1 = (
        A=(3.9551298987095, 4.7489988803935),
        B=(3.2585303811361, 4.3003415639903),
        C=(2.69180534989, 4.4066025073489),
        D=(3.2231100666833, 4.6781582514877),
        E=(2.0424329182538, 4.8198395092992)
    )
    pt2 = (
        A=(3.5, 3.5),
        B=(3.8384340736083, 3.2777159548861),
        C=(4.173119090818, 3.0606229707501),
        D=(4.0736181397556, 3.8475850382431),
        E=(4.390212074954, 3.6938108411468),
        F=(4.390212074954, 4.2546343834981),
        G=(4.7520337151806, 4.2998620885264),
        H=(4.7520337151806, 4.661683728753)
    )
    pt3 = (
        A=(3.114790793155, 2.9520764786821),
        B=(3.3952025643307, 2.8163933635972),
        C=(3.3137926952797, 2.4545717233705),
        D=(2.951971055053, 2.3098430672799),
        E=(2.9157888910304, 2.4726628053818),
        F=(3.0786086291324, 2.6535736254952),
        G=(3.901752860648, 2.6626191665008),
        H=(3.0, 2.0),
        I=(2.7258325299114, 1.8213838529739)
    )
    pt4 = (
        A=(4.5781396673455, 2.7004728027825),
        B=(4.6264236360138, 2.9231155471972),
        C=(4.6693427192745, 3.1618529478347),
        D=(4.70153203172, 3.3871781349531),
        E=(4.5325381413811, 2.437593417811),
        F=(4.4708419591939, 2.1988560171735),
        G=(4.4520648602673, 2.0352270122422),
        H=(4.3501320375233, 1.3082850395147)
    )
    pt5 = (
        A=(3.433718968984, 1.1294534677817),
        B=(3.7382811110409, 1.4137114670348),
        C=(3.8499538964617, 0.9162599683419),
        D=(3.6875207540314, 0.5304812550699),
        E=(3.0377881843101, 1.2614303960064),
        F=(3.7484331824428, -0.0481868148381),
        G=(4.0, 2.0)
    )
    pt6 = (
        A=(2.8956591846836, 0.4086563982472),
        B=(2.4286639001964, 0.90610789694),
        C=(2.0936455439339, 1.2309741818007),
        D=(1.5149774740259, 1.4035593956329),
        E=(0.824636618697, 1.586296680867),
        F=(1.3322401887918, 1.2715824674082),
        G=(1.6063461166429, 1.0177806823609),
        H=(2.0225810441206, 0.713218540304),
        I=(2.3169911147756, 0.4391126124529),
        J=(2.6317053282343, 0.1954628988074),
        K=(3.1697651125348, -0.1395554574552)
    )
    pt7 = (
        A=(1.0581342609406, 2.2766375361959),
        B=(0.9363094041178, 2.550743464047),
        C=(1.2916319031842, 2.3070937504015),
        D=(1.7281709734657, 1.9923795369428),
        E=(1.0, 2.0),
        F=(0.5809869050515, 2.1142043937655)
    )
    pt8 = (
        A=(1.5454336882315, 2.845153534702),
        B=(1.9921248299149, 2.5405913926451),
        C=(2.1545579723453, 2.936522177319),
        D=(2.1951662579528, 2.4187665358224),
        E=(2.4185118287945, 2.8756097489077),
        F=(2.4185118287945, 2.5101351784394),
        G=(2.4489680430002, 2.0431398939523),
        H=(2.6317053282343, 2.9162180345153)
    )
    pt9 = (
        A=(-0.5458930205588, 3.5557985328346),
        B=(0.0733833349568, 3.2512363907778),
        C=(0.3170330486022, 2.9771304629266),
        D=(0.0835354063587, 2.6421121066641),
        E=(0.0, 2.4187665358224),
        F=(0.3576413342098, 2.4695268928319),
        G=(-0.2210267356982, 2.9872825343285),
        H=(-0.4849805921475, 3.2410843193759),
        I=(0.5099224052383, 2.9162180345153)
    )
    pt10 = (
        A=(0.3576413342098, 4.1649228169483),
        B=(0.0, 4.0),
        C=(0.4794661910326, 3.6573192468536),
        D=(0.7028117618743, 3.3527571047967),
        E=(0.6520514048648, 4.4796370304071),
        F=(-0.3530036639228, 4.1243145313408),
        G=(0.0, 3.7689920322744)
    )
    pt11 = (
        A=(1.3931526172031, 4.0735541743313),
        B=(2.0022769013168, 3.718231675265),
        C=(1.3931526172031, 3.5354943900309),
        D=(2.1444059009434, 3.4339736760119),
        E=(1.3017839745861, 4.3172038879768),
        F=(1.3017839745861, 4.5507015302204),
        G=(1.7992354732789, 3.9923376031161),
        H=(1.6875626878581, 3.5151902472271),
        I=(1.4337609028107, 3.809600317882)
    )
    test_pts = [pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8, pt9, pt10, pt11]
    for i in eachindex(T)
        @test DT.isoriented(T[i], pts) == 1
        for r in eachindex(test_pts[i])
            push!(pts, [test_pts[i][r]...])
            @test DT.isintriangle(T[i], pts, length(pts)) == 1
            pop!(pts)
        end
    end
    push!(pts, [2.0, 5.7])
    @test DT.isintriangle(DT.BoundingTriangle, pts, length(pts)) == 1
    push!(pts, [-2.7, 29.5])
    @test DT.isintriangle(DT.BoundingTriangle, pts, length(pts)) == 1
    push!(pts, [2.422, 188.2])
    @test DT.isintriangle(DT.BoundingTriangle, pts, length(pts)) == 1
    push!(pts, [172.3, 178.0])
    @test DT.isintriangle(DT.BoundingTriangle, pts, length(pts)) == 1
    push!(pts, [49.1, 1720.0])
    @test DT.isintriangle(DT.BoundingTriangle, pts, length(pts)) == 1
    push!(pts, [-2.02, 13.4])
    @test DT.isintriangle(DT.BoundingTriangle, pts, length(pts)) == 1
    push!(pts, [2.0, 5.7])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex), pts, length(pts)) == 1
    push!(pts, [-2.7, 29.5])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex), pts, length(pts)) == 1
    push!(pts, [2.422, 188.2])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex), pts, length(pts)) == 1
    push!(pts, [172.3, 178.0])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex), pts, length(pts)) == 1
    push!(pts, [49.1, 1720.0])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex), pts, length(pts)) == 1
    push!(pts, [-2.02, 13.4])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex), pts, length(pts)) == 1
    push!(pts, [2.0, 5.7])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex), pts, length(pts)) == 1
    push!(pts, [-2.7, 29.5])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex), pts, length(pts)) == 1
    push!(pts, [2.422, 188.2])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex), pts, length(pts)) == 1
    push!(pts, [172.3, 178.0])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex), pts, length(pts)) == 1
    push!(pts, [49.1, 1720.0])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex), pts, length(pts)) == 1
    push!(pts, [-2.02, 13.4])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex), pts, length(pts)) == 1
end

@testset "Test if edge is on the bounding triangle" begin
    @test DT.edge_on_bounding_triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)
    @test DT.edge_on_bounding_triangle(DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex)
    @test DT.edge_on_bounding_triangle(DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex)
    @test DT.edge_on_bounding_triangle(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    @test DT.edge_on_bounding_triangle(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
    @test DT.edge_on_bounding_triangle(DT.UpperBoundingIndex, DT.LowerRightBoundingIndex)
    @test !DT.edge_on_bounding_triangle(2, 7)
    @test !DT.edge_on_bounding_triangle(3, 9)
    @test !DT.edge_on_bounding_triangle(DT.LowerLeftBoundingIndex, 2)
    @test !DT.edge_on_bounding_triangle(20, DT.LowerLeftBoundingIndex)
    @test !DT.edge_on_bounding_triangle(10, DT.LowerRightBoundingIndex)
    @test !DT.edge_on_bounding_triangle(DT.LowerRightBoundingIndex, 18)
    @test !DT.edge_on_bounding_triangle(DT.UpperBoundingIndex, 8)
    @test !DT.edge_on_bounding_triangle(91, DT.UpperBoundingIndex)
end

@testset "Test if a given edge is legal" begin
    p1 = [-1.0, 4.0]
    p2 = [4.0, 6.0]
    p3 = [2.0, 2.0]
    p4 = [6.0, -3.0]
    p5 = [7.0, 3.0]
    p6 = [-2.5519459826976, 13.6106262700637]
    p7 = [-21.0502221073507, -23.3458204355075]
    p8 = [23.7055438911111, 19.1906513812123]
    p9 = [31.7813088302274, -22.9759380718838]
    pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9]
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 1, 3, 2))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 1, 4, 3))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 3, 4, 5))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 3, 5, 2))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 1, 7, 4))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 7, 9, 4))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 4, 9, 5))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 5, 9, 8))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 2, 5, 8))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 2, 8, 6))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 6, 1, 2))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 6, 7, 1))
    @test DT.islegal(DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, adj, pts)
    @test DT.islegal(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, adj, pts)
    @test DT.islegal(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, adj, pts)
    @test DT.islegal(DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex, adj, pts)
    @test DT.islegal(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, adj, pts)
    @test DT.islegal(DT.UpperBoundingIndex, DT.LowerRightBoundingIndex, adj, pts)
    @test all(DT.islegal(i, j, adj, pts) for (i, j) in ((1, 3), (3, 2), (2, 1), (1, 4), (4, 3), (3, 1), (3, 4), (4, 5), (5, 3), (3, 5), (5, 2), (2, 3)))
    @test_broken !DT.islegal(2, DT.LowerRightBoundingIndex, 17, 3, pts)
    @test_broken !DT.islegal(DT.LowerRightBoundingIndex, 2, 17, 3, pts)
    @test_broken DT.islegal(17, 5, DT.LowerRightBoundingIndex, 5, pts)
    @test_broken DT.islegal(20, 21, 17, DT.LowerRightBoundingIndex, pts)
    @test_broken !DT.islegal(2, DT.LowerLeftBoundingIndex, 1117, 36, pts)
    @test_broken !DT.islegal(DT.LowerLeftBoundingIndex, 223, 5117, 613, pts)
    @test_broken DT.islegal(147, 235, DT.LowerLeftBoundingIndex, 15, pts)
    @test_broken DT.islegal(2610, 1221, 1517, DT.LowerLeftBoundingIndex, pts)
    @test_broken !DT.islegal(2, DT.UpperBoundingIndex, 517, 63, pts)
    @test_broken !DT.islegal(DT.UpperBoundingIndex, 2, 117, 36, pts)
    @test_broken DT.islegal(157, 51, DT.UpperBoundingIndex, 35, pts)
    @test_broken DT.islegal(205, 2121, 174, DT.UpperBoundingIndex, pts)
    @test_broken !DT.islegal(-2, 3, 1, -5, pts)
    @test_broken DT.islegal(-5, 53, 51, -3, pts)
    @test_broken !DT.islegal(4, -2, 51, -3, pts)
    @test_broken DT.islegal(4, -12, 51, -7, pts)
    @test_broken !DT.islegal(-2, 3, -5, 17, pts)
    @test_broken DT.islegal(-5, 53, -3, 18, pts)
    @test_broken !DT.islegal(4, -2, -3, 172, pts)
    @test_broken DT.islegal(4, -12, -7, 80, pts)
end

@testset "Can we find the correct edge that a point is on for a given triangle?" begin
    p1 = [2.0, 3.5]
    p2 = [0.0, 0.0]
    p3 = [3.0, 0.0]
    p4 = [17.2, -2.5]
    p5 = [0.0, 3.0]
    T = DT.construct_triangle(NTuple{3,Int64}, 2, 3, 5)
    pts = [p1, p2, p3, p4, p5]
    push!(pts, [1.0, 0.0])
    @test DT.find_edge(T, pts, length(pts)) == (2, 3)
    push!(pts, [2.0, 0.0])
    @test DT.find_edge(T, pts, length(pts)) == (2, 3)
    push!(pts, [1.5, 0.0])
    @test DT.find_edge(T, pts, length(pts)) == (2, 3)
    push!(pts, [1.0, 0.0])
    @test DT.find_edge(T, pts, length(pts)) == (2, 3)
    push!(pts, [0.5, 0.0])
    @test DT.find_edge(T, pts, length(pts)) == (2, 3)
    push!(pts, [2.5, 0.5])
    @test DT.find_edge(T, pts, length(pts)) == (3, 5)
    push!(pts, [2.0, 1.0])
    @test DT.find_edge(T, pts, length(pts)) == (3, 5)
    push!(pts, [1.5, 1.5])
    @test DT.find_edge(T, pts, length(pts)) == (3, 5)
    push!(pts, [1.0, 2.0])
    @test DT.find_edge(T, pts, length(pts)) == (3, 5)
    push!(pts, [0.5, 2.5])
    @test DT.find_edge(T, pts, length(pts)) == (3, 5)
    push!(pts, [0.0, 2.5])
    @test DT.find_edge(T, pts, length(pts)) == (5, 2)
    push!(pts, [0.0, 2.2])
    @test DT.find_edge(T, pts, length(pts)) == (5, 2)
    push!(pts, [0.0, 2.0])
    @test DT.find_edge(T, pts, length(pts)) == (5, 2)
    push!(pts, [0.0, 1.5])
    @test DT.find_edge(T, pts, length(pts)) == (5, 2)
    push!(pts, [0.0, 0.8])
    @test DT.find_edge(T, pts, length(pts)) == (5, 2)
    push!(pts, [0.0, 0.2])
    @test DT.find_edge(T, pts, length(pts)) == (5, 2)
end

@testset "Testing if a boundary edge is legal" begin
    pts, T, DG, adj, adj2v = example_triangulation()
    hg = DT.HistoryGraph{NTuple{3,Int64}}()
    DT.add_triangle!(hg, T...)
    p1, p2, p3, p4, p5, p6 = pts
    p7 = @SVector[2.0, 1.0]
    DT.add_point!(pts, p7)
    DT.add_triangle!(6, 2, 3, T, adj, adj2v, DG)
    DT.split_triangle!(1, 3, 5, 7, T, adj, adj2v, DG)
    DT.flip_edge!(5, 1, 4, 7, T, adj, adj2v, DG, hg)
    for ((i, j), v) in adjacent(adj)
        if v ≠ DT.DefaultAdjacentValue
            @test DT.islegal(i, j, adj, pts)
        end
    end
end