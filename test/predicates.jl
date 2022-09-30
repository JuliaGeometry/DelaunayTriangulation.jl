@testset "Can we correctly identify the orientation of a triangle?" begin
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
    @test orient((4, 6, 7), pts) == 1
    @test orient((4, 7, 6), pts) == -1
    @test orient((4, 2, 3), pts) == -1
    @test orient((4, 7, 3), pts) == 1
    @test orient((5, 7, 9), pts) == 1
    @test orient((5, 9, 7), pts) == -1
    @test orient((3, 8, 5), pts) == -1
    @test orient((1, 2, 3), [[1.0, 2.0], [1.0, 5.0], [1.0, 8.0]]) == 0
end

@testset "Can we correct identify if a point is to the left of a line?" begin
    A = [4.6, 3.2]
    B = [3.2, 2.2]
    C = [3.4, 3.2]
    @test DT.leftofline(C, B, A) == 1
    @test DT.leftofline(C, A, B) == -1
    @test DT.leftofline([A, B], C, 2, 1) == 1
    @test DT.leftofline([A, B], C, 1, 2) == -1
    C = [5.8, 3.6]
    @test DT.leftofline(C, B, A) == -1
    A = [1.0, 7.0]
    B = [1.0, 1.0]
    C = [1.0, 5.0]
    @test DT.leftofline(C, B, A) == 0
    @test DT.leftofline(C, A, B) == 0
    @test DT.leftofline([A, B], C, 2, 1) == 0
    @test DT.leftofline([A, B], C, 1, 2) == 0
    A = [2.123933267613, 7.1892809338214]
    B = [-1.5542939635314, 3.3935384556756]
    C = [2.8172732249214, 5.085758012496]
    @test DT.leftofline(C, B, A) == -1
    @test DT.leftofline(C, A, B) == 1
    C = [-2.8172732249214, 5.085758012496]
    @test DT.leftofline(C, B, A) == 1
    @test DT.leftofline(C, A, B) == -1
    pts = [A, B]
    @test DT.leftofline(pts, C, 2, 1) == 1
    @test DT.leftofline(pts, C, 1, 2) == -1
end

@testset "Can we correctly test that a point is in a triangle, outside, or on the edge?" begin
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
    T1 = DT.TriangleType((4, 1, 6))
    T2 = DT.TriangleType((4, 2, 1))
    T3 = DT.TriangleType((3, 2, 4))
    T4 = DT.TriangleType((8, 1, 2))
    T5 = DT.TriangleType((8, 2, 3))
    T6 = DT.TriangleType((8, 3, 5))
    T7 = DT.TriangleType((5, 3, 7))
    T8 = DT.TriangleType((3, 4, 7))
    T9 = DT.TriangleType((5, 7, 9))
    T10 = DT.TriangleType((7, 6, 9))
    T11 = DT.TriangleType((7, 4, 6))
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
        @test orient(T[i], pts) == 1
        for p in test_pts[i]
            @test DT.intriangle(T[i], pts, collect(p)) == 1
        end
    end
end

@testset "Can we correctly compute whether a point is higher than another?" begin
    p = [2.0, 3.71]
    q = [2.0, 4.81]
    @test !DT.is_point_higher(p, q)
    q = [2.0, 3.60]
    @test DT.is_point_higher(p, q)
    q = [2.381, 3.71]
    @test DT.is_point_higher(p, q)
    p = [1.2999, 1.0]
    q = [1.981, 1.71]
    @test !DT.is_point_higher(p, q)
    @test DT.is_point_higher(q, p)
    p = [57.131, 4.0]
    q = [2.0, 3.1]
    @test DT.is_point_higher(p, q)
    @test DT.is_point_lower(q, p)
    p = [-2.31, 4.0]
    q = [5.0, 4.0]
    @test DT.is_point_higher(p, q)
    @test DT.is_point_lower(q, p)
end

@testset "Can we correctly sort points by height?" begin
    v = Vector{Vector{Float64}}(undef, 100)
    for _ in 1:500
        v .= [rand(2) for _ in 1:100]
        DT.partial_highest_point_sort!(v, 1)
        @test all(DT.is_point_higher(v[1], v[j]) for j in 2:lastindex(v))
    end
end

@testset "Can we correctly identify a given edge is on the bounding triangle?" begin
    @test DT.edge_on_large_triangle(1, 0)
    @test DT.edge_on_large_triangle(0, 1)
    @test DT.edge_on_large_triangle(0, -1)
    @test DT.edge_on_large_triangle(-1, 0)
    @test DT.edge_on_large_triangle(-1, 1)
    @test DT.edge_on_large_triangle(1, -1)
    @test !DT.edge_on_large_triangle(1, 5)
    @test !DT.edge_on_large_triangle(0, 2)
    @test !DT.edge_on_large_triangle(-1, 2)
    @test !DT.edge_on_large_triangle(0, -2)
end

@testset "Can we correctly test if a point is in a circle?" begin
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
    @test incircle(pts, 5, 7, 6, 9) == 1
    @test incircle(pts, 5, 7, 6, 3) == -1
    @test incircle(pts, 5, 7, 6, 3) == -1
    @test incircle(pts, 5, 7, 6, 6) == 0
    @test incircle(pts, 3, 2, 1, 4) == 1
    @test incircle(pts, 3, 2, 1, 6) == 1
    @test incircle(pts, 3, 2, 1, 7) == 1
    @test incircle(pts, 3, 2, 1, 5) == -1
    @test incircle(pts, 3, 2, 1, 8) == -1
end

@testset "Can we correctly test if an edge is legal?" begin
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
    𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, root = DT.initialise_triangulation()
    r = 2
    pᵣ = pts[r]
    Tᵢⱼₖ, flag = DT.locate_triangle(𝒟, pts, pᵣ, root)
    i, j, k = Tᵢⱼₖ
    DT.add_point!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, Tᵢⱼₖ, r)
    @test DT.is_legal(0, -1, 𝒜, pts)
    @test DT.is_legal(1, -1, 𝒜, pts)
    @test DT.is_legal(1, 0, 𝒜, pts)
    @test DT.is_legal(-1, 0, 𝒜, pts)
    @test DT.is_legal(-1, 1, 𝒜, pts)
    @test DT.is_legal(0, 1, 𝒜, pts)
    @test !DT.is_legal(7, 1, 6, 4, pts)
    @test DT.is_legal(7, 6, 9, 4, pts)
    @test DT.is_legal(7, 3, 4, 5, pts)
    @test DT.is_legal(7, 4, 6, 2, pts)
    @test DT.is_legal(k, i, 𝒜, pts)
end