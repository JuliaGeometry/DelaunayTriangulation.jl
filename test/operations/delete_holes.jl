using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using StableRNGs
using StaticArrays
using Preferences

@testset "Simple domain" begin
    p1 = (0.2, 0.2)
    p2 = (0.8, 0.2)
    p3 = (0.8, 0.8)
    p4 = (0.2, 0.8)
    rng = StableRNG(1919)
    existing_pts = [(rand(rng), rand(rng)) for _ in 1:15]
    push!(existing_pts, (0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (1.0, 1.0))
    nodes, points = convert_boundary_points_to_indices([p1, p2, p3, p4, p1], existing_points=existing_pts)
    tri = triangulate(points; boundary_nodes=nodes, delete_holes=false, delete_ghosts=false, rng)
    points_to_process = DT.find_all_points_to_delete(tri)
    @test points_to_process == Set((17, 1, DT.𝒢, 19, 12, 3, 8, 16, 15, 5, 11, 18, 9, 4))
    triangles_to_delete = DT.find_all_triangles_to_delete(tri, points_to_process)
    true_triangles_to_delete = Set((
        (1, 22, 4),
        (1, 4, 19),
        (1, 19, 17),
        (1, 17, 12),
        (1, 12, 22),
        (3, 8, 23),
        (3, 23, 17),
        (3, 17, 8),
        (4, 11, 19),
        (4, 22, 21),
        (4, 21, 11),
        (5, 18, 11),
        (5, 11, 21),
        (5, 21, 9),
        (5, 9, 15),
        (5, 15, 16),
        (8, 16, 20),
        (8, 17, 16),
        (8, 20, 23),
        (9, 20, 15),
        (9, 21, 20),
        (5, 16, 18),
        (11, 18, 19),
        (12, 23, 22),
        (12, 17, 23),
        (15, 20, 16),
        (17, 19, DT.𝒢),
        (19, 18, DT.𝒢),
        (18, 16, DT.𝒢),
        (16, 17, DT.𝒢)
    ))
    @test DT.compare_triangle_collections(triangles_to_delete, true_triangles_to_delete)

    p1 = (0.2, 0.2)
    p2 = (0.8, 0.2)
    p3 = (0.8, 0.8)
    p4 = (0.2, 0.8)
    rng = StableRNG(1919)
    existing_pts = [(rand(rng), rand(rng)) for _ in 1:15]
    push!(existing_pts, (0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (1.0, 1.0))
    nodes, points = convert_boundary_points_to_indices([p1, p2, p3, p4, p1], existing_points=existing_pts)
    tri = triangulate(points; boundary_nodes=nodes, rng)
    @test validate_triangulation(tri)
    validate_statistics(tri)
end

@testset "Multiply-connected" begin
    a = (0.0, 0.0)
    b = (5.0, 0.0)
    c = (5.0, 5.0)
    d = (0.0, 5.0)
    e = (1.0, 4.0)
    f = (0.0, 3.0)
    g = (1.5, 2.0)
    h = (4.0, 2.0)
    i = (4.0, 1.0)
    j = (1.5, 1.0)
    k = (2.0, 1.5)
    ℓ = (2.4, 1.0)
    m = (2.8, 3.6)
    n = (2.2, 2.8)
    o = (7.2, 2.4)
    p = (7.3, 3.9)
    q = (2.0, 1.4)
    r = (3.0, 1.6)
    s = (3.4, 1.4)
    t = (3.6, 2.2)
    u = (1.0, 2.8)
    v = (0.8, 1.8)
    w = (0.8, 4.0)
    z = (0.4, 4.4)
    c1 = (5.0, 4.6)
    d1 = (4.9, 4.8)
    a1 = (5.6, 1.8)
    b1 = (6.0, 0.6)
    e1 = (5.1, 4.4)
    f1 = (5.0, 4.3)
    pts = [[[d, e, f, a, b, f1, e1, c1, d1, c, d]], [[g, h, i, ℓ, k, j, g]]]
    existing_points = [m, n, o, p, q, r, s, t, u, v, w, z, a1, b1]
    nodes, points = convert_boundary_points_to_indices(pts; existing_points=existing_points)
    tri = triangulate(points; boundary_nodes=nodes, delete_holes=false, delete_ghosts=false)
    points_to_process = DT.find_all_points_to_delete(tri)
    @test points_to_process == Set((12, 11, 14, 3, 13, 4, 6, 7, DT.𝒢))
    triangles_to_delete = DT.find_all_triangles_to_delete(tri, points_to_process)
    true_triangles_to_delete = Set((
        (15, 17, 12),
        (12, 17, 11),
        (12, 11, 16),
        (11, 17, 16),
        (25, 30, 29),
        (29, 28, 6),
        (25, 29, 6),
        (25, 6, 26),
        (6, 7, 26),
        (6, 28, 7),
        (7, 28, 27),
        (26, 7, 27),
        (20, 19, 13),
        (13, 19, 14),
        (13, 14, 3),
        (20, 13, 3),
        (20, 3, 4),
        (21, 20, 4),
        (24, 23, 22),
        (24, 22, 21),
        (24, 21, 4),
        (15, 24, DT.𝒢),
        (24, 4, DT.𝒢),
        (4, 3, DT.𝒢),
        (3, 14, DT.𝒢),
        (14, 19, DT.𝒢),
        (19, 18, DT.𝒢),
        (18, 17, DT.𝒢),
        (17, 15, DT.𝒢),
        (12, 16, 15)
    ))
    @test DT.compare_triangle_collections(triangles_to_delete, true_triangles_to_delete)

    tri = triangulate(points; boundary_nodes=nodes)
    @test validate_triangulation(tri)
    validate_statistics(tri)
end

if load_preference(DelaunayTriangulation, "USE_EXACTPREDICATES", true)
    @testset "Interior holes that were already triangles" begin
        p1 = (0.0, 0.0)
        p2 = (1.0, 0.0)
        p3 = (1.0, 1.0)
        p4 = (0.0, 1.0)
        p5 = (0.2, 0.2)
        p6 = (0.8, 0.2)
        p7 = (0.8, 0.8)
        pts = [p1, p2, p3, p4, p5, p6, p7]
        tri = triangulate(pts, boundary_nodes=[[[1, 2, 3, 4, 1]], [[5, 7, 6, 5]]], delete_holes=false, delete_ghosts=false)
        @test num_ghost_triangles(tri) == 4
        @test num_ghost_edges(tri) == 4
        points_to_process = DT.find_all_points_to_delete(tri)
        @test points_to_process == Set((DT.𝒢,))
        triangles_to_delete = DT.find_all_triangles_to_delete(tri, points_to_process)
        @test DT.compare_triangle_collections(triangles_to_delete, Set((
            (1, 4, DT.𝒢),
            (4, 3, DT.𝒢),
            (3, 2, DT.𝒢),
            (2, 1, DT.𝒢),
            (5, 6, 7)
        )))
        tri = triangulate(pts; boundary_nodes=[[[1, 2, 3, 4, 1]], [[5, 7, 6, 5]]])
        @test validate_triangulation(tri)
        validate_statistics(tri)

        tri = triangulate(pts, boundary_nodes=[[[1, 2, 3, 4, 1]], [[5, 7], [7, 6, 5]]], delete_holes=false, delete_ghosts=false)
        points_to_process = DT.find_all_points_to_delete(tri)
        @test points_to_process == Set((DT.𝒢,))
        triangles_to_delete = DT.find_all_triangles_to_delete(tri, points_to_process)
        @test DT.compare_triangle_collections(triangles_to_delete, Set((
            (1, 4, DT.𝒢),
            (4, 3, DT.𝒢),
            (3, 2, DT.𝒢),
            (2, 1, DT.𝒢),
            (5, 6, 7)
        )))
        tri = triangulate(pts; boundary_nodes=[[[1, 2, 3, 4, 1]], [[5, 7], [7, 6, 5]]])
        @test validate_triangulation(tri)
        validate_statistics(tri)

        tri = triangulate(pts, boundary_nodes=[[[1, 2, 3, 4, 1]], [[5, 7, 6, 5]]], delete_holes=false, delete_ghosts=false)
        points_to_process = DT.find_all_points_to_delete(tri)
        @test points_to_process == Set((DT.𝒢,))
        triangles_to_delete = DT.find_all_triangles_to_delete(tri, points_to_process)
        @test DT.compare_triangle_collections(triangles_to_delete, Set((
            (1, 4, DT.𝒢),
            (4, 3, DT.𝒢),
            (3, 2, DT.𝒢),
            (2, 1, DT.𝒢),
            (5, 6, 7)
        )))
        tri = triangulate(pts; boundary_nodes=[[[1, 2, 3, 4, 1]], [[5, 7, 6, 5]]])
        @test validate_triangulation(tri)
        validate_statistics(tri)
    end
end

@testset "A previously broken example" begin
    a = 4 / 5
    t = LinRange(0, 2π, 6)
    x = @. a * (2cos(t) + cos(2t))
    y = @. a * (2sin(t) - sin(2t))
    points = [2.4 -0.152786404500042 -1.047213595499958 -1.047213595499958 -0.1527864045000426; 0.0 1.051462224238267 1.70130161670408 -1.70130161670408 -1.051462224238268]
    bn_nodes = [1, 2, 3, 4, 5, 1]
    tri = triangulate(points; boundary_nodes=bn_nodes, delete_ghosts=false, delete_holes=false)
    points_to_process = DT.find_all_points_to_delete(tri)
    @test points_to_process == Set((DT.𝒢,))
    triangles_to_delete = DT.find_all_triangles_to_delete(tri, points_to_process)
    true_triangles_to_delete = Set((
        (3, 2, 1),
        (4, 1, 5),
        (3, 1, DT.𝒢),
        (1, 4, DT.𝒢),
        (4, 3, DT.𝒢)
    ))
    @test DT.compare_triangle_collections(triangles_to_delete, true_triangles_to_delete)

    tri = triangulate(points; boundary_nodes=bn_nodes)
    @test validate_triangulation(tri)
    validate_statistics(tri)

    a = 4 / 5
    t = LinRange(0, 2π, 100)
    x = @. a * (2cos(t) + cos(2t))
    y = @. a * (2sin(t) - sin(2t)) 
    points = [2.4 2.390342532619651 2.361486661646358 2.313780262033189 2.247797423889998 2.16432999441036 2.064375923309772 1.949124594926212 1.819939377400058 1.678337662931399 1.525968712312314 1.364589651123909 1.196039993627106 1.022215093005981 0.8450389328830836 0.6664366846599798 0.4883074580921845 0.3124976685434042 0.1407754336467903 -0.02519360519360872 -0.1838686646186263 -0.333854062108855 -0.473917012361638 -0.6030027096667079 -0.7202464727404799 -0.824982776042241 -0.9167510419119772 -0.9952981201288541 -1.060577434849961 -1.112744832492037 -1.152151217102105 -1.179332111274063 -1.194994329879372 -1.2 -1.19534820274028 -1.182154550372295 -1.161629044931084 -1.13505259139796 -1.103752559560399 -1.069077803180918 -1.032373553014235 -0.994956601357817 -0.958091190190919 -0.9229660026456057 -0.8906726387622712 -0.8621859315182511 -0.8383464283873544 -0.8198453276894746 -0.8072121183067701 -0.8008051266355886 -0.8008051266355887 -0.8072121183067701 -0.8198453276894746 -0.8383464283873544 -0.8621859315182508 -0.890672638762271 -0.9229660026456055 -0.9580911901909194 -0.9949566013578169 -1.032373553014235 -1.069077803180918 -1.103752559560399 -1.135052591397961 -1.161629044931084 -1.182154550372295 -1.19534820274028 -1.2 -1.194994329879372 -1.179332111274063 -1.152151217102105 -1.112744832492038 -1.060577434849961 -0.9952981201288543 -0.9167510419119774 -0.8249827760422405 -0.7202464727404799 -0.6030027096667084 -0.473917012361639 -0.3338540621088567 -0.1838686646186263 -0.02519360519360934 0.1407754336467914 0.3124976685434047 0.4883074580921845 0.6664366846599791 0.8450389328830821 1.022215093005981 1.196039993627106 1.364589651123908 1.525968712312312 1.678337662931397 1.819939377400058 1.949124594926211 2.064375923309773 2.16432999441036 2.247797423889998 2.313780262033188 2.361486661646357 2.390342532619651 0.3170611563256605 0.3222593190661471 -0.6183519451946706 -0.6497375598908914 -0.03601883360826612 0.6813132651727021 -0.6442990922286358 0.677301056434819 -0.03212472989883586 1.049402030184789 -0.6932272745230642 -0.3561747556617257 1.050855112561903 -0.3586613134669976 -0.6921937990949052 -0.7738275096391242 1.420372082515759 -0.6465445728766354 -0.7401785386842585 -0.8130282799705698 1.553206818654827 -0.4939148325245979 -0.5486837455157717 -0.3610262116797264 -0.73633174592415 -0.5446205627063301 -0.5160523518705824 1.230246578448748 1.093304308222101 0.8980444010246713 1.226019949979472 -0.4997667391156656 -0.726253210863807 -0.6697924507551998 -0.4891744968539103 -0.4634488008923348 -0.301600826224346 -0.2783275986241442 -0.1177692287900377 -0.09393337228059367 -0.2546427830837441 -0.0703514373220399 0.08836743861011277 -0.2313237479039875 0.1137973151090641 0.2724811315356115 0.2996278178466897 0.4682688930880384 0.1380153936605681 0.406551246715369 0.6046552182453121 -0.426784127732563 -0.3923263038801806 -0.208465408072328 -0.3607313375439971 -0.1857892526190982 0.4960419968524275 -0.1415763437967015 0.03591030790050782 0.8685411342185969 -0.1926902902309472 -0.8592155035553608 1.716817787430798 -0.8576022838754369 -0.1921452564055306 -0.6683334786985783 0.8654953794762578 0.6952267131341282 -0.6287806122165458 0.4999252549098672 0.3218569037134379 -0.5986272435711641 0.1449039324769594 2.203335189185611 -1.101667594592808 -1.101667594592807 -0.0286998603756622 -0.3499716813694415 -0.4450559779573338 -0.04783795748406422 0.09480479796065325 2.114352958860066 -1.054192316972474 -1.054192316972475 -0.3317870997991825 1.39665732177026 -0.6206155967973439 -0.7760417249729167 0.4988340023743217 0.1461550962121355 0.3242397406248529 0.9036042833648704 -0.5300851339677813 -0.3652221764727946 -0.9373851776825108 1.875039225311675 -0.9376540476291644 -0.5066666554551625 0.2368123037070368 -0.02440461564986861 0.150646683235148 1.991415825682697 -0.9957149117466657 -0.9957009139360333 0.5336878756578392 0.7210970221965526 -0.1928809320498361 0.7689122343492278 -0.6467071096582683 1.285560087566004 -0.6388529779077368 -0.1969978949942105 0.1875732054682595 -0.4975220279896986 -0.7186754363850769 1.577166091039188 -0.8584906546541111; 0.0 0.000204308591503799 0.001629535973135443 0.005472026448394285 0.01287939060935179 0.02492716987386698 0.04259671987083449 0.06675468790926908 0.09813443380113843 0.1373197116453414 0.1847308933321656 0.2406139730886964 0.3050325470249727 0.3778629130904994 0.4587923858949879 0.5473208683088024 0.6427656684861544 0.7442694978083142 0.8508115330838226 0.9612213760111374 1.074195695220302 1.188317291935738 1.302076290158998 1.413893116908774 1.522142908049466 1.625180951076655 1.721368758301788 1.809100352482636 1.886828342268937 1.953089366954632 2.006528498919954 2.04592220767042 2.070199511290803 2.078460969082653 2.069995202699299 2.044292671697285 2.001056472471559 1.94020997634528 1.86190117239507 1.766503632611802 1.654614070392519 1.527046517275516 1.384823196404126 1.229162223576608 1.061462317070302 0.8832847449107668 0.6963327821298029 0.5024289901161495 0.3034906647750211 0.1015038293221599 -0.1015038293221605 -0.3034906647750203 -0.50242899011615 -0.6963327821298022 -0.8832847449107647 -1.061462317070301 -1.229162223576608 -1.384823196404126 -1.527046517275515 -1.654614070392518 -1.766503632611802 -1.86190117239507 -1.940209976345281 -2.001056472471559 -2.044292671697285 -2.069995202699299 -2.078460969082653 -2.070199511290803 -2.04592220767042 -2.006528498919954 -1.953089366954633 -1.886828342268938 -1.809100352482636 -1.721368758301788 -1.625180951076654 -1.522142908049467 -1.413893116908775 -1.302076290158999 -1.18831729193574 -1.074195695220302 -0.9612213760111378 -0.8508115330838222 -0.7442694978083139 -0.6427656684861544 -0.5473208683088027 -0.4587923858949887 -0.3778629130904991 -0.3050325470249729 -0.2406139730886968 -0.184730893332166 -0.137319711645342 -0.09813443380113851 -0.06675468790926926 -0.04259671987083431 -0.02492716987386698 -0.01287939060935179 -0.005472026448394374 -0.001629535973135488 -0.000204308591503799 -0.5473096449852685 0.5447620596126443 0.01428450391194042 -0.4309321583994598 0.7536074241825893 0.3536290172290472 0.4171992061645116 -0.354274248521617 -0.7697113089269275 0.1945973624956351 0.8115101356751746 -1.00610749817081 -0.1925650703607583 1.006349758355742 -0.8137846879949843 1.19333488738332 0.07348683780306886 -1.266821725186389 1.366146471157314 -1.324086653415268 -0.04174935346989145 1.135404524959176 0.9456563651429131 0.8229544359251958 -0.995445054751826 -0.9480022448317108 -0.741992356781049 -0.13995947020735 0.002345879688797254 -0.08111410706343182 0.130762025431598 -1.127145434944549 0.9963834095129508 0.6153059195085407 0.5329795234649529 0.3351323238863241 0.4528041803207281 0.2552428653861011 0.3751002341424148 0.1760741054849401 0.05622243244883154 -0.02310920543998908 0.0916097709503603 -0.1431517599572751 -0.1033611549989752 0.01112539166303821 -0.1796378655622734 -0.06423821903495724 -0.2933071885887684 0.1346998725151015 0.06389662809144629 -0.06697977586584479 -0.263147702667072 -0.3418086892373844 -0.4602230597841197 -0.5395791962871652 -0.260267079463506 0.5739282454666178 0.5007002518322551 0.2757413735876992 -0.8939205770864943 -1.486342121175829 0.0009313928164653013 1.487273513992295 0.8909822595409715 -0.6270118467477362 -0.2692763527002932 -0.1472140205276644 0.2179450297267495 0.4510351900656056 0.3329413377177828 -0.2182665505265999 -0.6631900778915942 0.0 1.908144246886934 -1.908144246886936 -0.4106426344765784 -0.6482825657476072 0.1353078965823503 -0.2192301054498285 0.3004465131589842 4.336808689942018e-17 1.825914653945082 -1.825914653945081 0.6415639095768699 -0.08968357289590935 1.254408379505216 -1.164673062563948 -0.4517076066327135 0.6541117112279523 -0.36426907951982 0.09997393244042915 0.731016147780223 -0.8325312305826258 1.623909218280243 0.0001552321360775167 -1.623753986144165 -0.5474230356714256 0.1741645772010534 -0.5950229828132197 -0.4758032489975526 8.081639737527319e-6 -1.724612653719712 1.724620735359448 0.2655807939077015 0.1867114129881211 -0.717357771233564 0.01594342553679897 -1.111060401383298 -0.004534585080420719 1.115594986463718 0.7308601982314702 0.4541306336569979 -0.3789760643386818 -1.406227077779554 0.08072235390443261 1.325504723875121]
    bn_nodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 1]
    tri = triangulate(points; boundary_nodes=bn_nodes)
    @test validate_triangulation(tri)
    validate_statistics(tri)
end