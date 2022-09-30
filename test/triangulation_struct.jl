@testset "Can we correctly add and delete triangles?" begin
    𝒯 = DT.Triangles(DT.TriangleType[(1, 2, 3), (3, 2, 5), (2, 1, 9)])
    DT.delete_triangle!(𝒯, (1, 2, 3))
    @test 𝒯.triangles == DT.TriangleType[(3, 2, 5), (2, 1, 9)] |> Set
    DT.delete_triangle!(𝒯, (5, 3, 2))
    @test 𝒯.triangles == DT.TriangleType[(2, 1, 9)] |> Set
    DT.add_triangle!(𝒯, (2, 3, 10))
    @test 𝒯.triangles == DT.TriangleType[(2, 1, 9), (2, 3, 10)] |> Set
    DT.add_triangle!(𝒯, (2, 3, 11), (11, 3, 4), (2, 3, 1))
    @test 𝒯.triangles == DT.TriangleType[(2, 1, 9), (2, 3, 10), (2, 3, 11), (11, 3, 4), (2, 3, 1)] |> Set
    DT.delete_triangle!(𝒯, (2, 1, 9), (10, 2, 3))
    @test 𝒯.triangles == DT.TriangleType[(2, 3, 11), (11, 3, 4), (2, 3, 1)] |> Set
end