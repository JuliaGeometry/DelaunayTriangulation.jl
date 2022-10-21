@testset "Can we correctly reset the centroid?" begin
    DT.CentroidCoordinates.x = 0.291
    DT.CentroidCoordinates.y = 0.4991
    DT.CentroidCoordinates.n = 99
    DT.clear_centroid_coordinates!()
    @test DT.CentroidCoordinates.x == 0.0
    @test DT.CentroidCoordinates.y == 0.0
    @test DT.CentroidCoordinates.n == 0
end

@testset "Can we correctly update the centroid incrementally?" begin
    pts = rand(SVector{2,Float64}, 9371)
    idx = shuffle(1:length(pts))
    DT.clear_centroid_coordinates!()
    for (j, i) in enumerate(idx)
        DT.update_centroid_after_new_point!(pts, i)
        pxy = sum(pts[idx[1:j]]) / j
        @test DT.CentroidCoordinates.x ≈ pxy[1]
        @test DT.CentroidCoordinates.y ≈ pxy[2]
        @test DT.CentroidCoordinates.n == j
    end
    idx = shuffle(1:(length(pts)-1)) # Don't want to go down to dividing by zero
    for (j, i) in enumerate(idx)
        DT.update_centroid_after_deleted_point!(pts, i)
        pxy = (sum(pts) - sum(pts[idx[1:j]])) / (length(pts) - j)
        @test DT.CentroidCoordinates.x ≈ pxy[1]
        @test DT.CentroidCoordinates.y ≈ pxy[2]
        @test DT.CentroidCoordinates.n == length(pts) - j
    end
end