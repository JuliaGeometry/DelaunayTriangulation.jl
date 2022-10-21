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

@testset "Can we correctly access a set of points at the boundary index?" begin
    pts = rand(SVector{2,Float64}, 1827301)
    for i in shuffle(Int64[zeros(5000)..., (1:length(pts))...])
        if i == 0
            @test DT.get_point(pts, i) == (DT.CentroidCoordinates.x, DT.CentroidCoordinates.y)
        else
            @test DT.get_point(pts, i) == pts[i]
        end
    end
end

@testset "Can we correctly compute the centroid of some given points?" begin 
    DT.CentroidCoordinates.x = 0.291
    DT.CentroidCoordinates.y = 0.4991
    DT.CentroidCoordinates.n = 99
    pts = rand(SVector{2,Float64}, 1827301)
    DT.compute_centroid!(pts)
    @test DT.CentroidCoordinates.x ≈ sum(pts)[1]/length(pts)
    @test DT.CentroidCoordinates.y ≈ sum(pts)[2]/length(pts)
    @test DT.CentroidCoordinates.n == length(pts)
end