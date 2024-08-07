using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using ..DelaunayTriangulation: Certificate

global i = [
    Certificate.Inside,
    Certificate.Degenerate,
    Certificate.Outside,
    Certificate.On,
    Certificate.Left,
    Certificate.Right,
    Certificate.PositivelyOriented,
    Certificate.NegativelyOriented,
    Certificate.Collinear,
    Certificate.None,
    Certificate.Single,
    Certificate.Multiple,
    Certificate.Touching,
    Certificate.Legal,
    Certificate.Illegal,
    Certificate.Closer,
    Certificate.Further,
    Certificate.Equidistant,
    Certificate.Above,
    Certificate.Below,
    Certificate.EncroachmentFailure,
    Certificate.PrecisionFailure,
    Certificate.Obtuse,
    Certificate.Acute,
    Certificate.FailedInsertion,
    Certificate.Visible,
    Certificate.Invisible,
]
global j = [
    DT.is_inside,
    DT.is_degenerate,
    DT.is_outside,
    DT.is_on,
    DT.is_left,
    DT.is_right,
    DT.is_positively_oriented,
    DT.is_negatively_oriented,
    DT.is_collinear,
    DT.has_no_intersections,
    DT.has_one_intersection,
    DT.has_multiple_intersections,
    DT.is_touching,
    DT.is_legal,
    DT.is_illegal,
    DT.is_closer,
    DT.is_further,
    DT.is_equidistant,
    DT.is_above,
    DT.is_below,
    DT.is_encroachment_failure,
    DT.is_precision_failure,
    DT.is_obtuse,
    DT.is_acute,
    DT.is_failed_insertion,
    DT.is_visible,
    DT.is_invisible,
]

@testset "Classifiers" begin
    for a in eachindex(i)
        for b in eachindex(j)
            if a == b
                @test j[b](i[a])
            else
                @test !j[b](i[a])
            end
        end
    end
end

@testset "Conversion" begin
    for _ in 1:100
        certs = rand(i, 3)
        k = rand(1:3)
        @test DT.convert_certificate((-1, 0, 1)[k], certs...) == certs[k]
    end
end
