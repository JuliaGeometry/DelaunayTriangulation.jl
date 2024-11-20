@enum Certificate begin
    Inside
    Degenerate
    Outside
    On
    Left
    Right
    PositivelyOriented
    NegativelyOriented
    Collinear
    None
    Single
    Multiple
    Touching
    Legal
    Illegal
    Closer
    Further
    Equidistant
    Obtuse
    Acute
    SuccessfulInsertion
    FailedInsertion
    PrecisionFailure
    EncroachmentFailure
    Above
    Below
    Visible
    Invisible
end

gen_doc(name, cert, prefix = "is") = """
    $(prefix)_$name(cert::Certificate) -> Bool

Returns `true` if `cert == $cert`, and `false` otherwise.
"""

@doc gen_doc("inside", Inside) is_inside(cert::Certificate) = cert == Inside
@doc gen_doc("degenerate", Degenerate) is_degenerate(cert::Certificate) = cert == Degenerate
@doc gen_doc("outside", Outside) is_outside(cert::Certificate) = cert == Outside
@doc gen_doc("on", On) is_on(cert::Certificate) = cert == On
@doc gen_doc("left", Left) is_left(cert::Certificate) = cert == Left
@doc gen_doc("right", Right) is_right(cert::Certificate) = cert == Right
@doc gen_doc("positively_oriented", PositivelyOriented) is_positively_oriented(cert::Certificate) = cert == PositivelyOriented
@doc gen_doc("negatively_oriented", NegativelyOriented) is_negatively_oriented(cert::Certificate) = cert == NegativelyOriented
@doc gen_doc("collinear", Collinear) is_collinear(cert::Certificate) = cert == Collinear
@doc gen_doc("none", None) is_none(cert::Certificate) = cert == None
@doc gen_doc("single", Single) is_single(cert::Certificate) = cert == Single
@doc gen_doc("multiple", Multiple) is_multiple(cert::Certificate) = cert == Multiple
@doc gen_doc("touching", Touching) is_touching(cert::Certificate) = cert == Touching
@doc gen_doc("legal", Legal) is_legal(cert::Certificate) = cert == Legal
@doc gen_doc("illegal", Illegal) is_illegal(cert::Certificate) = cert == Illegal
@doc gen_doc("closer", Closer) is_closer(cert::Certificate) = cert == Closer
@doc gen_doc("further", Further) is_further(cert::Certificate) = cert == Further
@doc gen_doc("equidistant", Equidistant) is_equidistant(cert::Certificate) = cert == Equidistant
@doc gen_doc("obtuse", Obtuse) is_obtuse(cert::Certificate) = cert == Obtuse
@doc gen_doc("acute", Acute) is_acute(cert::Certificate) = cert == Acute
@doc gen_doc("successful_insertion", SuccessfulInsertion) is_successful_insertion(cert::Certificate) = cert == SuccessfulInsertion
@doc gen_doc("failed_insertion", FailedInsertion) is_failed_insertion(cert::Certificate) = cert == FailedInsertion
@doc gen_doc("precision_failure", PrecisionFailure) is_precision_failure(cert::Certificate) = cert == PrecisionFailure
@doc gen_doc("encroachment_failure", EncroachmentFailure) is_encroachment_failure(cert::Certificate) = cert == EncroachmentFailure
@doc gen_doc("above", Above) is_above(cert::Certificate) = cert == Above
@doc gen_doc("below", Below) is_below(cert::Certificate) = cert == Below
@doc gen_doc("visible", Visible) is_visible(cert::Certificate) = cert == Visible
@doc gen_doc("invisible", Invisible) is_invisible(cert::Certificate) = cert == Invisible
@doc gen_doc("no_intersections", None, "has") has_no_intersections(cert::Certificate) = is_none(cert)
@doc gen_doc("one_intersection", Single, "has") has_one_intersection(cert::Certificate) = is_single(cert)
@doc gen_doc("multiple_intersections", Multiple, "has") has_multiple_intersections(cert::Certificate) = is_multiple(cert)

function convert_certificate(cert, cert1::Certificate, cert2::Certificate, cert3::Certificate)::Certificate
    return cert < 0 ? cert1 : cert > 0 ? cert3 : cert2
end