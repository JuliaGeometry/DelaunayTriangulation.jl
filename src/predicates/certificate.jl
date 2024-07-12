"""
    Certificate

This is an `Enum` that defines certificates returned from predicates. The instances, and associated certifiers, are:

- `Inside`: `is_inside`
- `Degenerate`: `is_degenerate`
- `Outside`: `is_outside`
- `On`: `is_on`
- `Left`: `is_left`
- `Right`: `is_right`
- `PositivelyOriented`: `is_positively_oriented`
- `NegativelyOriented`: `is_negatively_oriented`
- `Collinear`: `is_collinear`
- `None`: `is_none` or `has_no_intersections`
- `Single`: `is_single` or `has_one_intersection`
- `Multiple`: `is_multiple` or `has_multiple_intersections`
- `Touching`: `is_touching`
- `Legal`: `is_legal`
- `Illegal`: `is_illegal`
- `Closer`: `is_closer`
- `Further`: `is_further`
- `Equidistant`: `is_equidistant`
- `Obtuse`: `is_obtuse`
- `Acute`: `is_acute`
- `SuccessfulInsertion`: `is_successful_insertion`
- `FailedInsertion`: `is_failed_insertion`
- `PrecisionFailure`: `is_precision_failure`
- `EncroachmentFailure`: `is_encroachment_failure`
- `Above`: `is_above`
- `Below`: `is_below`
- `Visible`: `is_visible`
- `Invisible`: `is_invisible`
"""
EnumX.@enumx Certificate begin
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


for inst in instances(Certificate.T)
    name = String(Symbol(inst))
    @eval begin
        @doc """
            is_$($(lowercase(name)))(cert::Certificate) -> Bool

        Returns `true` if `cert` is `$($name)`, and `false` otherwise.
        """ ($(Symbol("is_$(lowercase(name))")))(cert::Certificate.T) = cert == $inst
    end
end
"""
    is_positively_oriented(cert::Certificate) -> Bool

Returns `true` if `cert` is `PositivelyOriented`, and `false` otherwise.
"""
is_positively_oriented(cert::Certificate.T) = is_positivelyoriented(cert)


"""
    is_negatively_oriented(cert::Certificate) -> Bool

Returns `true` if `cert` is `NegativelyOriented`, and `false` otherwise.
"""
is_negatively_oriented(cert::Certificate.T) = is_negativelyoriented(cert)


"""
    has_no_intersections(cert::Certificate) -> Bool

Returns `true` if `cert` is `None`, and `false` otherwise. Synonymous with `is_none`.
"""
has_no_intersections(cert::Certificate.T) = is_none(cert)


"""
    has_one_intersection(cert::Certificate) -> Bool

Returns `true` if `cert` is `Single`, and `false` otherwise. Synonymous with `is_single`.
"""
has_one_intersection(cert::Certificate.T) = is_single(cert)


"""
    has_multiple_intersections(cert::Certificate) -> Bool

Returns `true` if `cert` is `Multiple`, and `false` otherwise. Synonymous with `is_multiple`.
"""
has_multiple_intersections(cert::Certificate.T) = is_multiple(cert)
is_successful_insertion(cert::Certificate.T) = is_successfulinsertion(cert)
is_failed_insertion(cert::Certificate.T) = is_failedinsertion(cert)
is_precision_failure(cert::Certificate.T) = is_precisionfailure(cert)
is_encroachment_failure(cert::Certificate.T) = is_encroachmentfailure(cert)


"""
    convert_certificate(cert::I, Cert1, Cert2, Cert3) -> Certificate 

Given `cert ∈ (-1, 0, 1)`, return `Cert1`, `Cert2` or `Cert3` depending on if `cert == -1`,
`cert == 0` or `cert == 1`, respectively.
"""
@inline function convert_certificate(cert::I, Cert1, Cert2, Cert3)::Certificate.T where {I}
    if cert == I(-1)
        return Cert1
    elseif cert == I(0)
        return Cert2
    else # if cert == I(1)
        return Cert3
    end
end


const Cert = Certificate