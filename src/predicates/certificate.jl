"""
    Certificate

An `Enum` type that represents results from a geometric predicate. Below we provide a list of available certificates, 
along with the function that can be used for testing if a given `Certificate` matches that `certificate`.

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
- `Right`: `is_right`
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
end

for inst in instances(Certificate.T)
    name = String(Symbol(inst))
    @eval ($(Symbol("is_$(lowercase(name))")))(cert::Certificate.T) = cert == $inst
end
is_positively_oriented(cert::Certificate.T) = is_positivelyoriented(cert)
is_negatively_oriented(cert::Certificate.T) = is_negativelyoriented(cert)
has_no_intersections(cert::Certificate.T) = is_none(cert)
has_one_intersection(cert::Certificate.T) = is_single(cert)
has_multiple_intersections(cert::Certificate.T) = is_multiple(cert)

@inline function convert_certificate(cert::I, Cert1, Cert2, Cert3)::Certificate.T where {I}
    if cert == I(-1)
        return Cert1
    elseif cert == I(0)
        return Cert2
    elseif cert == I(1)
        return Cert3
    end
    throw(ArgumentError("The provided certificate value, $cert, must be one of ($(I(-1)), $(I(0)), $(I(1))."))
end

const Cert = Certificate
