"""
    cross(p::SphericalPoint, q::SphericalPoint) -> NTuple{3, Number}

Computes the cross product of two points on the sphere, treating them 
as vectors from the origin to the points. The result is returned 
as an `NTuple{3, Number}`.
"""
function cross(p::SphericalPoint, q::SphericalPoint)
    _px, _py, _pz = __getxyz(p)
    _qx, _qy, _qz = __getxyz(q)
    return (
        _py * _qz - _pz * _qy,
        _pz * _qx - _px * _qz,
        _px * _qy - _py * _qx
    )
end

_hypot(x, y) = hypot(x, y)
_hypot(x, y, z) = hypot(hypot(x, y), z)

function dot(p::SphericalPoint, q::SphericalPoint)
    _px, _py, _pz = __getxyz(p)
    _qx, _qy, _qz = __getxyz(q)
    return _px * _qx + _py * _qy + _pz * _qz
end

function scalar_triple_product(p::SphericalPoint, q::SphericalPoint, r::SphericalPoint)
    return dot(p, SphericalPoint(cross(q, r)))
end

function perturbed_north_pole(::Type{T}) where {T}
    dx, dy = 1e-8, 1e-8
    return SphericalPoint(T(dx), T(dy), T(sqrt(1 - dx^2 - dy^2)))
end

function fix_north_pole(p, q, r)
    T = number_type(p)
    np = north_pole(T)
    if p == np
        return perturbed_north_pole(T), q, r
    elseif q == np
        return p, perturbed_north_pole(T), r
    elseif r == np
        return p, q, perturbed_north_pole(T)
    else
        return p, q, r
    end
end