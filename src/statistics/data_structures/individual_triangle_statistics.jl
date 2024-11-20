# Some formulas in this file come from https://perso.uclouvain.be/jean-francois.remacle/LMECA2170/robnotes.pdf.

struct IndividualTriangleStatistics{T}
    area::T
    lengths::NTuple{3, T}
    circumcenter::NTuple{2, T}
    circumradius::T
    angles::NTuple{3, T}
    radius_edge_ratio::T
    edge_midpoints::NTuple{3, NTuple{2, T}}
    aspect_ratio::T
    inradius::T
    perimeter::T
    centroid::NTuple{2, T}
    offcenter::NTuple{2, T}
    sink::NTuple{2, T}
end
function IndividualTriangleStatistics(p, q, r, sink = (NaN, NaN))
    F = number_type(p)
    ℓmin², ℓmed², ℓmax² = squared_triangle_lengths(p, q, r)
    ℓmin, ℓmed, ℓmax = sqrt(ℓmin²), sqrt(ℓmed²), sqrt(ℓmax²)
    A = triangle_area(p, q, r)
    circumcenter = triangle_circumcenter(p, q, r, A)
    circumradius = triangle_circumradius(A, ℓmin², ℓmed², ℓmax²)
    radius_edge_ratio = triangle_radius_edge_ratio(circumradius, ℓmin)
    edge_midpoints = triangle_edge_midpoints(p, q, r)
    perimeter = triangle_perimeter(ℓmin, ℓmed, ℓmax)
    inradius = triangle_inradius(A, perimeter)
    aspect_ratio = triangle_aspect_ratio(inradius, circumradius)
    centroid = triangle_centroid(p, q, r)
    angles = triangle_angles(p, q, r)
    offcenter = triangle_offcenter(p, q, r, circumcenter)
    sx, sy = getxy(sink)
    return IndividualTriangleStatistics(
        F(A),
        F.((ℓmin, ℓmed, ℓmax)),
        F.(circumcenter),
        F(circumradius),
        F.(angles),
        F(radius_edge_ratio),
        (F.(edge_midpoints[1]), F.(edge_midpoints[2]), F.(edge_midpoints[3])),
        F(aspect_ratio),
        F(inradius),
        F(perimeter),
        F.(centroid),
        F.(offcenter),
        (F(sx), F(sy)),
    )
end