"""
    is_exterior_curve(tri::Triangulation, curve_index) -> Bool 

Returns `true` if the `curve_index`th curve in `tri` is an exterior curve, and `false` otherwise.
"""
is_exterior_curve(tri::Triangulation, curve_index) = curve_index ∈ get_exterior_curve_indices(tri)

"""
    is_interior_curve(tri::Triangulation, curve_index) -> Bool

Returns `true` if the `curve_index`th curve in `tri` is an interior curve, and `false` otherwise.
"""
is_interior_curve(tri::Triangulation, curve_index) = !is_exterior_curve(tri, curve_index)

"""
    num_exterior_curves(tri::Triangulation) -> Integer

Returns the number of exterior curves in `tri`.
"""
num_exterior_curves(tri::Triangulation) = (length ∘ get_exterior_curve_indices)(tri)

"""
    is_disjoint(tri::Triangulation) -> Bool

Returns `true` if `tri` has disjoint exterior boundary curves, and `false` otherwise.
"""
is_disjoint(tri::Triangulation) = num_exterior_curves(tri) > 1