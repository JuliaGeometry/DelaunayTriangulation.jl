############################################
##
## CONSTANTS 
##
############################################
"""The index used to represent the lower-right vertex of the bounding triangle, p₋₁."""
const LowerRightBoundingIndex = -1
"""The index used to represent the upper vertex of the bounding triangle, p₋₂."""
const UpperBoundingIndex = -2
"""The index used to represent the lower-left of the bounding triangle, p₋₃."""
const LowerLeftBoundingIndex = -3
"""The triangle representing the bounding triangle for any triangulation."""
const BoundingTriangle = (LowerRightBoundingIndex,
    UpperBoundingIndex,
    LowerLeftBoundingIndex)
"""The index used to represent the unbounded part of the triangulation, i.e. the boundary."""
const BoundaryIndex = 0
"""The first index where physical points begin."""
const FirstPointIndex = 1
"""The multiplier for shifting the bounding triangle."""
const BoundingTriangleShift = 27.3919105998678246734637673334353499992
"""The minimum width and height of the bounding box."""
const MinWidthHeight = 0.00182829999105935935923828882837273741
"""Representation of an empty set for the default value in Adjacent."""
const DefaultAdjacentValue = -4