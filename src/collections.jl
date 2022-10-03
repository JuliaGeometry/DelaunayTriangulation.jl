############################################
##
## PRIMITIVE COLLECTIONS 
##
############################################
"""
    struct Triangles{I,T<:AbstractTriangle{I}}

Structure for a collection of `Triangle`s, currently implemented using a `Set`. 
The triangles can be accessed via [`triangles`](@ref).
"""
struct Triangles{I,T<:AbstractTriangle{I}}
    triangles::Set{T}
    Triangles(triangles::Set{T}) where {I,T<:AbstractTriangle{I}} = new{I,T}(triangles)
    Triangles{I,T}(triangles::Set{F}) where {I,T,F} = new{I,T}(convert(Set{T}, triangles))
    Triangles(triangles::Base.AbstractVecOrTuple{T}) where {I,T<:AbstractTriangle{I}} = new{I,T}(Set{T}(triangles))
    Triangles(triangles...) = Triangles(triangles)
end
triangles(tris::Triangles) = tris.triangles
Base.iterate(T::Triangles) = Base.iterate(triangles(T))
Base.iterate(T::Triangles, state) = Base.iterate(triangles(T), state)
Base.eltype(::Triangles{I,T}) where {I,T<:AbstractTriangle{I}} = T
Base.eltype(::Type{Triangles{I,T}}) where {I,T<:AbstractTriangle{I}} = T
Base.length(T::Triangles) = length(T.triangles)
Base.enumerate(T::Triangles) = Base.enumerate(triangles(T))

"""
    add_triangle!(T::Triangles, V::AbstractTriangle...)

Adds the triangles in `V` into the collection of triangles `T`. No 
checks are made for triangles that are equal under circular shifts.
"""
function add_triangle!(T::Triangles, V::AbstractTriangle...)
    push!(triangles(T), V...)
    return nothing
end

"""
    delete_triangle!(T::Triangles, V::AbstractTriangle...)

Deletes the triangle `V` from the collection of triangles `T`. All circular 
shifts of `V`'s indices are also deleted.
"""
function delete_triangle!(T::Triangles, V::AbstractTriangle)
    delete!(triangles(T), V)
    delete!(triangles(T), shift_triangle_1(V))
    delete!(triangles(T), shift_triangle_2(V))
    return nothing
end
@doc (@doc delete_triangle!(::Triangles, ::AbstractTriangle))
function delete_triangle!(T::Triangles, V::AbstractTriangle...)
    for V in V
        delete_triangle!(T, V)
    end
    return nothing
end

"""
    Points{T,A,P<:AbstractPoint{T,A}}

Structure for a collection of `Point`s, currently implemented using a `Vector`.
The points can be accessed using [`points`](@ref). The struct contains the following 
fields: 

- `points::Vector{P}`: The vector of points. 
- `xcentroid::T`: The `x`-coordinate of the centroid, `(xmax + xmin)/2`.
- `ycentroid::T`: The `y`-coordinate of the centroid, `(ymax + ymin)/2`.
- `xmin::T`: The minimum `x`-coordinate in the collection of points.
- `xmax::T`: The maximum `x`-coordinate in the collection of points. 
- `ymin::T`: The minimum `y`-coordinate in the collection of points. 
- `ymax::T`: The maximum `y`-coordinate in the collection of points.
- `width::T`: The width of the collection of points, defined by `xmax - xmin`.
- `height::T`: The height of the collection of points, defined by `ymax - ymin`.
- `max_width_height`: The maximum of `width` and `height`.
- `lower_left_bounding_triangle_coords`: The coordinate of the lower-left vertex of the bounding triangle. 
- `lower_right_bounding_triangle_coords`: The coordinate of the lower-right vertex of the bounding triangle. 
- `upper_bounding_triangle_coords`: The coordinate of the upper vertex of the bounding triangle. 

(Strictly speaking, the 'centroid' is the centre of the box that tightly bounds the points.)

If you wish to provide specific values for the bounding triangle's coordinates, then the constructor `Points` 
accepts keyword arguments `pLR`, `pLL`, and `pU` for the lower-right, lower-left, and upper vertices, respectively.
"""
struct Points{T,A,P<:AbstractPoint{T,A}}
    points::Vector{P}
    xcentroid::T
    ycentroid::T
    xmin::T
    xmax::T
    ymin::T
    ymax::T
    width::T
    height::T
    max_width_height::T
    lower_left_bounding_triangle_coords::P
    lower_right_bounding_triangle_coords::P
    upper_bounding_triangle_coords::P
    function Points(points::AbstractVector{P};
        pLL=nothing, pLR=nothing, pU=nothing) where {T,A,P<:AbstractPoint{T,A}}
        xmin = typemax(T)
        xmax = typemin(T)
        ymin = typemax(T)
        ymax = typemin(T)
        for pt in points
            if getx(pt) < xmin
                xmin = getx(pt)
            end
            if getx(pt) > xmax
                xmax = getx(pt)
            end
            if gety(pt) < ymin
                ymin = gety(pt)
            end
            if gety(pt) > ymax
                ymax = gety(pt)
            end
        end
        width = max(xmax - xmin, MinWidthHeight)
        height = max(ymax - ymin, MinWidthHeight)
        xcentroid = (xmax + xmin) / 2
        ycentroid = (ymax + ymin) / 2
        max_width_height = max(width, height)
        lower_left_bounding_triangle_coords = isnothing(pLL) ? P(xcentroid - BoundingTriangleShift * max_width_height,
            ycentroid - max_width_height) : pLL
        lower_right_bounding_triangle_coords = isnothing(pLR) ? P(xcentroid + BoundingTriangleShift * max_width_height,
            ycentroid - max_width_height) : pLR
        upper_bounding_triangle_coords = isnothing(pU) ? P(xcentroid, ycentroid + BoundingTriangleShift * max_width_height) : pU
        return new{T,A,P}(points, xcentroid, ycentroid, xmin, xmax, ymin,
            ymax, width, height, max_width_height, lower_left_bounding_triangle_coords,
            lower_right_bounding_triangle_coords, upper_bounding_triangle_coords)
    end
    Points(points::NTuple{N,P};
        pLL=nothing, pLR=nothing, pU=nothing) where {N,T,A,P<:AbstractPoint{T,A}} = Points(collect(points); pLL, pLR, pU)
    Points{T,A,P}(points, xcentroid, ycentroid, xmin, xmax, ymin, ymax, width, height,
        max_width_height, lower_left_bounding_triangle_coords, lower_right_bounding_triangle_coords,
        upper_bounding_triangle_coords) where {T,A,P} = new{T,A,P}(points, xcentroid, ycentroid, xmin, xmax, ymin, ymax, width, height,
        max_width_height, lower_left_bounding_triangle_coords, lower_right_bounding_triangle_coords,
        upper_bounding_triangle_coords)
end
points(pts::Points) = pts.points
Base.length(pts::Points) = length(points(pts))
Points(pts::AbstractVector; pLL=nothing, pLR=nothing, pU=nothing) = Points(Point.(pts); pLL, pLR, pU)
Points(pts::Points) = pts
Points(pts...; pLL=nothing, pLR=nothing, pU=nothing) = Points(collect(pts); pLL, pLR, pU)
Base.iterate(pts::Points) = Base.iterate(points(pts))
Base.iterate(pts::Points, state) = Base.iterate(points(pts), state)
Base.eltype(::Points{T,A,P}) where {T,A,P<:AbstractPoint{T,A}} = P
Base.eltype(::Type{Points{T,A,P}}) where {T,A,P<:AbstractPoint{T,A}} = P
width(pts::Points) = pts.width
height(pts::Points) = pts.height
xmin(pts::Points) = pts.xmin
ymin(pts::Points) = pts.ymin
xmax(pts::Points) = pts.xmax
ymax(pts::Points) = pts.ymax
xcentroid(pts::Points) = pts.xcentroid
ycentroid(pts::Points) = pts.ycentroid
max_width_height(pts::Points) = pts.max_width_height
lower_left_bounding_triangle_coords(pts::Points) = pts.lower_left_bounding_triangle_coords
lower_right_bounding_triangle_coords(pts::Points) = pts.lower_right_bounding_triangle_coords
upper_bounding_triangle_coords(pts::Points) = pts.upper_bounding_triangle_coords
Base.size(pts::Points) = Base.size(points(pts))
Base.eachindex(pts::Points) = Base.eachindex(points(pts))
Random.shuffle!(pts::Points) = Random.shuffle!(points(pts))
Base.lastindex(pts::Points) = Base.lastindex(points(pts))
Base.firstindex(pts::Points) = Base.firstindex(points(pts))


"""
    get_point(pts::Points, i)

Obtains the `i`th point in `pts`. We also allow for indices 

- `i = LowerLeftBoundingIndex`: This returns the coordinates for the lower-left vertex of the bounding triangle.
- `i = LowerRightBoundingIndex`: This returns the coordinates for the lower-right vertex of the bounding triangle.
- `i = UpperBoundingIndex`: This returns the coordinates for the upper vertex of the bounding triangle.

These latter coordinates are stored in `Points` rather than computed lazily.
"""
function get_point(pts::Points, i)
    if i ≥ FirstPointIndex
        return points(pts)[i]
    elseif i == LowerRightBoundingIndex
        return lower_right_bounding_triangle_coords(pts)
    elseif i == LowerLeftBoundingIndex
        return lower_left_bounding_triangle_coords(pts)
    elseif i == UpperBoundingIndex
        return upper_bounding_triangle_coords(pts)
    end
    throw(BoundsError(pts, i))
end

"""
    add_point!(pts::Points, p...)
    
Adds the points in `p` to the collection of points `pts`.
"""
function add_point!(pts::Points, p...)
    push!(points(pts), p...)
    return nothing
end