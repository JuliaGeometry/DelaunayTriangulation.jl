"""
    SphericalLine(p, q) 

A line on the surface of a sphere defined by the points `p` and `q`. The line is parametrised by the angle `θ` from `p` to `q`.
You can evaluate this struct like a function.
"""
struct SphericalLine{P,T} <: AbstractParametricCurve
    p::SphericalPoint{P}
    q::SphericalPoint{P}
    upper::T
    orth::NTuple{3,T}
end
function SphericalLine(p::SphericalPoint, q::SphericalPoint)
    # https://comp.soft-sys.matlab.narkive.com/lSloCHDO/line-between-two-points-on-a-sphere
    pq_cross = cross(p, q)
    pq_orth = cross(SphericalPoint(pq_cross), p)
    pq_orth_norm = _hypot(pq_orth...)
    pq_orth = pq_orth ./ pq_orth_norm
    pq_dot = dot(p, q)
    pq_cross_norm = _hypot(pq_cross...)
    ϕ = mod2pi(atan(pq_cross_norm, pq_dot))
    return SphericalLine(p, q, ϕ, pq_orth)
end
function (L::SphericalLine)(t)
    s, c = sincos(t)
    p = getp(L.p)
    orth = L.orth
    coords = p .* c .+ orth .* s
    return SphericalPoint(coords)
end

"""
    _to_coords(L::SphericalLine[, n = 20]) -> Vector{NTuple{3, Number}}

Returns the coordinates of the points on the line `L` as a vector of `NTuple{3, Number}`.
The line is evaluated at `n` equally spaced points.
"""
function _to_coords(L::SphericalLine, n)
    upper = L.upper
    t = LinRange(zero(upper), upper, n)
    points = L.(t)
    return points
end

"""
    SphericalTriangle(p, q, r)

A triangle on the surface of a sphere defined by the points `p`, `q`, and `r`. The triangle is stored 
as three [`SphericalLine`](@ref)s.
"""
struct SphericalTriangle{P,T}
    L1::SphericalLine{P,T}
    L2::SphericalLine{P,T}
    L3::SphericalLine{P,T}
end
function SphericalTriangle(p::SphericalPoint, q::SphericalPoint, r::SphericalPoint)
    return SphericalTriangle(
        SphericalLine(p, q),
        SphericalLine(q, r),
        SphericalLine(r, p)
    )
end

"""
    _to_coords(T::SphericalTriangle, add_end=true, n=20) -> Vector{NTuple{3, Number}}

Returns the coordinates of the points on the triangle `T` as a vector of `NTuple{3, Number}`. If `add_end` is `true`,
the first point is repeated at the end of the vector.

The triangle is evaluated at `n` equally spaced points on each line.
"""
function _to_coords(T::SphericalTriangle, add_end=true, n=20)
    coords1 = _to_coords(T.L1, n)
    _coords2 = _to_coords(T.L2, n)
    coords2 = view(_coords2, 2:length(_coords2))
    _coords3 = _to_coords(T.L3, n)
    coords3 = view(_coords3, 2:(length(_coords3)-1))
    coords = vcat(coords1, coords2, coords3)
    add_end && push!(coords, coords[1])
    return coords
end

"""
    _to_mesh(points::Vector{<:SphericalPoint}[, offset=0]) -> Matrix{Int}

Given a closed curve on the surface of a [`UnitSphere`](@ref) defined by `points`, returns a 
triangulation of that curve appropriate for use with Makie.jl's `mesh`.

The optional argument `offset` is the index of the first point in the curve.
"""
function _to_mesh(points::Vector{<:SphericalPoint}, offset=0)
    @assert points[begin] != points[end]
    tri = spherical_triangulate(points)
    triangles = Matrix{Int}(undef, num_solid_triangles(tri), 3)
    for (i, T) in enumerate(each_solid_triangle(tri))
        u, v, w = triangle_vertices(T)
        triangles[i, :] .= (u, v, w) .+ offset
    end
    return triangles
end

"""
    SphericalTriangle(tri::Triangulation, T)

Given a Delaunay triangulation `tri` and a triangle `T`, returns the `SphericalTriangle` corresponding to `T`.
"""
function SphericalTriangle(tri::Triangulation, T)
    points = get_points(tri)
    V = number_type(tri)
    T = sort_triangle(T)
    u, v, w = triangle_vertices(T)
    p, q, r = _getindex(points, u), _getindex(points, v), _getindex(points, w) # can't do get_point since DelaunayTriangulation projects it
    p, q, r = fix_north_pole(p, q, r)
    ST = SphericalTriangle(p, q, r)
    return ST
end

"""
    _to_triangles(tri::SphericalTriangulation[, n = 20]) -> Vector{<:SphericalPoint}, Matrix{Int}

Given a spherical triangulation `tri`, returns the points and triangles of the triangulation as a vector of `SphericalPoint`s
and a matrix of indices. This format is appropriate for use with Makie.jl's `mesh`. `n` is the number of points to use for discretising 
each line of each triangle.
"""
function _to_triangles(tri::SphericalTriangulation, n = 20)
    points = get_points(tri)
    mesh_points = eltype(points)[]
    vertices = Matrix{Int}(undef, 0, 3)
    offset = 0
    for T in each_triangle(tri)
        ST = SphericalTriangle(tri, T)
        coords = _to_coords(ST, false, n)
        _vertices = _to_mesh(coords, offset)
        offset += length(coords)
        vertices = vcat(vertices, _vertices)
        append!(mesh_points, coords)
    end
    return mesh_points, vertices
end

"""
    _to_lines(tri::SphericalTriangulation[, n = 20]) -> Vector{<:SphericalPoint}

Given a spherical triangulation `tri`, returns the points on the lines of the triangulation as a vector of `SphericalPoint`s.
Each triangle's lines are separated by a `SphericalPoint(NaN, NaN, NaN)`.
"""
function _to_lines(tri::SphericalTriangulation, n = 20)
    points = get_points(tri)
    nt = num_triangles(tri)
    lines = Vector{eltype(points)}(undef, 3n*(nt+1)) # 3n points per triangle, plus one for the NaN
    offset = 0
    for T in each_triangle(tri)
        ST = SphericalTriangle(tri, T)
        coords = _to_coords(ST, true, n)
        lines[offset+1:offset+length(coords)] .= coords
        lines[offset+length(coords)+1] = SphericalPoint(NaN, NaN, NaN)
        offset += length(coords) + 1
    end
    return lines
end

"""
    SphericalPolygon(lines::Vector{<:SphericalLine})

A polygon on the surface of a sphere defined by the lines `lines`. The polygon is stored as a vector of `SphericalLine`s.
"""
struct SphericalPolygon{P,T} 
    lines::Vector{SphericalLine{P,T}}
end
function SphericalPolygon(points::Vector{<:SphericalPoint})
    n = length(points)
    @assert points[begin] != points[end]
    lines = map(eachindex(points)) do i
        SphericalLine(points[i], points[mod1(i+1, n)])
    end
    return SphericalPolygon(lines)
end

"""
    SphericalPolygon(vorn::SphericalTessellation, v)

Given a `SphericalTessellation` `vorn` and a polygon index `v`, returns the [`SphericalPolygon`](@ref) corresponding to the Voronoi polygon of `v`.
"""
function SphericalPolygon(vorn::SphericalTessellation, v)
    polygon = get_polygon(vorn, v)
    points = [get_polygon_point(vorn, i) for i in view(polygon, 1:(length(polygon)-1))]
    return SphericalPolygon(invproject.(points))
end

"""
    _to_coords(P::SphericalPolygon, add_end=true, n=20) -> Vector{NTuple{3, Number}}

Returns the coordinates of the points on the polygon `P` as a vector of `NTuple{3, Number}`. The polygon is evaluated at `n` equally spaced points on each line.
The first point is repeated at the end of the vector if `add_end` is `true`.
"""
function _to_coords(P::SphericalPolygon, add_end=true, n=20)
    coords = Vector{eltype(P.lines)}()
    for L in P.lines
        _coords = _to_coords(L, n)
        coords = vcat(coords, _coords)
    end
    add_end && push!(coords, coords[1])
    return coords
end

"""
    _to_lines(vorn::SphericalTessellation, n=20) -> Vector{<:SphericalPoint}

Given a `SphericalTessellation` `vorn`, returns the points on the lines of the tessellation as a vector of `SphericalPoint`s.
Each polygon's lines are separated by a `SphericalPoint(NaN, NaN, NaN)`. 
"""
function _to_lines(vorn::SphericalTessellation, n = 20)
    tri = get_triangulation(vorn)
    spoints = get_points(tri)
    polygon_points = get_polygon_points(vorn)
    nt = num_polygons(vorn)
    lines = Vector{eltype(spoints)}(undef, 0)
    sizehint!(lines, (6*3)*n*(nt+1))
    for v in each_polygon_index(vorn)
        spoly = SphericalPolygon(vorn, v)
        coords = _to_coords(spoly, true, n)
        append!(lines, coords)
        push!(lines, SphericalPoint(NaN, NaN, NaN))
    end
    return lines
end