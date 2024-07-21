"""
    getx(p) -> Number

Get the x-coordinate of `p`.

# Examples
```jldoctest
julia> using DelaunayTriangulation

julia> p = (0.3, 0.7);

julia> getx(p)
0.3
```
"""
getx(p) = p[1]

"""
    gety(p) -> Number

Get the y-coordinate of `p`.

# Examples
```jldoctest
julia> using DelaunayTriangulation

julia> p = (0.9, 1.3);

julia> gety(p)
1.3
```
"""
gety(p) = p[2]

"""
    getz(p) -> Number 

Get the z-coordinate of `p`.
"""
getz(p) = p[3]

"""
    getxy(p) -> NTuple{2, Number}

Get the coordinates of `p` as a `Tuple`.

# Examples
```jldoctest
julia> using DelaunayTriangulation

julia> p = [0.9, 23.8];

julia> getxy(p)
(0.9, 23.8)
```
"""
getxy(p) = (getx(p), gety(p))

"""
    getxyz(p) -> NTuple{3, Number}

Given a three-dimensional `p`, returns its coordinates as a `Tuple`.
"""
getxyz(p) = (getx(p), gety(p), getz(p))

"""
    _getx(p) -> Float64

Get the x-coordinate of `p` as a `Float64`.

# Examples
```jldoctest
julia> using DelaunayTriangulation

julia> p = (0.37, 0.7);

julia> DelaunayTriangulation._getx(p)
0.37

julia> p = (0.37f0, 0.7f0);

julia> DelaunayTriangulation._getx(p)
0.3700000047683716
```
"""
_getx(p) = Float64(getx(p))

"""
    _gety(p) -> Float64

Get the y-coordinate of `p` as a `Float64`.

# Examples
```jldoctest
julia> using DelaunayTriangulation

julia> p = (0.5, 0.5);

julia> DelaunayTriangulation._gety(p)
0.5

julia> p = (0.5f0, 0.5f0);

julia> DelaunayTriangulation._gety(p)
0.5
```
"""
_gety(p) = Float64(gety(p))

"""
    _getz(p) -> Float64 

Get the z-coordinate of `p` as a `Float64.`
"""
_getz(p) = Float64(getz(p))

"""
    _getxy(p) -> NTuple{2, Float64}

Get the coordinates of `p` as a `Tuple` of `Float64`s.

# Examples
```jldoctest
julia> using DelaunayTriangulation

julia> p = [0.3, 0.5];

julia> DelaunayTriangulation._getxy(p)
(0.3, 0.5)

julia> p = [0.3f0, 0.5f0];

julia> DelaunayTriangulation._getxy(p)
(0.30000001192092896, 0.5)
```
"""
_getxy(p) = (_getx(p), _gety(p))

"""
    _getxyz(p) -> NTuple{3, Float64}

Given a three-dimemsional `p`, returns its coordinates as a `Tuple` of `Float64`s.
"""
_getxyz(p) = (_getx(p), _gety(p), _getz(p))

@doc """
    getpoint(points, vertex) -> NTuple{2, Number}

Get the point associated with `vertex` in `points`, returned as a `Tuple` of the coordinates. 
If `vertex` is not an integer, then `vertex` is returned so that points and vertices can be easily mixed.

# Examples
```jldoctest
julia> using DelaunayTriangulation

julia> points = [(0.3, 0.7), (1.3, 5.0), (5.0, 17.0)];

julia> DelaunayTriangulation.getpoint(points, 2)
(1.3, 5.0)

julia> points = [0.3 1.3 5.0; 0.7 5.0 17.0];

julia> DelaunayTriangulation.getpoint(points, 2)
(1.3, 5.0)

julia> DelaunayTriangulation.getpoint(points, (17.3, 33.0))
(17.3, 33.0)
```
"""
getpoint
getpoint(points, i::Integer) = getxy(points[i])
getpoint(points::AbstractMatrix, i::Integer) = (points[1, i], points[2, i])
getpoint(points, p) = p # so that we can mix points and vertices

@doc """
    get_point(points, vertices...) -> NTuple{length(vertices), NTuple{2, Number}}

Get the points associated with `vertices` in `points`.

# Examples
```jldoctest 
julia> using DelaunayTriangulation

julia> points = [(1.0, 2.0), (3.0, 5.5), (1.7, 10.3), (-5.0, 0.0)];

julia> get_point(points, 1)
(1.0, 2.0)

julia> get_point(points, 1, 2, 3, 4)
((1.0, 2.0), (3.0, 5.5), (1.7, 10.3), (-5.0, 0.0))

julia> points = [1.0 3.0 1.7 -5.0; 2.0 5.5 10.3 0.0];

julia> get_point(points, 1)
(1.0, 2.0)

julia> get_point(points, 1, 2, 3, 4)
((1.0, 2.0), (3.0, 5.5), (1.7, 10.3), (-5.0, 0.0))

julia> typeof(ans)
NTuple{4, Tuple{Float64, Float64}}
```
"""
get_point
get_point(points, i) = getpoint(points, i)
get_point(points, i::Vararg{Any, N}) where {N} = ntuple(j -> get_point(points, i[j]), Val(N))

@doc """
    each_point_index(points) -> Iterator

Returns an iterator over each point index in `points`.

# Examples
```jldoctest 
julia> using DelaunayTriangulation

julia> points = [(1.0, 2.0), (-5.0, 2.0), (2.3, 2.3)];

julia> DelaunayTriangulation.each_point_index(points)
Base.OneTo(3)

julia> points = [1.0 -5.0 2.3; 2.0 2.0 2.3];

julia> DelaunayTriangulation.each_point_index(points)
Base.OneTo(3)
```
"""
each_point_index
each_point_index(points) = eachindex(points)
each_point_index(points::AbstractMatrix) = axes(points, 2)

@doc """
    each_point(points) -> Iterator

Returns an iterator over each point in `points`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> points = [(1.0, 2.0), (5.0, 13.0)];

julia> DelaunayTriangulation.each_point(points)
2-element Vector{Tuple{Float64, Float64}}:
 (1.0, 2.0)
 (5.0, 13.0)

julia> points = [1.0 5.0 17.7; 5.5 17.7 0.0];

julia> DelaunayTriangulation.each_point(points)
3-element ColumnSlices{Matrix{Float64}, Tuple{Base.OneTo{Int64}}, SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}:
 [1.0, 5.5]
 [5.0, 17.7]
 [17.7, 0.0]
```
"""
each_point
each_point(points) = points
each_point(points::AbstractMatrix) = eachcol(points)

@doc """
    num_points(points) -> Integer

Returns the number of points in `points`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> points = [(1.0, 1.0), (2.3, 1.5), (0.0, -5.0)];

julia> DelaunayTriangulation.num_points(points)
3

julia> points = [1.0 5.5 10.0 -5.0; 5.0 2.0 0.0 0.0];

julia> DelaunayTriangulation.num_points(points)
4
```
"""
num_points
num_points(points) = length(points)
num_points(points::AbstractMatrix) = size(points, 2)

"""
    points_are_unique(points) -> Bool

Returns `true` if all points in `points` are unique.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> points = [1.0 2.0 3.0 4.0 5.0; 0.0 5.5 2.0 1.3 17.0]
2×5 Matrix{Float64}:
 1.0  2.0  3.0  4.0   5.0
 0.0  5.5  2.0  1.3  17.0

julia> DelaunayTriangulation.points_are_unique(points)
true

julia> points[:, 4] .= points[:, 1];

julia> DelaunayTriangulation.points_are_unique(points)
false
```
"""
function points_are_unique(points)
    n = num_points(points)
    T = number_type(points)
    seen = Set{NTuple{2, T}}()
    sizehint!(seen, n)
    for i in each_point_index(points)
        p = get_point(points, i)
        p ∈ seen && return false
        push!(seen, p)
    end
    return n == length(seen)
end

"""
    lexicographic_order(points) -> Vector{Int}

Returns a set of indices that give the lexicographic ordering of `points`,
meaning the indices so that the points are sorted first by `x` and then by `y`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> points = [(1.0, 5.0), (0.0, 17.0), (0.0, 13.0), (5.0, 17.3), (3.0, 1.0), (5.0, -2.0)]
6-element Vector{Tuple{Float64, Float64}}:
 (1.0, 5.0)
 (0.0, 17.0)
 (0.0, 13.0)
 (5.0, 17.3)
 (3.0, 1.0)
 (5.0, -2.0)

julia> DelaunayTriangulation.lexicographic_order(points)
6-element Vector{Int64}:
 3
 2
 1
 5
 6
 4

julia> hcat(points, points[ans])
6×2 Matrix{Tuple{Float64, Float64}}:
 (1.0, 5.0)   (0.0, 13.0)
 (0.0, 17.0)  (0.0, 17.0)
 (0.0, 13.0)  (1.0, 5.0)
 (5.0, 17.3)  (3.0, 1.0)
 (3.0, 1.0)   (5.0, -2.0)
 (5.0, -2.0)  (5.0, 17.3)
```
"""
lexicographic_order(points) = sortperm(collect(each_point(points)))

@doc """
    push_point!(points, x, y)
    push_point!(points, p) = push_point!(points, getx(p), gety(p))

Pushes the point `p = (x, y)` into `points`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> points = [(1.0, 3.0), (5.0, 1.0)]
2-element Vector{Tuple{Float64, Float64}}:
 (1.0, 3.0)
 (5.0, 1.0)

julia> DelaunayTriangulation.push_point!(points, 2.3, 5.3)
3-element Vector{Tuple{Float64, Float64}}:
 (1.0, 3.0)
 (5.0, 1.0)
 (2.3, 5.3)

julia> DelaunayTriangulation.push_point!(points, (17.3, 5.0))
4-element Vector{Tuple{Float64, Float64}}:
 (1.0, 3.0)
 (5.0, 1.0)
 (2.3, 5.3)
 (17.3, 5.0)
```
"""
push_point!
push_point!(points::AbstractVector{T}, x, y) where {F, T <: NTuple{2, F}} = push!(points, (F(x), F(y)))
push_point!(points::AbstractVector{T}, x, y) where {F <: Number, T <: AbstractVector{F}} = push!(points, F[x, y])
push_point!(points::AbstractMatrix{T}, x, y) where {T <: Number} = append!(points, (x, y)) # ElasticArrays 
push_point!(points, p) = push_point!(points, getx(p), gety(p))

@doc """
    pop_point!(points)

Pops the last point from `points`.

# Examples
```jldoctest
julia> using DelaunayTriangulation

julia> points = [(1.0, 2.0), (1.3, 5.3)]
2-element Vector{Tuple{Float64, Float64}}:
 (1.0, 2.0)
 (1.3, 5.3)

julia> DelaunayTriangulation.pop_point!(points) # returns the popped vector
(1.3, 5.3)

julia> points
1-element Vector{Tuple{Float64, Float64}}:
 (1.0, 2.0)
```
"""
pop_point!
pop_point!(points) = pop!(points)
pop_point!(points::AbstractMatrix) = resize!(points, (2, size(points, 2) - 1))

"""
    mean_points(points[, vertices = each_point_index(points)]) -> NTuple{2, Number}

Returns the mean of the points in `points` indexed by `vertices`, 
given as a `Tuple` of the form `(mean_x, mean_y)`.

# Examples
```jldoctest
julia> using DelaunayTriangulation

julia> points = [(1.0, 2.0), (2.3, 5.0), (17.3, 5.3)]
3-element Vector{Tuple{Float64, Float64}}:
 (1.0, 2.0)
 (2.3, 5.0)
 (17.3, 5.3)

julia> DelaunayTriangulation.mean_points(points)
(6.866666666666667, 4.1000000000000005)

julia> (1.0 + 2.3 + 17.3)/3, (2.0 + 5.0 + 5.3)/3
(6.866666666666667, 4.1000000000000005)

julia> points = [1.0 2.3 17.3; 2.0 5.0 5.3]
2×3 Matrix{Float64}:
 1.0  2.3  17.3
 2.0  5.0   5.3

julia> DelaunayTriangulation.mean_points(points)
(6.866666666666667, 4.1000000000000005)

julia> DelaunayTriangulation.mean_points(points, (1, 3))
(9.15, 3.65)

julia> (1.0 + 17.3)/2, (2.0 + 5.3)/2
(9.15, 3.65)
```
"""
function mean_points(points, vertices = each_point_index(points))
    F = number_type(points)
    x = zero(F)
    y = zero(F)
    n = zero(F)
    for v in vertices
        p = get_point(points, v)
        px, py = getxy(p)
        x += px
        y += py
        n += 1
    end
    return (x / n, y / n)
end

@doc """
    set_point!(points, i, x, y)
    set_point!(points, i, p) = set_point!(points, i, getx(p), gety(p))

Sets the point at index `i` in `points` to `(x, y)`.

# Examples
```jldoctest
julia> using DelaunayTriangulation

julia> points = [(1.0, 3.0), (5.0, 17.0)]
2-element Vector{Tuple{Float64, Float64}}:
 (1.0, 3.0)
 (5.0, 17.0)

julia> DelaunayTriangulation.set_point!(points, 1, 0.0, 0.0)
(0.0, 0.0)

julia> points
2-element Vector{Tuple{Float64, Float64}}:
 (0.0, 0.0)
 (5.0, 17.0)

julia> points = [1.0 2.0 3.0; 4.0 5.0 6.0]
2×3 Matrix{Float64}:
 1.0  2.0  3.0
 4.0  5.0  6.0

julia> DelaunayTriangulation.set_point!(points, 2, (17.3, 0.0))
2-element view(::Matrix{Float64}, :, 2) with eltype Float64:
 17.3
  0.0

julia> points
2×3 Matrix{Float64}:
 1.0  17.3  3.0
 4.0   0.0  6.0
```
"""
set_point!
set_point!(points::AbstractVector, i, x, y) = points[i] = (x, y)
set_point!(points::AbstractVector{<:Vector}, i, x, y) = points[i] = [x, y]
set_point!(points::AbstractMatrix, i, x, y) = points[:, i] .= (x, y)
set_point!(points, i, p) = set_point!(points, i, getx(p), gety(p))

"""
    find_point_index(points, x, y) -> Integer
    find_point_index(points, p) -> Integer

Returns an index of a point in `points` that is equal to `p = (x, y)`.
If no such point exists, then `$∅` is returned.
"""
function find_point_index(points, x, y)
    for i in each_point_index(points)
        p = get_point(points, i)
        getxy(p) == (x, y) && return i
    end
    return ∅
end
find_point_index(points, p) = find_point_index(points, getx(p), gety(p))

"""
    find_duplicate_points(points) -> Dict{Point, Vector{Int}}

Returns a `Dict` of `points` that maps each duplicate point to a `Vector` of the indices of the duplicate points.
"""
function find_duplicate_points(points)
    T = number_type(points)
    n = num_points(points)
    dup_seen = Dict{NTuple{2, T}, Vector{Int}}()
    sizehint!(dup_seen, n)
    for i in each_point_index(points)
        p = get_point(points, i)
        existing_inds = get!(Vector{Int}, dup_seen, p)
        push!(existing_inds, Int(i))
    end
    filter!(dup_seen) do (_, ivec)
        length(ivec) > 1
    end
    return dup_seen
end
