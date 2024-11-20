# Interface
@inline getx(p) = p[1]

@inline gety(p) = p[2]

@inline getz(p) = p[3]

@inline getpoint(points, i::Integer) = getxy(points[i])
@inline getpoint(points::AbstractMatrix, i::Integer) = (points[1, i], points[2, i])
@inline getpoint(_, p) = p # So that we can mix points and vertices

@inline get_point(points, i) = getpoint(points, i)

@inline each_point_index(points) = eachindex(points)
@inline each_point_index(points::AbstractMatrix) = axes(points, 2)

@inline each_point(points) = points
@inline each_point(points::AbstractMatrix) = eachcol(points)

@inline num_points(points) = length(points)
@inline num_points(points::AbstractMatrix) = size(points, 2)

@inline push_point!(points, p) = push_point!(points, getx(p), gety(p))

@inline pop_point!(points) = pop!(points)
@inline pop_point!(points::AbstractMatrix) = resize!(points, (2, num_points(points) - 1))

@inline set_point!(points, i, p) = set_point!(points, i, getx(p), gety(p))

@inline find_point_index(points, p) = find_point_index(points, getx(p), gety(p))

@inline is_planar(points) = all(is_point2, each_point(points))
@inline is_planar(points::AbstractMatrix) = size(points, 1) == 2
@inline is_planar(points::AbstractVector{T}) where {T} = T <: NTuple{2,<:Number}

# Derived 
@inline getxy(p) = (getx(p), gety(p))
@inline getxyz(p) = (getx(p), gety(p), getz(p))

@inline _getx(p) = Float64(getx(p))
@inline _gety(p) = Float64(gety(p))
@inline _getz(p) = Float64(getz(p))

@inline _getxy(p) = (_getx(p), _gety(p))
@inline _getxyz(p) = (_getx(p), _gety(p), _getz(p))

@inline is_point2(p) = eltype(p) <: Number && length(p) == 2
@inline is_point3(p) = eltype(p) <: Number && length(p) == 3

@inline function get_point(points, i, j::Vararg{Any,N}) where {N}
    pᵢ = getpoint(points, i)
    tup = ntuple(k -> (@inline; get_point(points, j[k])), Val(N))
    return (pᵢ, tup...)
end

function points_are_unique(points)
    n = num_points(points)
    T = number_type(points)
    seen = Set{NTuple{2,T}}()
    sizehint!(seen, n)
    for i in each_point_index(points)
        p = get_point(points, i)
        p ∈ seen && return false
        push!(seen, p)
    end
    return n == length(seen)
end

function lexicographic_order(points)
    # equivalent to sortperm(collect(each_point(points)))
    itr = each_point_index(points)
    lt(i, j) = get_point(points, i) < get_point(points, j)
    return sortperm(itr; lt) # Tuples are compared lexicographically
end

@inline push_point!(points::AbstractVector{T}, x, y) where {F,T<:NTuple{2,F}} = push!(points, (F(x), F(y)))
@inline push_point!(points::AbstractVector{T}, x, y) where {F,T<:Vector{F}} = push!(points, [F(x), F(y)])
@inline push_point!(points::AbstractVector{T}, x, y) where {F,T<:AbstractVector{F}} = push!(points, T((F(x), F(y)))) # StaticArrays
@inline push_point!(points::AbstractMatrix{T}, x, y) where {T<:Number} = append!(points, (x, y)) # ElasticArrays 

@inline set_point!(points::AbstractVector, i, x, y) = (points[i] = (x, y); (x, y))
@inline set_point!(points::AbstractVector{<:Vector}, i, x, y) = (points[i] = [x, y]; (x, y))
@inline set_point!(points::AbstractMatrix, i, x, y) = (points[1, i] = x; points[2, i] = y; (x, y))

function mean_points(points, vertices=each_point_index(points))
    F = number_type(points)
    m, n = (zero(F), zero(F)), zero(Int)
    for v in vertices
        p = get_point(points, v)
        m = m .+ p
        n += 1
    end
    return m ./ n
end

function find_point_index(points, x, y)
    itr = each_point_index(points)
    for i in itr
        p = get_point(points, i)
        p == (x, y) && return i
    end
    return eltype(itr)(∅)
end

function find_duplicate_points(points)
    T = number_type(points)
    n = num_points(points)
    dup_seen = Dict{NTuple{2,T},Vector{Int}}()
    sizehint!(dup_seen, n)
    for i in each_point(points)
        p = get_point(points, i)
        existing_inds = get!(Vector{Int}, dup_seen, p)
        push!(existing_inds, i)
    end
    filter!(dup_seen) do (_, ivec)
        length(ivec) > 1
    end
    return dup_seen
end

@inline function swap_points!(points, i, j)
    pᵢ, pⱼ = get_point(points, i, j)
    set_point!(points, i, pⱼ)
    set_point!(points, j, pᵢ)
    return points
end

struct PointsWrapper{P,T} <: AbstractVector{NTuple{2,T}}
    points::P
    @inline PointsWrapper(points::P) where {P} = new{P,number_type(P)}(points)
end
@inline get_points(pw::PointsWrapper) = pw.points
@inline num_points(pw::PointsWrapper) = num_points(unwrap(pw))
@inline get_point(pw::PointsWrapper, i) = get_point(unwrap(pw), i)
@inline set_point!(pw::PointsWrapper, i, p) = set_point!(unwrap(pw), i, p)
@inline each_point_index(pw::PointsWrapper) = each_point_index(unwrap(pw))
@inline Base.size(pw::PointsWrapper) = (num_points(pw),)
@inline Base.getindex(pw::PointsWrapper, i::Int) = getxy(get_point(pw, i))
@inline Base.IndexStyle(::Type{<:PointsWrapper}) = Base.IndexLinear()
@inline Base.setindex!(pw::PointsWrapper, p, i::Int) = set_point!(pw, i, p)
@inline Base.axes(pw::PointsWrapper) = (each_point_index(pw),)

@inline unwrap(pw::PointsWrapper) = get_points(pw)
@inline unwrap(pw) = pw
