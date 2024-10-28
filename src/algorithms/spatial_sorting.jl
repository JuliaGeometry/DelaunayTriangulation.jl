#=
This implementation is basically a refined version from the implementation
in https://github.com/JuliaGeometry/DistMesh.jl/blob/master/src/hilbertsort.jl,
taken with permission from @sjkelly. Some changes have been made to make the algorithm
easier to understand from the code itself. The license for
the original code from DistMesh.jl is reprinted below (note that DistMesh.jl's repository has a different license
than the one for this code itself at https://github.com/JuliaGeometry/DistMesh.jl/blob/4f0a64a44bd807e2aef7bdfdadcf3a33386619ee/src/hilbertsort.jl#L2).

### DistMesh.jl's hilbertsort.jl license 

This file is licensed under the MIT "Expat" License.

implementing scale-free Hilbert ordering. Real all about it here:
http://doc.cgal.org/latest/Spatial_sorting/index.html

original implementation in:
https://github.com/JuliaGeometry/GeometricalPredicates.jl
The GeometricalPredicates.jl package is licensed under the MIT "Expat" License:
Copyright (c) 2014: Ariel Keselman.

modifications for StaticArrays: sjkelly (Under terms of MIT License)
=#

"""
    separate!(v, p; carry = nothing) -> Integer

Separates the vector `v` in-place based on the predicate `p`, returning the index 
`i` so that `v[1:i]` satisfies `p` and `v[i+1:end]` does not. This is equivalent to 
returning 

    [filter(p, v); filter(!p, v)]
    i = count(p, v)

although the order of the elements is not necessarily preserved.

If `carry` is not nothing, then `carry` is permuted as `v` is. 
"""
function separate!(v, p::F; carry=nothing) where {F} # from Andy Dienes on Slack
    isempty(v) && return firstindex(v)
    i = firstindex(v)
    j = lastindex(v)
    while i < j
        pi = p(v[i])
        pj = p(v[j])
        if !pi || pj
            v[i], v[j] = v[j], v[i]
            if !isnothing(carry)
                carry[i], carry[j] = carry[j], carry[i]
            end
            !pi && (j -= 1)
            pj && (i += 1)
        else
            i += 1
            j -= 1
        end
    end
    return i - !p(v[i])
end

"""
    triseparate!(v, p; carry = nothing) -> (Integer, Integer)

Separates the vector `v` in-place based on the three-valued predicate `p`, returning
two indices `(i, j)` such that:

    v[1:i-1]    satisfies `p(v[i]) < 0`
    v[i:j-1]    satisfies `p(v[i]) == 0`
    v[j:end]    satisfies `p(v[i]) > 0`

This is equivalent to returning:

    [filter(x -> p(x) < 0, v); 
     filter(x -> p(x) == 0, v); 
     filter(x -> p(x) > 0, v)]
     
    i = count(x -> p(x) < 0, v) + 1
    j = count(x -> p(x) == 0, v) + i

However, the order of the elements is not necessarily preserved.

If `carry` is not nothing, then `carry` is permuted in the same way as `v`.
"""
function triseparate!(v, p::F; carry=nothing) where {F}
    isempty(v) && return (firstindex(v), firstindex(v))

    i = firstindex(v)  # Start index for -1 partition
    j = firstindex(v)  # Start index for 0 partition
    k = lastindex(v)   # Start index for 1 partition

    while j ≤ k
        pj = p(v[j])

        if pj < 0
            v[i], v[j] = v[j], v[i]
            if !isnothing(carry)
                carry[i], carry[j] = carry[j], carry[i]
            end
            i += 1
            j += 1
        elseif pj > 0
            v[j], v[k] = v[k], v[j]
            if !isnothing(carry)
                carry[j], carry[k] = carry[k], carry[j]
            end
            k -= 1
        else  # pj == 0
            j += 1
        end
    end
    return i, j
end


"""
    get_median!(v, p; rng=Random.default_rng(), carry=nothing)

Returns the median value in `v`, which gets 
sorted in place, where the median is defined relative to the 
ordering defined by the binary three-valued predicate `p`. The predicate `p` 
gets used as

    p(x, pivot)

where `pivot` is a randomly chosen element from `v` and `x` is 
an element from `v`. The predicate should return `-1`, `0`, or `1`.

The `rng` keyword argument is used for providing a random number 
generator to use in the quickselect algorithm. If `carry`
is not nothing, then `carry` is permuted as `v` is. 

!!! note "Even number of elements"

    Instead of averaging between median elements in the case of an 
    even number of elements, this function returns the lower median.
    In particular, the median returned is always the `(ℓ + 1) ÷ 2` 
    order statistic, where `ℓ` is the length of `v`.
"""
function get_median!(v, p::F; rng=Random.default_rng(), carry=nothing) where {F} # Could also implement http://erdani.org/research/sea2017.pdf
    ℓ = length(v)
    return quickselect!(v, p, (ℓ + 1) ÷ 2; rng=rng, carry=carry)
end

"""
    QuickSelectPredicate{F, T}

Struct to store a predicate for quickselect.

# Fields
- `p::F`: The binary three-valued predicate function, which should return `-1`, `0`, or `1`.
- `pivot::T`: The pivot value to use in the predicate.

Given a value `x` and a `QuickSelectPredicate` `q`, the predicate 
can be invoked as 

    q(x)

which is equivalent to `q.p(x, q.pivot)`.
"""
struct QuickSelectPredicate{F,T}
    p::F
    pivot::T
end
(q::QuickSelectPredicate)(x) = q.p(x, q.pivot)

"""
    quickselect!(v, p, k; rng=Random.default_rng(), carry=nothing) 

Returns the `k`-th order statistic in `v`, which gets
sorted in place, where the order statistic is defined relative to the
ordering defined by the binary predicate `p`. The predicate `p` gets
used as

    p(x, pivot)

where `pivot` is a randomly chosen element from `v` and `x` is an
element from `v`.

The `rng` keyword argument is used for providing a random number
generator to use in the quickselect algorithm. If `carry`
is not nothing, then `carry` is permuted as `v` is.
"""
function quickselect!(v, p::F, k; rng=Random.default_rng(), carry=nothing) where {F}
    length(v) == 1 && return only(v) # assumes k == 1 
    @show v, k
    pivot = rand(rng, v)
    i, j = triseparate!(v, QuickSelectPredicate(p, pivot); carry)
    lows = view(v, firstindex(v):i-1)
    highs = view(v, j:lastindex(v))
    pivots = view(v, i:j-1)
    if k ≤ length(lows)
        return quickselect!(lows, p, k; rng, carry=isnothing(carry) ? carry : view(carry, firstindex(v):i))
    elseif k ≤ length(lows) + length(pivots)
        return first(pivots)
    else
        return quickselect!(highs, p, k - length(lows) - length(pivots); rng, carry=isnothing(carry) ? carry : view(carry, i+1:lastindex(v)))
    end
end

"""
    abstract type AbstractCoordinate

Abstract type for coordinate types.
"""
abstract type AbstractCoordinate end

"""
    struct X <: AbstractCoordinate

Represents the `X` coordinate.
"""
struct X <: AbstractCoordinate end

"""
    struct Y <: AbstractCoordinate

Represents the `Y` coordinate.
"""
struct Y <: AbstractCoordinate end

"""
    next(c::AbstractCoordinate) -> AbstractCoordinate

Returns the next coordinate type, so `next(X()) == Y()` and `next(Y()) == X()`.
"""
next(::X) = Y()
next(::Y) = X()

"""
    get(c::AbstractCoordinate) -> Function

Returns the function to get the coordinate value from a point,
so `get(X())(p) === getx(p)` and `get(Y())(p) === gety(p)`.
"""
get(::X) = getx
get(::Y) = gety

"""
    abstract type AbstractDirection

Abstract type for direction types.
"""
abstract type AbstractDirection end

"""
    struct Left <: AbstractDirection

Represents the `Left` direction, which is equivalent to the `<` operator.
"""
struct Left <: AbstractDirection end

"""
    struct NotLeft <: AbstractDirection

Represents the `NotLeft` direction, which is equivalent to the `≥` operator.
"""
struct NotLeft <: AbstractDirection end # Right or Equal 

"""
    op(d::AbstractDirection) -> Function

Returns the function to compare two values in the direction `d`,
so `op(Left())(a, b)` is equivalent to `a < b` and `op(NotLeft())(a, b)` is equivalent to `a ≥ b`.
"""
op(::Left) = <
op(::NotLeft) = ≥

"""
    reverse(d::AbstractDirection) -> AbstractDirection

Returns the reverse direction of `d`, so `reverse(Left()) == NotLeft()` and `reverse(NotLeft()) == Left()`.
"""
reverse(::Left) = NotLeft()
reverse(::NotLeft) = Left()

"""
    struct HilbertOp{D<:AbstractDirection,C<:AbstractCoordinate,T}

Represents an operation to compare a point in a Hilbert sort. 

# Fields 
- `dir::AbstractDirection`: The direction to compare the point.
- `coord::AbstractCoordinate`: The coordinate being compared.
- `pivot::T`: The value to compare the coordinate to.

To use this struct to make a comparison, use the method 

    (o::HilbertOp)(p)

which calls 

    op(o.dir)(get(o.coord)(p), o.pivot)
"""
struct HilbertOp{D<:AbstractDirection,C<:AbstractCoordinate,T}
    dir::D
    coord::C
    pivot::T
end
@inline (o::HilbertOp)(p) = op(o.dir)(get(o.coord)(p), o.pivot)

"""
    abstract type AbstractSplitter

Abstract type for splitter types.
"""
abstract type AbstractSplitter end

"""
    struct Middle <: AbstractSplitter

Represents the `Middle` splitter, which splits the bounding boxes in the middle
during a Hilbert sort.
"""
struct Middle <: AbstractSplitter end

"""
    struct Median <: AbstractSplitter

Represents the `Median` splitter, which splits the bounding boxes at the median of the points 
in each dimension during a Hilbert sort.

Not yet implemented.
"""
struct Median <: AbstractSplitter end

"""
    HilbertDirections{C<:AbstractCoordinate, XD<:AbstractDirection, YD<:AbstractDirection}

Represents the directions for a Hilbert sort, defining the quadrant order.

# Fields
- `coord::C`: The coordinate to start the quadrant ordering.
- `xdir::XD`: The direction to order the initial quadrant's x-coordinate.
- `ydir::YD`: The direction to order the initial quadrant's y-coordinate.
"""
struct HilbertDirections{C<:AbstractCoordinate,XD<:AbstractDirection,YD<:AbstractDirection}
    coord::C
    xdir::XD
    ydir::YD
end

"""
    HilbertSorter{P, T, I, D<:HilbertDirections, S<:AbstractSplitter, R<:Random.AbstractRNG}

Struct to store the data for a Hilbert sort.

# Fields
- `points::P`: The points to sort.
- `perm::Vector{Int}`: The permutation order of the points.
- `bbox::NTuple{4, T}`: The bounding box of the points.
- `first::I`: The first index of the points to sort.
- `last::I`: The last index of the points to sort.
- `directions::D`: The directions for the Hilbert sort.
- `splitter::S`: The splitter type for the Hilbert sort.
- `capacity::Int`: The capacity of the sorter. This defines the maximum number of points in a single Hilbert square.
- `rng::R`: The random number generator to use for the sort. Only used if `splitter === Median()`.
"""
struct HilbertSorter{P,T,I,D<:HilbertDirections,S<:AbstractSplitter,R<:Random.AbstractRNG}
    points::P
    perm::Vector{Int}
    bbox::NTuple{4,T}
    first::I
    last::I
    directions::D
    splitter::S
    capacity::Int
    rng::R
end

"""
    HilbertSorter(points; capacity=8, splitter=Median(), rng=Random.default_rng(), perm=convert(Vector{Int}, collect(each_point_index(points))))) -> HilbertSorter

Initialises a `HilbertSorter` object for the points `points`.

# Arguments
- `points`: The points to sort.

# Keyword Arguments
- `capacity=8`: The capacity of the sorter. This defines the maximum number of points in a single Hilbert square.
- `splitter=Median()`: The splitter type for the Hilbert sort, defining how to split the bounding boxes.
- `rng=Random.default_rng()`: The random number generator to use for the sort. Only used if `splitter === Median()`.
- `perm=convert(Vector{Int}, collect(each_point_index(points)))`: The vector to mutate to eventually store the permutation order of the points. (The reason for allowing this to be provided is to simplify [`hilbert_sort_rounds`](@ref).)
"""
function HilbertSorter(points; capacity=8, splitter=Median(), rng=Random.default_rng(), perm=convert(Vector{Int}, collect(each_point_index(points))))
    points = PointsWrapper(copy(points))
    T = number_type(points)
    indices = each_point_index(points)
    i, j = first(indices), last(indices)
    bbox = T.(bounds(bounding_box(points)))
    directions = HilbertDirections(X(), Left(), Left())
    return HilbertSorter(points, perm, bbox, i, j, directions, splitter, capacity, rng)
end

"""
    get_split(sorter::HilbertSorter) -> NTuple{2,T}

Returns the split values for the sorter.
"""
@inline get_split(sorter::HilbertSorter) = get_split(sorter, sorter.splitter)

"""
    get_split(sorter::HilbertSorter, ::Middle) -> NTuple{2,T}

Returns the split values for the sorter using the `Middle` splitter,
which splits the bounding box in the middle.
"""
@inline get_split(sorter, ::Middle) = get_middle_split(sorter.bbox)

"""
    get_split(sorter::HilbertSorter, ::Median) -> NTuple{2,T}

Returns the split values for the sorter using the `Median` splitter,
which splits the bounding box at the median of the points in each dimension.
"""
@inline get_split(sorter, ::Median) = get_median_split!(sorter.points, sorter.perm, sorter.first, sorter.last, sorter.directions.coord, sorter.rng)

"""
    get_middle_split(bbox) -> NTuple{2,T}

Returns the split values for the bounding box `bbox` using the `Middle` splitter,
which splits the bounding box in the middle.
"""
@inline function get_middle_split(bbox)
    xmin, xmax, ymin, ymax = bbox
    return midpoint(xmin, xmax), midpoint(ymin, ymax)
end

"""
    struct HilbertMedianComparator{C}

Comparator for the median of a Hilbert sort.

# Fields
- `coord::C`: The coordinate to compare.

To use this struct to compare two points `p` and `q`, use the method

    (o::HilbertMedianComparator)(p, q)
"""
struct HilbertMedianComparator{C}
    coord::C
end
@inline function (o::HilbertMedianComparator)(p, q)
    x, y = get(o.coord)(p), get(o.coord)(q)
    return x < y ? -1 : x > y ? 1 : 0
end

"""
    get_median_split!(points, perm, i, j, coord) -> NTuple{2,T}

Returns the split values for the points `points` in the range `i:j` using the `Median` splitter,
which splits the bounding box at the median of the points in each dimension. The `perm`
array is used to store the permutation order of the points. The `coord` argument defines 
the order of the dimensions, with `xsplit` being the median of the `coord` dimension and
`ysplit` being the median of the next dimension.
"""
@inline function get_median_split!(points, perm, i, j, coord, rng)
    _points = view(points, i:j)
    _perm = view(perm, i:j)
    xsplit = get_median!(_points, HilbertMedianComparator(coord); rng, carry=_perm)
    ysplit = get_median!(_points, HilbertMedianComparator(next(coord)); rng, carry=_perm)
    return get(coord)(xsplit), get(next(coord))(ysplit)
end

"""
    span(sorter::HilbertSorter) -> Int

Returns the number of points in the sorter.
"""
span(sorter::HilbertSorter) = sorter.last - sorter.first + 1

"""
    check_bbox(sorter::HilbertSorter) -> Bool

Checks if the bounding box of the sorter is too small. See 
also [`check_precision`](@ref).
"""
function check_bbox(sorter::HilbertSorter)
    xmin, xmax, ymin, ymax = sorter.bbox
    A = abs(xmax - xmin) * abs(ymax - ymin)
    return check_precision(A)
end

function stop_sorter(sorter::HilbertSorter)
    return (span(sorter) ≤ sorter.capacity) || check_bbox(sorter)
end

"""
    quadrant_sort!(sorter::HilbertSorter, xsplit, ysplit) -> NTuple{5,Int}

Sorts the data in `sorter` into quadrants based on the `xsplit` and `ysplit` values
and in the directions defined by the `sorter.directions`. Returns indices 
`(i₁, i₂, i₃, i₄, i₅)` such that 

- Points in `sorter.points[i₁:i₂-1]` are in the first quadrant.
- Points in `sorter.points[i₂:i₃-1]` are in the second quadrant.
- Points in `sorter.points[i₃:i₄-1]` are in the third quadrant.
- Points in `sorter.points[i₄:i₅]` are in the fourth quadrant.

Here, the quadrants are enumerated in the order defined by `sorter.directions` 
(so what is "first", for example, depends on the `sorter.directions`).
"""
function quadrant_sort!(sorter::HilbertSorter, xsplit, ysplit)
    points = sorter.points
    directions = sorter.directions
    op1 = HilbertOp(directions.xdir, directions.coord, xsplit)
    op2 = HilbertOp(directions.ydir, next(directions.coord), ysplit)
    op3 = HilbertOp(reverse(directions.ydir), next(directions.coord), ysplit)
    i₁, i₅ = sorter.first, sorter.last
    i₃ = separate!(view(points, i₁:i₅), op1; carry=view(sorter.perm, i₁:i₅)) + i₁
    i₂ = separate!(view(points, i₁:i₃-1), op2; carry=view(sorter.perm, i₁:i₃-1)) + i₁
    i₄ = separate!(view(points, i₃:i₅), op3; carry=view(sorter.perm, i₃:i₅)) + i₃
    return i₁, i₂, i₃, i₄, i₅
end

function to_quadrant1(sorter::HilbertSorter, xsplit, ysplit, i₁, i₂, _, _, _)
    xmin, _, ymin, _ = sorter.bbox
    rotated_bbox = (ymin, ysplit, xmin, xsplit) # scaled, rotated, and reflected
    directions = sorter.directions
    new_directions = HilbertDirections(next(directions.coord), directions.ydir, directions.xdir)
    return HilbertSorter(sorter.points, sorter.perm, rotated_bbox, i₁, i₂ - 1, new_directions, sorter.splitter, sorter.capacity, sorter.rng)
end
function to_quadrant2(sorter::HilbertSorter, xsplit, ysplit, _, i₂, i₃, _, _)
    xmin, _, _, ymax = sorter.bbox
    rotated_bbox = (xmin, xsplit, ysplit, ymax) # scaled
    directions = sorter.directions
    new_directions = HilbertDirections(directions.coord, directions.ydir, directions.xdir)
    return HilbertSorter(sorter.points, sorter.perm, rotated_bbox, i₂, i₃ - 1, new_directions, sorter.splitter, sorter.capacity, sorter.rng)
end
function to_quadrant3(sorter::HilbertSorter, xsplit, ysplit, _, _, i₃, i₄, _)
    _, xmax, _, ymax = sorter.bbox
    rotated_bbox = (xsplit, xmax, ysplit, ymax) # scaled
    return HilbertSorter(sorter.points, sorter.perm, rotated_bbox, i₃, i₄ - 1, sorter.directions, sorter.splitter, sorter.capacity, sorter.rng)
end
function to_quadrant4(sorter::HilbertSorter, xsplit, ysplit, _, _, _, i₄, i₅)
    _, xmax, ymin, _ = sorter.bbox
    rotated_bbox = (ysplit, ymin, xmax, xsplit) # scaled, rotated, and reflected
    directions = sorter.directions
    new_directions = HilbertDirections(next(directions.coord), reverse(directions.ydir), reverse(directions.xdir))
    return HilbertSorter(sorter.points, sorter.perm, rotated_bbox, i₄, i₅, new_directions, sorter.splitter, sorter.capacity, sorter.rng)
end

"""
    hilbert_sort(points; capacity = 8, splitter = Middle(), rng=Random.default_rng()) -> Tuple{P, Vector{Int}, Int)

Sorts the points `points` using a Hilbert sort. 

# Arguments 
- `points`: The points to sort.

# Keyword Arguments
- `capacity = 8`: The capacity of the sorter. This defines the maximum number of points in a single Hilbert square.
- `splitter = Middle()`: The splitter type for the Hilbert sort, defining how to split the bounding boxes.
- `rng = Random.default_rng()`: The random number generator to use for the sort. Only used if `splitter === Median()`.

!!! warning "Median"

    Median splitting is not yet implemented.

# Outputs 
- `spoints`: The sorted points. These are not aliased with the input `points`.
- `perm`: The permutation order of the points, so that `[get_point(points, i) for i in perm] == each_point(spoints)`.
- `order`: The order of the Hilbert sort. This defines the order of the Hilbert curve that sorts the points.
"""
function hilbert_sort(points; capacity=8, splitter=Middle(), rng=Random.default_rng())
    if splitter === Median()
        throw(ArgumentError("Median splitting is not yet implemented."))
        #=
        The reason it's not yet implemented is because 
            E = (12.8728090546081, 1.057828011921)
            F = (2.9914997728536, 4.9162440171775)
            G = (4.9677616292045, 10.8920834399529)
            H = (2.7091766505178, 21.1498235514886)
            I = (11.0847626131478, 18.9853462802471)
            J = (17.0263810513976, 26.8690055903223)
            K = (23.2289344966548, 24.9012989801027)
            L = (22.9722771127131, 22.7197112165985)
            M = (25.0683124149035, 22.9763686005402)
            N = (25.1110886455604, 16.6027102326552)
            O = (24.9399837229326, 15.1055421596621)
            P = (31.014208476219, 9.3307510209744)
            Q = (25.0683124149035, 6.8497296428715)
            R = (28.6938091798674, 0.903833581556)
            S = (22.9722771127131, 0.7755048895852)
            T = (30.6719986309634, 31.1038524253599)
            points = [E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T]
            spoints, perm, order = DelaunayTriangulation.hilbert_sort(points; capacity=1, splitter=DelaunayTriangulation.Median())
        leads to a stack overflow for some reason. Maybe I am not understanding the actual definition of median splitting.
        =#
    end
    sorter = HilbertSorter(points; capacity, splitter, rng)
    order = hilbert_sort!(sorter)
    return unwrap(sorter.points), sorter.perm, order
end

"""
    hilbert_sort!(sorter::HilbertSorter[, depth = 0]) -> Int 

Sorts `sorter` using a Hilbert sort, counting the depth recursively. The returned depth 
is the maximum depth reached.
"""
function hilbert_sort!(sorter::HilbertSorter, depth=0)
    stop_sorter(sorter) && return depth
    depth += 1
    xsplit, ysplit = get_split(sorter)
    i₁, i₂, i₃, i₄, i₅ = quadrant_sort!(sorter, xsplit, ysplit)
    depth1 = hilbert_sort!(to_quadrant1(sorter, xsplit, ysplit, i₁, i₂, i₃, i₄, i₅), depth)
    depth2 = hilbert_sort!(to_quadrant2(sorter, xsplit, ysplit, i₁, i₂, i₃, i₄, i₅), depth)
    depth3 = hilbert_sort!(to_quadrant3(sorter, xsplit, ysplit, i₁, i₂, i₃, i₄, i₅), depth)
    depth4 = hilbert_sort!(to_quadrant4(sorter, xsplit, ysplit, i₁, i₂, i₃, i₄, i₅), depth)
    return max(depth1, depth2, depth3, depth4)
end

"""
    assign_rounds(indices; rng=Random.default_rng()) -> Vector{Vector{Int}}

Given `indices`, assign rounds to traverse the indices. In particular,
there are `ceil(Int, log2(n)) + 1` rounds, where `n = length(indices)`, and 

- Indices are assigned to the final round with probability 1/2.
- To the next round with probability 1/2,
- and so on.

The returned value is a `Vector` that maps a round number to the indices 
participating in that round. You can use the `rng` keyword to control the 
random number generator.
"""
function assign_rounds(indices; rng=Random.default_rng())
    n = length(indices)
    nr = ceil(Int, log2(n)) + 1
    participant_map = Vector{Vector{Int}}(undef, nr)
    entered = falses(n)
    for i in 1:(nr-1)
        participants = Int[]
        sizehint!(participants, ceil(Int, (1 / 2)^i * n))
        for j in indices
            entered[j] && continue
            if rand(rng) < 1 / 2
                push!(participants, j)
                entered[j] = true
            end
        end
        participant_map[i] = participants
    end
    # Put all the remaining participants into the last round
    participant_map[end] = filter(i -> !entered[i], indices)
    # Reverse the rounds so that those in the first round are first
    reverse!(participant_map)
    return participant_map
end

"""
    hilbert_sort_rounds(points; capacity = 2, splitter = Middle(), rng=Random.default_rng())

Sorts the points `points` using a Hilbert sort with rounds. In particular, the points are placed 
into rounds using [`assign_rounds`](@ref) and then points are sorted using a Hilbert sort
within-rounds. 

# Arguments
- `points`: The points to sort.

# Keyword Arguments
- `capacity = 2`: The capacity of the sorter. This defines the maximum number of points in a single Hilbert square.
- `splitter = Middle()`: The splitter type for the Hilbert sort, defining how to split the bounding boxes.
- `rng = Random.default_rng()`: The random number generator to use for the sort and for assigning rounds.

# Outputs
- `round_perms`: These are the permutation vectors associated with each round. The ordering of the permutation vectors defines the order of the points along the Hilbert curve. 
"""
function hilbert_sort_rounds(points; capacity=2, splitter=Middle(), rng=Random.default_rng())
    if splitter === Median()
        throw(ArgumentError("Median splitting is not yet implemented."))
    end
    indices = each_point_index(points)
    rounds = assign_rounds(indices; rng)
    pw = PointsWrapper(copy(points))
    for r in rounds 
        isempty(r) && continue 
        # We don't use a view here since `r` is a vector of indices,
        # so a view is a bit slower than just copying the points.
        sorter = HilbertSorter(pw[r]; capacity, splitter, rng, perm=r) 
        hilbert_sort!(sorter)
    end
    return rounds
end
