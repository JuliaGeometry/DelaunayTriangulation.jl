const HILBERT_ORDER = 32 # The grid used for the Hilbert curve has 2^HILBERT_ORDER cells along each axis

"""
    hilbert_index(x, y, order) -> UInt64

Returns the position of the cell `(x, y)` along a Hilbert curve through a `2^order × 2^order` grid of cells, where `x` and `y`
are integers in `0:(2^order - 1)` and `order ≤ 32`. The curve starts at the cell `(0, 0)` and ends at the cell `(2^order - 1, 0)`,
and consecutive positions along the curve are adjacent cells.
"""
function hilbert_index(x, y, order)
    n = one(UInt64) << order
    x, y = UInt64(x), UInt64(y)
    d = zero(UInt64)
    s = n >> 1
    while s > 0
        rx = (x & s) > 0
        ry = (y & s) > 0
        d += s * s * UInt64((3 * rx) ⊻ ry)
        # Rotate the quadrant so that the curve inside it has the same orientation as the curve through the whole grid
        if !ry
            if rx
                x = n - 1 - x
                y = n - 1 - y
            end
            x, y = y, x
        end
        s >>= 1
    end
    return d
end

"""
    hilbert_cell(x, xmin, scale) -> UInt64

Returns the index of the cell along an axis of the grid used for the Hilbert curve in [`brio_order`](@ref) that the coordinate `x` lies in,
where `xmin` is the smallest coordinate along that axis and `scale` is the number of cells per unit length. Coordinates
that are not finite are put into the first cell.
"""
function hilbert_cell(x, xmin, scale)
    t = (Float64(x) - xmin) * scale
    return isfinite(t) ? unsafe_trunc(UInt64, clamp(t, 0.0, Float64(typemax(UInt32)))) : zero(UInt64)
end

"""
    hilbert_bounds(points, vertices) -> (Float64, Float64, Float64, Float64)

Returns the bounding box `(xmin, xmax, ymin, ymax)` of the points of `vertices`, ignoring coordinates that are not finite. If there are no
finite coordinates along an axis, the bounds along that axis are both zero.
"""
function hilbert_bounds(points, vertices)
    xmin, xmax, ymin, ymax = Inf, -Inf, Inf, -Inf
    for i in vertices
        x, y = getxy(get_point(points, i))
        x, y = Float64(x), Float64(y)
        if isfinite(x)
            xmin, xmax = min(xmin, x), max(xmax, x)
        end
        if isfinite(y)
            ymin, ymax = min(ymin, y), max(ymax, y)
        end
    end
    xmin > xmax && ((xmin, xmax) = (0.0, 0.0))
    ymin > ymax && ((ymin, ymax) = (0.0, 0.0))
    return xmin, xmax, ymin, ymax
end

"""
    brio_order(points, vertices, rng::Random.AbstractRNG=Random.default_rng()) -> Vector

Returns the `vertices` in a biased randomised insertion order (BRIO), as described by [Amenta, Choi, and Rote (2003)](https://doi.org/10.1145/777792.777824),
for inserting their `points` into a triangulation.

# Arguments
- `points`: The points.
- `vertices`: The vertices to order.
- `rng::Random.AbstractRNG=Random.default_rng()`: The random number generator to use.

# Output
- `order`: The vertices in the order to insert them in.

# Extended help
Each vertex is put into one of `⌈log₂ n⌉ + 1` rounds, where `n = length(vertices)`: the last round with probability `1/2`, otherwise the round before
with probability `1/2`, and so on, with the first round taking the vertices that remain. The rounds are then inserted one after the other, and the
vertices of each round are sorted along a Hilbert curve (see [`hilbert_index`](@ref)) through the bounding box of the points. The direction of the
curve alternates between rounds, so that each round starts close to where the previous round ended.

Consecutive vertices are then nearly always close to each other, so that the point location for each new point can start from the last
inserted point, and the data that each insertion uses is likely to have been used recently. The rounds keep enough of the randomisation
that the expected running time of the Bowyer-Watson algorithm relies on.
"""
function brio_order(points, vertices, rng::Random.AbstractRNG = Random.default_rng())
    I = eltype(vertices)
    n = length(vertices)
    n == 0 && return I[]
    num_rounds = ceil(Int, log2(n)) + 1
    xmin, xmax, ymin, ymax = hilbert_bounds(points, vertices)
    num_cells = Float64(typemax(UInt32)) # 2^HILBERT_ORDER - 1
    xscale = xmax > xmin ? num_cells / (xmax - xmin) : 0.0
    yscale = ymax > ymin ? num_cells / (ymax - ymin) : 0.0
    entries = Vector{Tuple{Int, UInt64, I}}(undef, n)
    for (k, i) in enumerate(vertices)
        # The number of trailing zeros of a random integer is 0 with probability 1/2, 1 with probability 1/4, and so on
        round = max(1, num_rounds - trailing_zeros(rand(rng, UInt64)))
        x, y = getxy(get_point(points, i))
        key = hilbert_index(hilbert_cell(x, xmin, xscale), hilbert_cell(y, ymin, yscale), HILBERT_ORDER)
        entries[k] = (round, iseven(num_rounds - round) ? key : ~key, i)
    end
    sort!(entries)
    order = Vector{I}(undef, n)
    for (k, (_, _, i)) in enumerate(entries)
        order[k] = i
    end
    return order
end

"""
    spatial_num_samples(n) -> Integer

Returns `1`. This is the default number of samples used in the point location step of [`triangulate`](@ref) when the points are inserted
in the order from [`brio_order`](@ref): the last inserted point is then nearly always close to the next point, so that sampling more points
to start the point location from costs more than it saves.
"""
spatial_num_samples(::I) where {I} = one(I)
