abstract type AbstractSplitter end

struct Middle <: AbstractSplitter end

@inline get_split(sorter, ::Middle) = get_middle_split(sorter.bbox)
@inline function get_middle_split(bbox)
    xmin, xmax, ymin, ymax = bbox
    return midpoint(xmin, xmax), midpoint(ymin, ymax)
end

struct Median <: AbstractSplitter end

@inline get_split(sorter, ::Median) = get_median_split!(sorter.points, sorter.perm, sorter.first, sorter.last, sorter.directions.coord, sorter.rng)

struct HilbertMedianComparator{C}
    coord::C
end
@inline function (o::HilbertMedianComparator)(p, q)
    x, y = get(o.coord)(p), get(o.coord)(q)
    return x < y ? -1 : x > y ? 1 : 0
end

@inline function get_median_split!(points, perm, i, j, coord, rng)
    _points = view(points, i:j)
    _perm = view(perm, i:j)
    xsplit = get_median!(_points, HilbertMedianComparator(coord); rng, carry=_perm)
    ysplit = get_median!(_points, HilbertMedianComparator(next(coord)); rng, carry=_perm)
    return get(coord)(xsplit), get(next(coord))(ysplit)
end