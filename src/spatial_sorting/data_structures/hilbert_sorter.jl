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

function HilbertSorter(points; capacity=8, splitter=Median(), rng=Random.default_rng(), perm=convert(Vector{Int}, collect(each_point_index(points))))
    points = PointsWrapper(copy(points))
    T = number_type(points)
    indices = each_point_index(points)
    i, j = first(indices), last(indices)
    bbox = T.(bounds(bounding_box(points)))
    directions = HilbertDirections(X(), Left(), Left())
    return HilbertSorter(points, perm, bbox, i, j, directions, splitter, capacity, rng)
end

@inline get_split(sorter::HilbertSorter) = get_split(sorter, sorter.splitter)

span(sorter::HilbertSorter) = sorter.last - sorter.first + 1

function check_bbox(sorter::HilbertSorter)
    xmin, xmax, ymin, ymax = sorter.bbox
    A = abs(xmax - xmin) * abs(ymax - ymin)
    return check_precision(A)
end

function stop_sorter(sorter::HilbertSorter)
    return (span(sorter) ≤ sorter.capacity) || check_bbox(sorter)
end