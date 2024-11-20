abstract type AbstractBoundingShape end

struct BoundingInterval <: AbstractBoundingShape
    a::Float64
    b::Float64
end

bounds(I::BoundingInterval) = (I.a, I.b)

const InvalidBoundingInterval = BoundingInterval(NaN, NaN)

width(I::BoundingInterval) = I.b - I.a

Base.isempty(I::BoundingInterval) = isnan(I.a) || isnan(I.b) || width(I) < 0

midpoint(I::BoundingInterval) = midpoint(I.a, I.b)

function Base.:(∩)(I::BoundingInterval, J::BoundingInterval)
    a′ = max(I.a, J.a)
    b′ = min(I.b, J.b)
    if a′ ≤ b′
        return BoundingInterval(a′, b′)
    else
        return InvalidBoundingInterval
    end
end

function Base.:(∪)(I::BoundingInterval, J::BoundingInterval)
    a′ = min(I.a, J.a)
    b′ = max(I.b, J.b)
    return BoundingInterval(a′, b′)
end

Base.in(a::Float64, I::BoundingInterval) = I.a ≤ a ≤ I.b

Base.in(I::BoundingInterval, J::BoundingInterval) = I.a ∈ J && I.b ∈ J

struct BoundingBox <: AbstractBoundingShape
    x::BoundingInterval
    y::BoundingInterval
end

bounds(r::BoundingBox) = (bounds(r.x)..., bounds(r.y)...)

const InvalidBoundingBox = BoundingBox(InvalidBoundingInterval, InvalidBoundingInterval)

BoundingBox(a, b, c, d) = BoundingBox(BoundingInterval(a, b), BoundingInterval(c, d))

BoundingBox(p::NTuple{2, <:Number}) = BoundingBox(getx(p), getx(p), gety(p), gety(p))

width(r::BoundingBox) = width(r.x)

height(r::BoundingBox) = width(r.y)

get_area(r::BoundingBox) = width(r) * height(r)

Base.isempty(r::BoundingBox) = isempty(r.x) || isempty(r.y)

midpoint(r::BoundingBox) = (midpoint(r.x), midpoint(r.y))

function Base.:(∩)(r1::BoundingBox, r2::BoundingBox)
    I = r1.x ∩ r2.x
    J = r1.y ∩ r2.y
    if isempty(I) || isempty(J)
        return InvalidBoundingBox
    else
        return BoundingBox(I, J)
    end
end

function Base.:(∪)(r1::BoundingBox, r2::BoundingBox)
    I = r1.x ∪ r2.x
    J = r1.y ∪ r2.y
    return BoundingBox(I, J)
end

Base.in(r1::BoundingBox, r2::BoundingBox) = (r1.x ∈ r2.x) && (r1.y ∈ r2.y)

Base.in(p::NTuple{2, <:Number}, r::BoundingBox) = BoundingBox(p) ∈ r

is_touching(r1::BoundingBox, r2::BoundingBox) = r1.x.a == r2.x.a || r1.x.b == r2.x.b || r1.y.a == r2.y.a || r1.y.b == r2.y.b # only considers interior touching

get_bl_corner(r::BoundingBox) = (r.x.a, r.y.a)

get_tr_corner(r::BoundingBox) = (r.x.b, r.y.b)

function diametral_circle(p, q)
    center = midpoint(p, q)
    radius = dist(p, q) / 2
    return center, radius
end

function bounding_box(center::NTuple{2, <:Number}, radius::Number)
    cx, cy = getxy(center)
    return BoundingBox(cx - radius, cx + radius, cy - radius, cy + radius)
end

function bounding_box(points)
    xmin = Inf
    xmax = -Inf
    ymin = Inf
    ymax = -Inf
    for p in each_point(points)
        px, py = f64_getxy(p)
        xmin = min(xmin, px)
        xmax = max(xmax, px)
        ymin = min(ymin, py)
        ymax = max(ymax, py)
    end
    return BoundingBox(xmin, xmax, ymin, ymax)
end

function expand(box::BoundingBox, perc = 0.1)
    x = box.x
    y = box.y
    a, b = x.a, x.b
    c, d = y.a, y.b
    Δx = (b - a) * perc
    Δy = (d - c) * perc
    return BoundingBox(a - Δx, b + Δx, c - Δy, d + Δy)
end

function bounding_box(p::NTuple, q::NTuple, r::NTuple)
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    xmin, _, xmax = min_med_max(px, qx, rx)
    ymin, _, ymax = min_med_max(py, qy, ry)
    return BoundingBox(xmin, xmax, ymin, ymax)
end

struct DiametralBoundingBox
    bounding_box::BoundingBox
    edge::NTuple{2, Int}
end

function bounding_box(points, i, j)
    p, q = get_point(points, i, j)
    c, r = diametral_circle(p, q)
    bbox = bounding_box(c, r)
    edge = (Int(i), Int(j))
    return DiametralBoundingBox(bbox, edge)
end

get_bounding_box(id_bounding_box::DiametralBoundingBox) = id_bounding_box.bounding_box

get_edge(id_bounding_box::DiametralBoundingBox) = id_bounding_box.edge