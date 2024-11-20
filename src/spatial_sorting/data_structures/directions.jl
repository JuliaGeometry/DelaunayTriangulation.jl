abstract type AbstractDirection end

struct Left <: AbstractDirection end

struct NotLeft <: AbstractDirection end # Right or Equal 

op(::Left) = <
op(::NotLeft) = ≥

reverse(::Left) = NotLeft()
reverse(::NotLeft) = Left()

struct HilbertOp{D<:AbstractDirection,C<:AbstractCoordinate,T}
    dir::D
    coord::C
    pivot::T
end
@inline (o::HilbertOp)(p) = op(o.dir)(get(o.coord)(p), o.pivot)

struct HilbertDirections{C<:AbstractCoordinate,XD<:AbstractDirection,YD<:AbstractDirection}
    coord::C
    xdir::XD
    ydir::YD
end