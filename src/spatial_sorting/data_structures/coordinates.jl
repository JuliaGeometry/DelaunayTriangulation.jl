abstract type AbstractCoordinate end

struct X <: AbstractCoordinate end

struct Y <: AbstractCoordinate end

@inline next(::X) = Y()
@inline next(::Y) = X()

@inline get(::X) = getx
@inline get(::Y) = gety
