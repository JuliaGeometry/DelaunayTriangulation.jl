@doc """
    abstract type AbstractPredicateKernel 

Abstract type for defining a method for computing predicates. The subtypes are:

- `Fast`: Uses the determinant definitions of the predicates, with no adaptivity or exact arithmetic.
- `Exact`: Uses ExactPredicates.jl.
- `Adaptive`: Uses AdaptivePredicates.jl.

Please see the documentation for more information on the differences between these predicate types.
"""
AbstractPredicateKernel

"""
    Fast()

Pass this to predicates to declare that determinant definitions of predicates
should be used, avoiding adaptivity and exact arithmetic.

See also [`Exact`](@ref) and [`Adaptive`](@ref).
"""
struct Fast <: AbstractPredicateKernel end 

"""
    Exact()

Pass this to predicates to use ExactPredicates.jl for computing predicates.

See also [`Fast`](@ref) and [`Adaptive`](@ref).
"""
struct Exact <: AbstractPredicateKernel end 

"""
    Adaptive()

Pass this to predicates to use AdaptivePredicates.jl for computing predicates.

See also [`Fast`](@ref) and [`Exact`](@ref).
"""
struct Adaptive <: AbstractPredicateKernel end 

include("predicates/certificate.jl")
include("predicates/predicate_definitions.jl")
include("predicates/predicates.jl")
include("predicates/boundaries_and_ghosts.jl")