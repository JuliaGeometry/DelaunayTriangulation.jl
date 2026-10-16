@doc """
    abstract type AbstractPredicateKernel 

Abstract type for defining a method for computing predicates. The subtypes are:

- `FastKernel`: Uses the determinant definitions of the predicates, with no adaptivity or exact arithmetic.
- `ExactKernel`: Uses ExactPredicates.jl.
- `AdaptiveKernel`: Uses AdaptivePredicates.jl.

Please see the documentation for more information on the differences between these predicate types.
"""
AbstractPredicateKernel

"""
    FastKernel()

Pass this to predicates to declare that determinant definitions of predicates
should be used, avoiding adaptivity and exact arithmetic.

See also [`ExactKernel`](@ref) and [`AdaptiveKernel`](@ref).
"""
struct FastKernel <: AbstractPredicateKernel end

"""
    ExactKernel()

Pass this to predicates to use ExactPredicates.jl for computing predicates.

See also [`FastKernel`](@ref) and [`AdaptiveKernel`](@ref).
"""
struct ExactKernel <: AbstractPredicateKernel end

"""
    AdaptiveKernel()

Pass this to predicates to use AdaptivePredicates.jl for computing predicates.

See also [`FastKernel`](@ref) and [`ExactKernel`](@ref).
"""
struct AdaptiveKernel <: AbstractPredicateKernel end

include("predicates/certificate.jl")
include("predicates/predicate_definitions.jl")
include("predicates/predicates.jl")
include("predicates/boundaries_and_ghosts.jl")
