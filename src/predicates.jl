@doc """
    abstract type AbstractPredicateType 

Abstract type for defining a method for computing predicates. The subtypes are:

- `Fast`: Uses the determinant definitions of the predicates, with no adaptivity or exact arithmetic.
- `Exact`: Uses ExactPredicates.jl.
- `Adaptive`: Uses AdaptivePredicates.jl.

Please see the documentation for more information on the differences between these predicate types.
"""
AbstractPredicateType

"""
    Fast()

Pass this to predicates to declare that determinant definitions of predicates
should be used, avoiding adaptivity and exact arithmetic.
"""
struct Fast <: AbstractPredicateType end 

"""
    def_alg222()

Pass this to predicates to use ExactPredicates.jl for computing predicates.
"""
struct Exact <: AbstractPredicateType end 

"""
    Adaptive()

Pass this to predicates to use AdaptivePredicates.jl for computing predicates.
"""
struct Adaptive <: AbstractPredicateType end 

def_alg222() = throw("...")

include("predicates/certificate.jl")
include("predicates/predicate_definitions.jl")
include("predicates/predicates.jl")
include("predicates/boundaries_and_ghosts.jl")