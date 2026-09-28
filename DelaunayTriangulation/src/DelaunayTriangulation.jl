module DelaunayTriangulation

include("setup.jl")

import ExactPredicates as EP
import AdaptivePredicates as AP
import EnumX
import Random

abstract type AbstractPredicateKernel end # needs to be defined early for use in data_structures.jl
const PredicateCacheType = Union{Nothing, <:Tuple}

include("data_structures.jl")
include("predicates.jl")
include("geometric_primitives.jl")
include("utils.jl")
include("algorithms.jl")
include("validation.jl")
include("exports.jl")
include("public.jl")

end
