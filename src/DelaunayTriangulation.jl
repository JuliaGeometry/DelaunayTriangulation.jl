module DelaunayTriangulation

include("setup.jl")

@static if USE_EXACTPREDICATES
    using ExactPredicates
end
using EnumX
using Random

include("data_structures.jl")
include("geometric_primitives.jl")
include("predicates.jl")
include("utils.jl")
include("algorithms.jl")
include("validation.jl")
include("exports.jl")


end