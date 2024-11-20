module DelaunayTriangulation

const DefaultAdjacentValue = 0 
const ∅ = DefaultAdjacentValue
const GhostVertex = -1
const 𝒢 = GhostVertex

@inline ε(::Type{Float64}) = 1.4901161193847656e-8 # sqrt(eps(Float64))
@inline ε(::Type{Float32}) = 0.00034526698f0 # sqrt(eps(Float32))
@inline ε(::Type{T}) where {T} = sqrt(eps(T))
@inline ε(x) = ε(number_type(x))

const INF_WARN = Ref(true)
toggle_inf_warn!() = (INF_WARN[] = !INF_WARN[])

@eval macro $(Symbol("const"))(field)
    if VERSION >= v"1.8.0-DEV.1148"
        return Expr(:const, esc(field))
    else
        return esc(field)
    end
end

macro optarg1(arg, expr)
    #=
    Transforms 
    @optarg1 arg function f(x, args...; kwargs...)
        body 
    end
    into 
    function f(x, args...; kwargs...)
        body 
    end
    @inline function f(args...; kwargs...)
        return f(arg, args...; kwargs...)
    end
    =#
    @assert expr.head == :function || expr.head == :(=) "Expected function definition"
    call = expr.args[1]
    fname = call.args[1]
    haskw = Base.isexpr(call.args[2], :parameters)
    if !haskw
        remaining_args = call.args[3:end]
        newcall = Expr(:call, fname, remaining_args...)
        newval = Expr(:call, fname, arg, remaining_args...)
    else
        remaining_args = call.args[4:end]
        newcall = Expr(:call, fname, Expr(:parameters, Expr(:(...), :(kwargs))), remaining_args...)
        newval = Expr(:call, fname, Expr(:parameters, Expr(:(...), :(kwargs))), arg, remaining_args...)
    end
    newbody = Expr(:block, Expr(:return, newval))
    newf = Expr(:function, newcall, newbody)
    inlinedf = :(@inline($newf))
    return quote
        $(esc(expr))
        $(esc(inlinedf))
    end
end

import ExactPredicates
import AdaptivePredicates
import Random
import PrecompileTools

abstract type AbstractPredicateKernel end
const PredicateCache = Union{Nothing,<:Tuple}

export getx, gety, getz, getxy, getxyz, get_point
export each_point_index, each_point, num_points
export edge_vertices
export num_edges, each_edge
export get_boundary_nodes
export geti, getj, getk, triangle_vertices
include("primitives/points.jl")
include("primitives/boundary_nodes.jl")
include("primitives/edges.jl")
include("primitives/triangle.jl")

include("triangulation/data_structures/adjacent.jl")
include("triangulation/data_structures/adjacent2vertex_candidates.jl")
include("triangulation/data_structures/boundary_edge_map.jl")
include("triangulation/data_structures/ghost_vertex_map.jl")
include("triangulation/data_structures/ghost_vertex_ranges.jl")
include("triangulation/data_structures/point_location_diagnostics.jl")
include("triangulation/data_structures/predicate_diagnostics.jl")
include("triangulation/data_structures/triangle_counts.jl")
include("triangulation/data_structures/triangulation.jl")

end