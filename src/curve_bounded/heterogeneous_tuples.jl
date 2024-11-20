const ANY32{N} = Tuple{Any,Any,Any,Any,Any,Any,Any,Any,
    Any,Any,Any,Any,Any,Any,Any,Any,
    Any,Any,Any,Any,Any,Any,Any,Any,
    Any,Any,Any,Any,Any,Any,Any,Any,
    Any,Vararg{Any,N}}

@inline function eval_fnc_at_het_tuple_element(f::F, tup::T, idx) where {F,T}
    return _eval_fnc_at_het_tuple_element(f, idx, tup...)
end
@inline function _eval_fnc_at_het_tuple_element(f::F, idx, el::E, tup...) where {F,E}
    idx == 1 && return _eval_fnc_at_het_tuple_element(f, 1, el)
    return _eval_fnc_at_het_tuple_element(f, idx - 1, tup...)
end
@inline function _eval_fnc_at_het_tuple_element(f::F, idx, el::E) where {F,E}
    return f(el)
end

function eval_fnc_at_het_tuple_element(f::F, tup::ANY32{N}, idx) where {F,N}
    return f(tup[idx])
end

@inline function eval_fnc_at_het_tuple_two_elements(f::F, tup::T, idx1, idx2) where {F, T <: Tuple}
    return _eval_fnc_at_het_tuple_two_elements(f, idx2, tup, idx1, tup...)
end
@inline function _eval_fnc_at_het_tuple_two_elements(f::F, idx2, next_tup::T, idx1, el::E, tup...) where {F, E, T <: Tuple}
    idx1 == 1 && return _eval_fnc_at_het_tuple_two_elements(f, idx2, next_tup, 1, el)
    return _eval_fnc_at_het_tuple_two_elements(f, idx2, next_tup, idx1 - 1, tup...)
end
@inline function _eval_fnc_at_het_tuple_two_elements(f::F, idx2, next_tup::T, idx1, el::E) where {F, E, T <: Tuple}
    return _eval_fnc_at_het_tuple_two_elements(f, idx2, el, next_tup...)
end
@inline function _eval_fnc_at_het_tuple_two_elements(f::F, idx2, el::E, el2::V, tup...) where {F, E, V}
    idx2 == 1 && return _eval_fnc_at_het_tuple_two_elements(f, 1, el, el2)
    return _eval_fnc_at_het_tuple_two_elements(f, idx2 - 1, el, tup...)
end
@inline function _eval_fnc_at_het_tuple_two_elements(f::F, idx2, el::E, el2::V) where {F, E, V}
    return f(el, el2)
end

function eval_fnc_at_het_tuple_two_elements(f::F, tup::ANY32{N}, idx1, idx2) where {F, N}
    return f(tup[idx1], tup[idx2])
end

@inline function eval_fnc_at_het_tuple_element_with_arg(f::F, tup::T, arg, idx) where {F, T}
    return _eval_fnc_at_het_tuple_element_with_arg(f, idx, arg, tup...)
end
@inline function _eval_fnc_at_het_tuple_element_with_arg(f::F, idx, arg, el::E, tup...) where {F, E}
    idx == 1 && return _eval_fnc_at_het_tuple_element_with_arg(f, 1, arg, el)
    return _eval_fnc_at_het_tuple_element_with_arg(f, idx - 1, arg, tup...)
end
@inline function _eval_fnc_at_het_tuple_element_with_arg(f::F, idx, arg, el::E) where {F, E}
    return f(el, arg...)
end

function eval_fnc_at_het_tuple_element_with_arg(f::F, tup::ANY32{N}, arg, idx) where {F, N}
    return f(tup[idx], arg...)
end

@inline function eval_fnc_at_het_tuple_element_with_arg_and_prearg(f::F, tup::T, prearg, arg, idx) where {F, T}
    return _eval_fnc_at_het_tuple_element_with_arg_and_prearg(f, idx, prearg, arg, tup...)
end
@inline function _eval_fnc_at_het_tuple_element_with_arg_and_prearg(f::F, idx, prearg, arg, el::E, tup...) where {F, E}
    idx == 1 && return _eval_fnc_at_het_tuple_element_with_arg_and_prearg(f, 1, prearg, arg, el)
    return _eval_fnc_at_het_tuple_element_with_arg_and_prearg(f, idx - 1, prearg, arg, tup...)
end
@inline function _eval_fnc_at_het_tuple_element_with_arg_and_prearg(f::F, idx, prearg, arg, el::E) where {F, E}
    return f(prearg, el, arg...)
end

function eval_fnc_at_het_tuple_element_with_arg_and_prearg(f::F, tup::ANY32{N}, prearg, arg, idx) where {F, N}
    return f(prearg, tup[idx], arg...)
end

@inline function eval_fnc_in_het_tuple(tup::T, arg::A, idx) where {T, A}
    return eval_fnc_at_het_tuple_element_with_arg(self_eval, tup, arg, idx)
end

@inline self_eval(f, args...) = f(args...)