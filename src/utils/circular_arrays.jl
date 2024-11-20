@inline is_circular(A) == isempty(A) || (A[begin] == A[end])

function circular_equality(A, B, by::F=isequal) where {F}
    @assert is_circular(A) && is_circular(B) "A and B must be circular"
    length(A) ≠ length(B) && return false
    isempty(A) && return true # isempty(B) is true as well because of the previous assertion 
    _A = @views A[begin:(end-1)]
    _B = @views B[begin:(end-1)]
    same_idx = findfirst(Base.Fix1(by, _A[begin]), _B)
    same_idx === nothing && return false
    n = length(_A)
    for (i, a) in pairs(_A)
        j = mod1(i + same_idx - 1, n)
        b = _B[j]
        !by(a, b) && return false
    end
    return true
end

@inline nextindex_circular(C, i) = i == lastindex(C) ? firstindex(C) : i + 1

@inline previndex_circular(C, i) = i == firstindex(C) ? lastindex(C) : i - 1
