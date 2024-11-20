function separate!(v, p::F; carry=nothing) where {F} # from Andy Dienes on Slack
    isempty(v) && return firstindex(v)
    i = firstindex(v)
    j = lastindex(v)
    while i < j
        pi = p(v[i])
        pj = p(v[j])
        if !pi || pj
            v[i], v[j] = v[j], v[i]
            if !isnothing(carry)
                carry[i], carry[j] = carry[j], carry[i]
            end
            !pi && (j -= 1)
            pj && (i += 1)
        else
            i += 1
            j -= 1
        end
    end
    return i - !p(v[i])
end

function triseparate!(v, p::F; carry=nothing) where {F}
    isempty(v) && return (firstindex(v), firstindex(v))

    i = firstindex(v)  # Start index for -1 partition
    j = firstindex(v)  # Start index for 0 partition
    k = lastindex(v)   # Start index for 1 partition

    while j ≤ k
        pj = p(v[j])

        if pj < 0
            v[i], v[j] = v[j], v[i]
            if !isnothing(carry)
                carry[i], carry[j] = carry[j], carry[i]
            end
            i += 1
            j += 1
        elseif pj > 0
            v[j], v[k] = v[k], v[j]
            if !isnothing(carry)
                carry[j], carry[k] = carry[k], carry[j]
            end
            k -= 1
        else  # pj == 0
            j += 1
        end
    end
    return i, j
end

function get_median!(v, p::F; rng=Random.default_rng(), carry=nothing) where {F} # Could also implement http://erdani.org/research/sea2017.pdf
    ℓ = length(v)
    return quickselect!(v, p, (ℓ + 1) ÷ 2; rng=rng, carry=carry)
end

struct QuickSelectPredicate{F,T}
    p::F
    pivot::T
end
(q::QuickSelectPredicate)(x) = q.p(x, q.pivot)

function quickselect!(v, p::F, k; rng=Random.default_rng(), carry=nothing) where {F}
    length(v) == 1 && return only(v) # assumes k == 1 
    @show v, k
    pivot = rand(rng, v)
    i, j = triseparate!(v, QuickSelectPredicate(p, pivot); carry)
    lows = view(v, firstindex(v):i-1)
    highs = view(v, j:lastindex(v))
    pivots = view(v, i:j-1)
    if k ≤ length(lows)
        return quickselect!(lows, p, k; rng, carry=isnothing(carry) ? carry : view(carry, firstindex(v):i))
    elseif k ≤ length(lows) + length(pivots)
        return first(pivots)
    else
        return quickselect!(highs, p, k - length(lows) - length(pivots); rng, carry=isnothing(carry) ? carry : view(carry, i+1:lastindex(v)))
    end
end