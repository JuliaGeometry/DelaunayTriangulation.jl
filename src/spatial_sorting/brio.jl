function assign_rounds(indices; rng=Random.default_rng())
    n = length(indices)
    nr = ceil(Int, log2(n)) + 1
    participant_map = Vector{Vector{Int}}(undef, nr)
    entered = falses(n)
    for i in 1:(nr-1)
        participants = Int[]
        sizehint!(participants, ceil(Int, (1 / 2)^i * n))
        for j in indices
            entered[j] && continue
            if rand(rng) < 1 / 2
                push!(participants, j)
                entered[j] = true
            end
        end
        participant_map[i] = participants
    end
    # Put all the remaining participants into the last round
    participant_map[end] = filter(i -> !entered[i], indices)
    # Reverse the rounds so that those in the first round are first
    reverse!(participant_map)
    return participant_map
end

function hilbert_sort_rounds(points; capacity=2, splitter=Middle(), rng=Random.default_rng())
    if splitter === Median()
        throw(ArgumentError("Median splitting is not yet implemented."))
    end
    indices = each_point_index(points)
    rounds = assign_rounds(indices; rng)
    pw = PointsWrapper(copy(points))
    for r in rounds 
        isempty(r) && continue 
        # We don't use a view here since `r` is a vector of indices,
        # so a view is a bit slower than just copying the points.
        sorter = HilbertSorter(pw[r]; capacity, splitter, rng, perm=r) 
        hilbert_sort!(sorter)
    end
    return rounds
end