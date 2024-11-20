# Based on https://github.com/JuliaCollections/DataStructures.jl/blob/master/src/priorityqueue.jl, "Introduction to Algorithms" (4e) by Cormen et al. (2022),
#    and https://stevenschmatz.gitbooks.io/data-structures-and-algorithms/content/281/lecture_11.html.

struct MaxPriorityQueue{K, V} <: AbstractDict{K, V}
    data::Vector{Pair{K, V}}
    map::Dict{K, Int}
end
function MaxPriorityQueue{K, V}() where {K, V}
    return MaxPriorityQueue{K, V}(Pair{K, V}[], Dict{K, Int}())
end
function MaxPriorityQueue(data::Dict{K, V}) where {K, V}
    queue = MaxPriorityQueue{K, V}()
    for pair in data
        push!(queue, pair)
    end
    return queue
end

@inline hparent(k) = k >> 1 # ⌊k/2⌋
@inline hleft(k) = k << 1 # 2k
@inline hright(k) = (k << 1) + 1 # 2k + 1
@inline hchildren(k) = (hleft(k), hright(k))
@inline is_root(k) = isone(k)
@inline Base.length(queue::MaxPriorityQueue) = length(queue.data)
@inline has_children(queue::MaxPriorityQueue, k) = hleft(k) ≤ length(queue)
@inline Base.isempty(queue::MaxPriorityQueue) = isempty(queue.data)
@inline Base.haskey(queue::MaxPriorityQueue, key) = haskey(queue.map, key)
@inline Base.first(queue::MaxPriorityQueue) = first(queue.data)

@inline function swap!(queue::MaxPriorityQueue, i, j)
    # i and j are indices in the data vector
    xⱼ = queue.data[j]
    queue.map[xⱼ.first] = i
    queue.data[i] = xⱼ
    return queue
end

function fix_up!(queue::MaxPriorityQueue, k::Int)
    node = queue.data[k]
    while !is_root(k)
        j = hparent(k)
        parent = queue.data[j]
        if parent.second < node.second
            swap!(queue, k, j)
            k = j
        else
            break
        end
    end
    queue.map[node.first] = k
    queue.data[k] = node
    return queue
end

function fix_down!(queue::MaxPriorityQueue, k::Int)
    node = queue.data[k]
    while has_children(queue, k)
        ℓ, r = hchildren(k)
        if r > length(queue) || queue.data[r].second < queue.data[ℓ].second
            j = ℓ
        else
            j = r
        end
        child = queue.data[j]
        if node.second < child.second
            swap!(queue, k, j)
            k = j
        else
            break
        end
    end
    queue.map[node.first] = k
    queue.data[k] = node
    return queue
end

function Base.setindex!(queue::MaxPriorityQueue{K, V}, priority, key) where {K, V}
    new_pair = Pair{K, V}(key, priority)
    if haskey(queue, key)
        idx = queue.map[key]
        orig_priority = queue.data[idx].second
        queue.data[idx] = new_pair
        if orig_priority < priority
            fix_up!(queue, idx)
        else
            fix_down!(queue, idx)
        end
    else
        push!(queue, new_pair)
    end
    return priority
end

function Base.push!(queue::MaxPriorityQueue, pair::Pair)
    key = pair.first
    if haskey(queue, key)
        throw(ArgumentError("Key $key already exists in queue."))
    end
    push!(queue.data, pair)
    queue.map[key] = length(queue)
    fix_up!(queue, length(queue))
    return queue
end

@inline function Base.popfirst!(queue::MaxPriorityQueue)
    root = queue.data[begin]
    last = pop!(queue.data)
    if !isempty(queue)
        queue.data[begin] = last # by putting the last node into the root, we violate the heap order and thus fix_down! makes all the changes needed
        queue.map[last.first] = firstindex(queue.data)
        fix_down!(queue, firstindex(queue.data))
    end
    delete!(queue.map, root.first)
    return root
end

@inline function Base.getindex(queue::MaxPriorityQueue, key)
    return queue.data[queue.map[key]].second
end

@inline function Base.iterate(queue::MaxPriorityQueue)
    isempty(queue) && return nothing
    state = MaxPriorityQueue(copy(queue.data), copy(queue.map))
    return popfirst!(state), state
end

@inline function Base.iterate(queue::MaxPriorityQueue, state::MaxPriorityQueue)
    isempty(state) && return nothing
    return popfirst!(state), state
end
