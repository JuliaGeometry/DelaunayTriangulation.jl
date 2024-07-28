# Based on https://github.com/JuliaCollections/DataStructures.jl/blob/master/src/priorityqueue.jl, "Introduction to Algorithms" (4e) by Cormen et al. (2022),
#    and https://stevenschmatz.gitbooks.io/data-structures-and-algorithms/content/281/lecture_11.html.

"""
    MaxPriorityQueue{K, V} 

Struct for a max priority queue. 

# Fields 
- `data::Vector{Pair{K, V}}`: The data of the queue, stored in a vector of key-value pairs mapping elements to their priority.
- `map::Dict{K, Int}`: A dictionary mapping elements to their index in the data vector.
"""
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

"""
    hparent(k) -> Int 

Returns the index of the parent of the element at index `k` in a heap.
"""
hparent(k) = k >> 1 # ⌊k/2⌋

"""
    hleft(k) -> Int

Returns the index of the left child of the element at index `k` in a heap.
"""
hleft(k) = k << 1 # 2k

"""
    hright(k) -> Int

Returns the index of the right child of the element at index `k` in a heap.
"""
hright(k) = (k << 1) + 1 # 2k + 1

"""
    hchildren(k) -> Tuple{Int, Int}

Returns the indices of the children of the element at index `k` in a heap.
"""
hchildren(k) = (hleft(k), hright(k))

"""
    is_root(k) -> Bool

Returns `true` if the element at index `k` is the root of a heap.
"""
is_root(k) = isone(k)

"""
    Base.length(queue::MaxPriorityQueue) -> Int

Returns the number of elements in `queue`.
"""
Base.length(queue::MaxPriorityQueue) = length(queue.data)

"""
    has_children(queue::MaxPriorityQueue, k) -> Bool

Returns `true` if the element at index `k` has children in the `queue`.
"""
has_children(queue::MaxPriorityQueue, k) = hleft(k) ≤ length(queue)

"""
    isempty(queue::MaxPriorityQueue) -> Bool

Returns `true` if the `queue` is empty.
"""
Base.isempty(queue::MaxPriorityQueue) = isempty(queue.data)

"""
    haskey(queue::MaxPriorityQueue, key) -> Bool

Returns `true` if the `queue` has an element with key `key`.
"""
Base.haskey(queue::MaxPriorityQueue, key) = haskey(queue.map, key)

"""
    first(queue::MaxPriorityQueue) -> Pair{K, V}

Returns the element with the highest priority in a `queue`, without removing it from `queue`.
"""
Base.first(queue::MaxPriorityQueue) = first(queue.data)

"""
    swap!(queue::MaxPriorityQueue, i, j) 

Swaps the elements at indices `i` and `j` in `queue`.
"""
function swap!(queue::MaxPriorityQueue, i, j)
    # i and j are indices in the data vector
    xⱼ = queue.data[j]
    queue.map[xⱼ.first] = i
    queue.data[i] = xⱼ
    return queue
end

"""
    fix_up!(queue::MaxPriorityQueue, k)

Fixes the `queue` after increasing the value of one of its elements by percolating upwards.
"""
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

"""
    fix_down!(queue::MaxPriorityQueue, k)

Fixes the `queue` after decreasing the value of one of its elements by percolating downwards.
"""
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

"""
    setindex!(queue::MaxPriorityQueue, priority, key)
    queue[key] = priority

Sets the priority of the element with key `key` in a `queue` to `priority`, or adds the element to the `queue` if it is not already present.
"""
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

"""
    push!(queue::MaxPriorityQueue, pair)

Adds the key-value pair `pair` to the `queue`.
"""
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

"""
    popfirst!(queue::MaxPriorityQueue{K, V}) where {K, V} -> Pair{K, V}

Removes and returns the element with the highest priority from the `queue`.
"""
function Base.popfirst!(queue::MaxPriorityQueue)
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

"""
    getindex(queue::MaxPriorityQueue, key) 
    queue[key]

Returns the priority of the element with key `key` in a `queue`.
"""
function Base.getindex(queue::MaxPriorityQueue, key)
    return queue.data[queue.map[key]].second
end

function Base.iterate(queue::MaxPriorityQueue)
    isempty(queue) && return nothing
    state = MaxPriorityQueue(copy(queue.data), copy(queue.map))
    return popfirst!(state), state
end
function Base.iterate(queue::MaxPriorityQueue, state::MaxPriorityQueue)
    isempty(state) && return nothing
    return popfirst!(state), state
end
