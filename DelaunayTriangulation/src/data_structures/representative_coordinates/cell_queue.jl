"""
    CellQueue{T}

A struct representing the priority queue of [`Cell`](@ref)s, used for sorting the cells in a grid
according to their maximum distance.

# Fields
- `queue::MaxPriorityQueue{Cell{T},T}`: The priority queue of cells, sorting according to maximum distance.

# Constructors
    CellQueue{T}()

Constructs a new `CellQueue` with elements of type `Cell{T}`.
"""
struct CellQueue{T}
    queue::MaxPriorityQueue{Cell{T}, T}
    function CellQueue{T}() where {T}
        return new{T}(MaxPriorityQueue{Cell{T}, T}())
    end
end

"""
    insert_cell!(queue::CellQueue, cell::Cell)

Inserts a `cell` into the `queue`.
"""
function insert_cell!(queue::CellQueue, cell::Cell)
    if !haskey(queue.queue, cell)
        queue.queue[cell] = cell.max_dist
    end
    return queue
end

"""
    get_next_cell!(queue::CellQueue)

Returns the next cell in the queue.
"""
get_next_cell!(queue::CellQueue) = popfirst!(queue.queue).first

#=
    isempty(queue::CellQueue) -> Bool 

Returns `true` if the `queue` is empty, and `false` otherwise.
=#
Base.isempty(queue::CellQueue) = Base.isempty(queue.queue)
