"""
    CellQueue{T}

A struct representing the priority queue of [`Cell`](@ref)s, used for sorting the cells in a grid
according to their maximum distance.

# Fields
- `queue::PriorityQueue{Cell{T},T,typeof(Base.Order.Reverse)}`

The priority queue of cells.

# Constructors
    CellQueue{T}()

Constructs a new `CellQueue` with elements of type `Cell{T}`.
"""
struct CellQueue{T} # Could a heap be used for this? Duplicate keys could show up...
    queue::PriorityQueue{Cell{T},T,typeof(Base.Order.Reverse)}
    function CellQueue{T}() where {T}
        return new{T}(PriorityQueue{Cell{T},T,typeof(Base.Order.Reverse)}(Base.Order.Reverse))
    end
end
@inline function insert_cell!(queue::CellQueue, cell::Cell)
    return cell ∉ keys(queue.queue) && enqueue!(queue.queue, cell, cell.max_dist)
end
@inline get_next_cell!(queue::CellQueue) = dequeue!(queue.queue)
@inline Base.isempty(queue::CellQueue) = Base.isempty(queue.queue)