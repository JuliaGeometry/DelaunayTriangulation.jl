struct CellQueue{T}
    queue::MaxPriorityQueue{Cell{T}, T}
    function CellQueue{T}() where {T}
        return new{T}(MaxPriorityQueue{Cell{T}, T}())
    end
end

function insert_cell!(queue::CellQueue, cell::Cell)
    if !haskey(queue.queue, cell)
        queue.queue[cell] = cell.max_dist
    end
    return queue
end

get_next_cell!(queue::CellQueue) = (popfirst!(queue.queue).first; return queue)

Base.isempty(queue::CellQueue) = isempty(queue.queue)