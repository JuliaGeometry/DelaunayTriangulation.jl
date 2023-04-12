function split_all_encroached_segments!(tri::Triangulation, queue, events, targets, rng=Random.default_rng(), maxiters=100)
    iters = 0
    while !encroachment_queue_is_empty(queue) && !compare_points(targets, num_points(tri))
        iters += 1
        @show iters
        e = encroachment_dequeue!(queue)
        if !edge_exists(tri, e) && !edge_exists(tri, reverse_edge(e))
            continue
        end
        split_subsegment!(tri, queue, events, targets, e, rng)
        iters ≥ maxiters && break
    end
    return nothing
end