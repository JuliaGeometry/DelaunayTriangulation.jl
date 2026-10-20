# Slow reference implementations for curves, and checks for boundary enrichment.

export slow_total_absolute_curvature, slow_arc_length, slow_bezier_eval, slow_eval_bspline, slow_get_segment
export closest_point_on_curve, flatten_boundary_nodes, maximum_total_variation, all_points_are_inside
export all_diametral_circles_are_empty

function closest_point_on_curve(c, p)
    t = LinRange(0, 1, 50000)
    t′ = argmin(t) do τ
        q = c(τ)
        norm(p .- q)
    end
    return t′, c(t′)
end

function slow_arc_length(c, t₁, t₂)
    t = LinRange(t₁, t₂, 15000)
    s = 0.0
    for i in 1:(length(t)-1)
        s += norm(c(t[i+1]) .- c(t[i]))
    end
    return s
end

function slow_total_absolute_curvature(c, t₁, t₂)
    t = LinRange(t₁, t₂, 1500)
    s = 0.0
    for i in 1:(length(t)-1)
        T₁ = DT.differentiate(c, t[i])
        T₂ = DT.differentiate(c, t[i+1])
        _rat = dot(T₁, T₂) / (norm(T₁) * norm(T₂))
        _dot = _rat > 1 ? 1 : _rat < -1 ? -1 : _rat
        θ = acos(_dot)
        s += abs(θ)
    end
    return s
end

function slow_get_segment(control_points, knots, alpha, tension, t)
    L = 0.0
    for i in 2:lastindex(control_points)
        L += norm(control_points[i] .- control_points[i-1])^alpha
    end
    is_closed = control_points[begin] == control_points[end]
    s = 0
    for outer s in eachindex(knots)
        knots[s] ≤ t ≤ knots[s+1] && break
    end
    points = NTuple{2,Float64}[]
    if s == 1
        if is_closed
            push!(points, control_points[end-1])
        else
            push!(points, DT.extend_left_control_point(control_points))
        end
    else
        push!(points, control_points[s-1])
    end
    push!(points, control_points[s], control_points[s+1])
    if s == length(knots) - 1
        if is_closed
            push!(points, control_points[begin+1])
        else
            push!(points, DT.extend_right_control_point(control_points))
        end
    else
        push!(points, control_points[s+2])
    end
    return DT.catmull_rom_spline_segment(points..., alpha, tension), s
end

function slow_eval_bspline(control_points, knots, t)
    if t == 0
        return control_points[begin] |> collect
    elseif t == 1
        return control_points[end] |> collect
    end
    val = Real[0, 0]
    order = length(knots) - length(control_points)
    a, b = knots[order], knots[length(control_points)+1]
    t = a + (b - a) * t
    for i in eachindex(control_points)
        val = val .+ control_points[i] .* slow_eval_bspline_basis(knots, i, order, t)
    end
    return val
end

function slow_eval_bspline_basis(knots, i, order, t)
    if order == 1
        return (knots[i] ≤ t < knots[i+1]) ? 1.0 : 0.0
    else
        coeff1 = (t - knots[i]) / (knots[i+order-1] - knots[i])
        coeff2 = (knots[i+order] - t) / (knots[i+order] - knots[i+1])
        coeff1 = isfinite(coeff1) ? coeff1 : 0.0
        coeff2 = isfinite(coeff2) ? coeff2 : 0.0
        return coeff1 * slow_eval_bspline_basis(knots, i, order - 1, t) + coeff2 * slow_eval_bspline_basis(knots, i + 1, order - 1, t)
    end
end

function slow_bezier_eval(points, t)
    points = collect.(points)
    n = length(points) - 1
    return collect(sum([binomial(n, i) .* (1 - t)^(n - i) .* t^i .* points[i+1] for i in 0:n]))
end

function flatten_boundary_nodes(points, boundary_nodes, segments=nothing)
    _points = NTuple{2,Float64}[]
    if DT.has_multiple_curves(boundary_nodes)
        nc = DT.num_curves(boundary_nodes)
        for i in 1:nc
            _boundary_nodes = get_boundary_nodes(boundary_nodes, i)
            ns = DT.num_sections(_boundary_nodes)
            for j in 1:ns
                __boundary_nodes = get_boundary_nodes(_boundary_nodes, j)
                if !(get_boundary_nodes(__boundary_nodes, 1) isa DT.AbstractParametricCurve)
                    n = DT.num_boundary_edges(__boundary_nodes)
                    for k in 1:(n+1)
                        push!(_points, get_point(points, get_boundary_nodes(__boundary_nodes, k)))
                    end
                else
                    append!(_points, get_boundary_nodes(__boundary_nodes, 1).(LinRange(0, 1, 250)))
                end
            end
            if !(get_boundary_nodes(get_boundary_nodes(_boundary_nodes, 1), 1) isa DT.AbstractParametricCurve)
                push!(_points, get_point(points, get_boundary_nodes(get_boundary_nodes(_boundary_nodes, 1), 1)))
            end
            push!(_points, (NaN, NaN))
        end
    elseif DT.has_multiple_sections(boundary_nodes)
        ns = DT.num_sections(boundary_nodes)
        for j in 1:ns
            __boundary_nodes = get_boundary_nodes(boundary_nodes, j)
            if !(get_boundary_nodes(__boundary_nodes, 1) isa DT.AbstractParametricCurve)
                n = DT.num_boundary_edges(__boundary_nodes)
                for k in 1:n
                    push!(_points, get_point(points, get_boundary_nodes(__boundary_nodes, k)))
                end
            else
                append!(_points, get_boundary_nodes(__boundary_nodes, 1).(LinRange(0, 1, 250)))
            end
        end
        if !(get_boundary_nodes(get_boundary_nodes(boundary_nodes, 1), 1) isa DT.AbstractParametricCurve)
            push!(_points, get_point(points, get_boundary_nodes(get_boundary_nodes(boundary_nodes, 1), 1)))
        else
            push!(_points, get_boundary_nodes(get_boundary_nodes(boundary_nodes, 1), 1).(0.0))
        end
    else
        if !(get_boundary_nodes(boundary_nodes, 1) isa DT.AbstractParametricCurve)
            n = DT.num_boundary_edges(boundary_nodes)
            for k in 1:n
                push!(_points, get_point(points, get_boundary_nodes(boundary_nodes, k)))
            end
            push!(_points, get_point(points, get_boundary_nodes(boundary_nodes, 1)))
        else
            append!(_points, get_boundary_nodes(boundary_nodes, 1).(LinRange(0, 1, 250)))
        end
    end
    if !isnothing(segments)
        push!(_points, (NaN, NaN))
        for e in segments
            i, j = e
            p, q = get_point(points, i, j)
            push!(_points, p, q, (NaN, NaN))
        end
    end
    push!(_points, (NaN, NaN))
    for p in DT.each_point(points)
        push!(_points, p)
        push!(_points, (NaN, NaN))
    end
    return _points
end

function maximum_total_variation(points, boundary_nodes, boundary_curves)
    isempty(boundary_curves) && return 0.0
    max_Δθ = -Inf
    if DT.has_multiple_curves(boundary_nodes)
        ctr = 1
        for i in 1:DT.num_curves(boundary_nodes)
            curve_nodes = get_boundary_nodes(boundary_nodes, i)
            for j in 1:DT.num_sections(curve_nodes)
                section_nodes = get_boundary_nodes(curve_nodes, j)
                curve = boundary_curves[ctr]
                DT.is_piecewise_linear(curve) && continue
                for k in 1:DT.num_boundary_edges(section_nodes)
                    u, v = get_boundary_nodes(section_nodes, k), get_boundary_nodes(section_nodes, k + 1)
                    p, q = get_point(points, u, v)
                    t₁, t₂ = DT.get_inverse(curve, p), DT.get_inverse(curve, q)
                    Δθ = slow_total_absolute_curvature(curve, t₁, t₂)
                    if Δθ > max_Δθ
                        max_Δθ = Δθ
                    end
                end
                ctr += 1
            end
        end
    elseif DT.has_multiple_sections(boundary_nodes)
        for i in 1:DT.num_sections(boundary_nodes)
            section_nodes = get_boundary_nodes(boundary_nodes, i)
            curve = boundary_curves[i]
            DT.is_piecewise_linear(curve) && continue
            for j in 1:DT.num_boundary_edges(section_nodes)
                u, v = get_boundary_nodes(section_nodes, j), get_boundary_nodes(section_nodes, j + 1)
                p, q = get_point(points, u, v)
                t₁, t₂ = DT.get_inverse(curve, p), DT.get_inverse(curve, q)
                Δθ = slow_total_absolute_curvature(curve, t₁, t₂)
                if Δθ > max_Δθ
                    max_Δθ = Δθ
                end
            end
        end
    else
        curve = boundary_curves[1]
        DT.is_piecewise_linear(curve) && return 0.0
        for i in 1:DT.num_boundary_edges(boundary_nodes)
            u, v = get_boundary_nodes(boundary_nodes, i), get_boundary_nodes(boundary_nodes, i + 1)
            p, q = get_point(points, u, v)
            t₁, t₂ = DT.get_inverse(curve, p), DT.get_inverse(curve, q)
            if t₂ == 0.0
                t₂ = 1.0
            end
            Δθ = slow_total_absolute_curvature(curve, t₁, t₂)
            if Δθ > max_Δθ
                max_Δθ = Δθ
            end
        end
    end
    return max_Δθ
end

function all_diametral_circles_are_empty(enricher::DT.BoundaryEnricher)
    points = get_points(enricher)
    tree = DT.BoundaryRTree(points)
    n = 0
    for e in DT.each_boundary_edge(enricher)
        n += 1
        i, j = DT.edge_vertices(e)
        insert!(tree, i, j)
    end
    if DT.has_segments(enricher)
        for e in DT.get_segments(enricher)
            n += 1
            i, j = DT.edge_vertices(e)
            insert!(tree, i, j)
        end
    end
    rat = n
    for i in DT.each_point_index(points)
        intersects = DT.get_intersections(tree, i)
        for box in intersects
            u, v = DT.get_edge(box)
            u, v = DT.reorient_edge(enricher, u, v)
            (u == i || v == i) && continue
            cert = DT.point_position_relative_to_diametral_circle(DT.ExactKernel(), get_point(points, u, v, i)...)
            flag = DT.is_inside(cert) && DT.is_visible(DT.test_visibility(DT.ExactKernel(), enricher, u, v, i))
            if flag
                p, q = get_point(points, u, v)
                t, Δθ, ct = DT.compute_split_position(enricher, u, v, DT.ExactKernel())
                flag = flag && !isnan(Δθ)
            end
            rat -= flag
        end
    end
    return rat / n # can't just look at rat/n == 1 since, due to precision issues, we might have some points just slightly in the circle
end

function all_points_are_inside(enricher::DT.BoundaryEnricher, orig_points, orig_boundary_nodes)
    orig_points, orig_boundary_nodes = deepcopy(orig_points), deepcopy(orig_boundary_nodes)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(orig_points, orig_boundary_nodes, Int)
    new_points, new_boundary_nodes = DT.polygonise(orig_points, new_boundary_nodes, boundary_curves)
    points = get_points(enricher)
    hierarchy = DT.get_polygon_hierarchy(enricher)
    for p in DT.each_point(points)
        δ = DT.distance_to_polygon(p, new_points, new_boundary_nodes)
        δ < -1.0e-4 && return false
    end
    return true
end
