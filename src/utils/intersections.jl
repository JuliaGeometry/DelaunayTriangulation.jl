function intersection_of_ray_with_bounding_box(p, q, a, b, c, d)
    px, py = getxy(p)
    qx, qy = getxy(q)
    pℓbx, pℓby = a, c
    pℓtx, pℓty = a, d
    prtx, prty = b, d
    prbx, prby = b, c
    θ = mod2pi(atan(qy - py, qx - px))
    θlb = mod2pi(atan(pℓby - py, pℓbx - px))
    θlt = mod2pi(atan(pℓty - py, pℓtx - px))
    θrt = mod2pi(atan(prty - py, prtx - px))
    θrb = mod2pi(atan(prby - py, prbx - px))
    if iszero(θ)
        return b, py
    elseif θ == oftype(θ, π / 2)
        return px, d
    elseif θ == oftype(θ, π)
        return a, py
    elseif θ == oftype(θ, 3π / 2)
        return px, c
    elseif θlb ≤ θ ≤ θrb # y = c, a ≤ x ≤ b
        # y = py + Rsinθ = c ⟹ R = (c - py) / sinθ 
        # x = px + Rcosθ = px + (c - py)cotθ
        return px + (c - py) * cot(θ), c
    elseif θrt ≤ θ ≤ θlt # y = d, a ≤ x ≤ b 
        # y = py + Rsinθ = d ⟹ R = (d - py) / sinθ
        # x = px + Rcosθ = px + (d - py)cotθ
        return px + (d - py) * cot(θ), d
    elseif θlt ≤ θ ≤ θlb # x = a, c ≤ y ≤ d
        # x = px + Rcosθ = a ⟹ R = (a - px) / cosθ
        # y = py + Rsinθ = py + (a - px)tanθ
        return a, py + (a - px) * tan(θ)
    else # x = b, c ≤ y ≤ d
        # x = px + Rcosθ = b ⟹ R = (b - px) / cosθ
        # y = py + Rsinθ = py + (b - px)tanθ
        return b, py + (b - px) * tan(θ)
    end
end

function segment_intersection_coordinates(a, b, c, d)
    ax, ay = getxy(a)
    bx, by = getxy(b)
    cx, cy = getxy(c)
    dx, dy = getxy(d)
    num = (cx - ax) * (dy - ay) - (cy - ay) * (dx - ax)
    den = (bx - ax) * (dy - cy) - (by - ay) * (dx - cx)
    α = num / den
    return ax + α * (bx - ax), ay + α * (by - ay)
end

@optarg1 DEFAULT_KERNEL function intersection_of_edge_and_bisector_ray(kernel::AbstractPredicateKernel, a, b, c; project = false)
    cert = point_position_relative_to_line(kernel, a, b, c)
    if !is_left(cert)
        if project 
            m = project_onto_line(a, b, c)
        else
            m = midpoint(a, b)
        end
        return cert, m
    else
        F = number_type(a)
        return cert, (F(NaN), F(NaN))
    end
end

@optarg1 DEFAULT_KERNEL function classify_and_compute_segment_intersection(kernel::AbstractPredicateKernel, a, b, c, d)
    F = number_type(a)
    if any(!isfinite, getxy(a)) || any(!isfinite, getxy(b)) || any(!isfinite, getxy(c)) || any(!isfinite, getxy(d))
        return None, None, None, (F(NaN), F(NaN))
    end
    cert = line_segment_intersection_type(kernel, a, b, c, d)
    cert_c = point_position_relative_to_line(kernel, a, b, c)
    cert_d = point_position_relative_to_line(kernel, a, b, d)
    if !is_none(cert)
        p = segment_intersection_coordinates(a, b, c, d)
        return cert, cert_c, cert_d, p
    else
        F = number_type(a)
        return cert, cert_c, cert_d, (F(NaN), F(NaN))
    end
end
