# Slow but simple implementations to compare the package's algorithms against.

export is_sink, _validate_offcenter, get_idx, is_interior_sink, _approx_ispow2
export _slow_compute_concentric_shell_ternary_split_position, slow_encroachment_test, compute_diametral_circle
export get_points_in_diametral_circle, slow_encroachment_test_diametral_lens, compute_diametral_lens
export get_points_in_diametral_lens

function slow_encroachment_test(tri::DT.Triangulation)
    E = DT.edge_type(tri)
    I = DT.integer_type(tri)
    ch = Channel{Pair{E,Tuple{Bool,I}}}(Inf) # https://discourse.julialang.org/t/can-dicts-be-threadsafe/27172/17
    @sync for i in collect(each_solid_vertex(tri))
        Base.Threads.@spawn begin
            for j in each_solid_vertex(tri)
                if i < j
                    e = DT.construct_edge(E, i, j)
                    p, q = get_point(tri, i, j)
                    r2 = 0.25((getx(p) - getx(q))^2 + (gety(p) - gety(q))^2)
                    m = (0.5(getx(p) + getx(q)), 0.5(gety(p) + gety(q)))
                    flag = false
                    k_flag = 0
                    for k in each_solid_vertex(tri)
                        if i == j || i == k || j == k
                            continue
                        end
                        r = get_point(tri, k)
                        if (getx(r) - getx(m))^2 + (gety(r) - gety(m))^2 ≤ r2
                            flag = true
                            k_flag = k
                            break
                        end
                    end
                    put!(ch, e => (flag, k_flag))
                end
            end
        end
    end
    not_in_dt_encroached_edges = Dict{E,Tuple{Bool,I}}()
    in_dt_encroached_edges = Dict{E,Tuple{Bool,I}}()
    while !isempty(ch)
        e, (b, k) = take!(ch)
        if DT.edge_exists(tri, e) || DT.edge_exists(tri, DT.reverse_edge(e))
            in_dt_encroached_edges[e] = (b, k)
        else
            not_in_dt_encroached_edges[e] = (b, k)
        end
    end
    return in_dt_encroached_edges, not_in_dt_encroached_edges
end

function slow_encroachment_test_diametral_lens(tri::DT.Triangulation, lens_angle)
    E = DT.edge_type(tri)
    I = DT.integer_type(tri)
    ch = Channel{Pair{E,Tuple{Bool,I}}}(Inf) # https://discourse.julialang.org/t/can-dicts-be-threadsafe/27172/17
    @sync for i in collect(each_solid_vertex(tri))
        Base.Threads.@spawn begin
            for j in each_solid_vertex(tri)
                if i < j
                    e = DT.construct_edge(E, i, j)
                    p, q = get_point(tri, i, j)
                    _, _, lens = compute_diametral_lens(p, q, lens_angle)
                    unique!(lens)
                    push!(lens, lens[begin])
                    boundary_nodes = [(1:length(lens)); 1]
                    if DT.polygon_features(lens, boundary_nodes)[1] < 0.0
                        reverse!(lens)
                    end
                    flag = false
                    k_flag = 0
                    for k in each_solid_vertex(tri)
                        if i == j || i == k || j == k
                            continue
                        end
                        r = get_point(tri, k)
                        if DT.distance_to_polygon(r, lens, boundary_nodes) ≥ 0
                            flag = true
                            k_flag = k
                            break
                        end
                    end
                    put!(ch, e => (flag, k_flag))
                end
            end
        end
    end
    not_in_dt_encroached_edges = Dict{E,Tuple{Bool,I}}()
    in_dt_encroached_edges = Dict{E,Tuple{Bool,I}}()
    while !isempty(ch)
        e, (b, k) = take!(ch)
        if DT.edge_exists(tri, e) || DT.edge_exists(tri, DT.reverse_edge(e))
            in_dt_encroached_edges[e] = (b, k)
        else
            not_in_dt_encroached_edges[e] = (b, k)
        end
    end
    return in_dt_encroached_edges, not_in_dt_encroached_edges
end

function _slow_compute_concentric_shell_ternary_split_position(p, q)
    m = t -> p .+ t .* (q .- p)
    L = t -> norm(m(t) .- p)
    R = t -> norm(m(t) .- q)
    tₙ = n -> exp2(n) / norm(q .- p)
    n_upp = ceil(log2(norm(q .- p)))
    obj = n -> (L(tₙ(n)) - R(tₙ(n)))^2
    n_range = -64:n_upp
    obj_vals = obj.(n_range)
    t_vals = tₙ.(n_range)
    return t_vals[argmin(obj_vals)]
end

_orient(a, b, c) = det([a[1] a[2] 1; b[1] b[2] 1; c[1] c[2] 1])

function _validate_offcenter(p, q, r, β)
    c = DT.triangle_circumcenter(p, q, r)
    offcenter = DT.triangle_offcenter(p, q, r, c, β)
    p, q, r = DT.make_shortest_edge_first(p, q, r, DT.squared_triangle_lengths_and_smallest_index(p, q, r)[4])
    if DT.squared_triangle_lengths_and_smallest_index(p, q, r)[1] ≈ DT.squared_triangle_lengths_and_smallest_index(p, q, r)[2]
        p, q, r = DT.select_shortest_edge_for_offcenter(p, q, r, c, norm(p .- q)^2)
    end
    m = (p .+ q) ./ 2

    # Check orientations 
    @test _orient(m, c, offcenter) ≈ 0.0 atol = 1.0e-6
    @test _orient(p, q, offcenter) > 0

    # Check the radius-edge ratios 
    ρ = DT.triangle_radius_edge_ratio(p, q, offcenter)
    @test ρ ≈ β
end

function get_both_offcenter_triangle_candidates(p, q, r, β)
    h1 = (2β + sqrt(4β^2 - 1)) / 2
    h2 = 1 / (2 * sqrt(4β^2 - 1))
    L1, L2, L3, idx = DT.squared_triangle_lengths_and_smallest_index(p, q, r)
    L1 = sqrt(L1)
    p, q, r = DT.make_shortest_edge_first(p, q, r, idx)
    c = DT.triangle_circumcenter(p, q, r)
    dir = c .- (p .+ q) ./ 2
    dir = dir ./ norm(dir)
    offcenter1 = (p .+ q) ./ 2 .+ dir .* h1 .* L1
    offcenter2 = (p .+ q) ./ 2 .+ dir .* h2 .* L1
    return [p, q, offcenter1], [p, q, offcenter2]
end

function reflect_point_across_line(p, q, u)
    orig_p = p
    q = collect(q .- p)
    u = collect(u .- p)
    p = zeros(2)
    pq_θ = -atan(q[2], q[1])
    cwrot_θ = [cos(pq_θ) -sin(pq_θ); sin(pq_θ) cos(pq_θ)]
    rotu = cwrot_θ * u
    # pq is now a horizontal line on the x-axis. So, to reflect rotu about this line, just do 
    rotu = [rotu[1], -rotu[2]]
    # Now, rotate back
    rotu = cwrot_θ' * rotu
    u = rotu .+ orig_p
    return Tuple(u)
end

function compute_diametral_circle(p, q)
    m = (p .+ q) ./ 2
    θ = LinRange(0, 2π, 250)
    rad = norm(p .- q) ./ 2
    circle = tuple.(m[1] .+ rad * cos.(θ), m[2] .+ rad * sin.(θ))
    lower_circle = circle[1:125, :] |> vec
    upper_circle = circle[126:end, :] |> vec
    unique!(circle)
    push!(circle, circle[begin])
    return upper_circle, lower_circle, circle
end

function compute_diametral_lens(p, q, lens_angle)
    orig_p = p
    orig_q = q
    q = collect(q .- p)
    p = zeros(2)
    pq_θ = -atan(q[2], q[1])
    cwrot_θ = [cos(pq_θ) -sin(pq_θ); sin(pq_θ) cos(pq_θ)]
    p = cwrot_θ * p
    q = cwrot_θ * q
    m = (p .+ q) ./ 2
    pq_normal = (q[2] - p[2], q[1] - p[1])
    pq_normal = pq_normal ./ norm(pq_normal)
    a = norm(p .- q) / 2
    c = a / cosd(lens_angle)
    b = c * sind(lens_angle)
    u = m .+ b .* pq_normal
    circumcenter = DT.triangle_circumcenter(p, q, u)
    circumradius = DT.triangle_circumradius(p, q, u)
    p_angle = atan(p[2] - circumcenter[2], p[1] - circumcenter[1])
    q_angle = atan(q[2] - circumcenter[2], q[1] - circumcenter[1])
    θ = LinRange(p_angle, q_angle, 250)
    lens = hcat(circumcenter[1] .+ circumradius * cos.(θ), circumcenter[2] .+ circumradius * sin.(θ))
    reflected_lens = reflect_point_across_line.(Ref(p), Ref(q), eachrow(lens))
    reflected_lens = hcat(first.(reflected_lens), last.(reflected_lens))
    combined_lens = vcat(lens, reverse(reflected_lens))
    rerotated_and_shifted_lens = (cwrot_θ' * lens')'
    rerotated_and_shifted_reflected_lens = (cwrot_θ' * reflected_lens')'
    rerotated_and_shifted_combined_lens = (cwrot_θ' * combined_lens')'
    rerotated_and_shifted_lens = tuple.(rerotated_and_shifted_lens[:, 1] .+ orig_p[1], rerotated_and_shifted_lens[:, 2] .+ orig_p[2])
    rerotated_and_shifted_reflected_lens = tuple.(rerotated_and_shifted_reflected_lens[:, 1] .+ orig_p[1], rerotated_and_shifted_reflected_lens[:, 2] .+ orig_p[2])
    rerotated_and_shifted_combined_lens = vcat(rerotated_and_shifted_lens, reverse(rerotated_and_shifted_reflected_lens))
    unique!(rerotated_and_shifted_combined_lens)
    push!(rerotated_and_shifted_combined_lens, rerotated_and_shifted_combined_lens[begin])
    return rerotated_and_shifted_lens, rerotated_and_shifted_reflected_lens, rerotated_and_shifted_combined_lens
end

function get_points_in_diametral_circle(p, q)
    tri = triangulate(rand(2, 50))
    args = DT.RefinementArguments(tri; use_lens=false)
    θ = LinRange(0, 2π, 250)
    r = LinRange(0, norm(p .- q) / 2, 250)
    m = (p .+ q) ./ 2
    points = [(m[1] + r * cos(θ), m[2] + r * sin(θ)) for r in r, θ in θ] |> vec
    classifications = DT.encroaches_upon.(Ref(p), Ref(q), points, Ref(args))
    return points[classifications], points[.!classifications]
end

function get_points_in_diametral_lens(p, q, lens_angle)
    tri = triangulate(rand(2, 50))
    args = DT.RefinementArguments(tri; use_lens=true, min_angle=lens_angle)
    θ = LinRange(0, 2π, 250)
    r = LinRange(0, norm(p .- q) / 2, 250)
    m = (p .+ q) ./ 2
    points = [(m[1] + r * cos(θ), m[2] + r * sin(θ)) for r in r, θ in θ] |> vec
    classifications = DT.encroaches_upon.(Ref(p), Ref(q), points, Ref(args))
    return points[classifications], points[.!classifications]
end

function is_sink(tri::Triangulation, T)
    _, _, t = DT.triangle_angles(get_point(tri, T...)...)
    return t ≤ π / 2 || DT.is_boundary_triangle(tri, T)
end

function is_interior_sink(tri::Triangulation, T)
    return is_sink(tri, T) && !DT.is_boundary_triangle(tri, T)
end

function get_idx(tri, p, q, r)
    i = findfirst(==(p), get_points(tri))
    j = findfirst(==(q), get_points(tri))
    k = findfirst(==(r), get_points(tri))
    return i, j, k
end

_approx_ispow2(x) = ispow2(x) || isapprox(log2(x), round(log2(x)), atol=1.0e-9, rtol=1.0e-9)
