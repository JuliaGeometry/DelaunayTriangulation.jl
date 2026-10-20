# Sampling uniformly random points from a triangle or a triangulation, e.g. `rand(rng, tri, 100)`.

struct TriangleSampler{T}
    p::NTuple{2,T}
    q::NTuple{2,T}
    r::NTuple{2,T}
end

Random.eltype(::Type{TriangleSampler{T}}) where {T} = NTuple{2,T}

function Random.rand(rng::Random.AbstractRNG, v::Random.SamplerTrivial{<:TriangleSampler{T}}) where {T}
    # https://blogs.sas.com/content/iml/2020/10/19/random-points-in-triangle.html
    itr = v[]
    p, q, r = itr.p, itr.q, itr.r
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    ax, ay = qx - px, qy - py
    bx, by = rx - px, ry - py
    u₁ = rand(rng, T)
    u₂ = rand(rng, T)
    if u₁ + u₂ > 1
        u₁ = one(T) - u₁
        u₂ = one(T) - u₂
    end
    wx = u₁ * ax + u₂ * bx
    wy = u₁ * ay + u₂ * by
    return (px + wx, py + wy)
end

Random.eltype(::Type{T}) where {P,T<:Triangulation{P}} = NTuple{2,DT.number_type(P)}

function Random.Sampler(::Type{<:Random.AbstractRNG}, tri::Triangulation, ::Random.Repetition)
    V = DT.triangle_type(tri)
    F = DT.number_type(tri)
    triangles = Vector{V}()
    areas = Vector{F}()
    sizehint!(triangles, DT.num_solid_triangles(tri))
    sizehint!(areas, DT.num_solid_triangles(tri))
    for T in DT.each_solid_triangle(tri)
        push!(triangles, T)
        i, j, k = triangle_vertices(T)
        p, q, r = get_point(tri, i, j, k)
        push!(areas, DT.triangle_area(p, q, r))
    end
    A = sum(areas)
    areas ./= A
    return Random.SamplerSimple(tri, (; tri, triangles, areas))
end

function Random.rand(rng::Random.AbstractRNG, spl::Random.SamplerSimple{<:Triangulation})
    tri = spl.data.tri
    triangles = spl.data.triangles
    areas = spl.data.areas
    multi = Multinomial(1, areas)
    T = rand(rng, multi)
    idx = findfirst(isone, T)
    i, j, k = triangle_vertices(triangles[idx])
    p, q, r = get_point(tri, i, j, k)
    return rand(rng, TriangleSampler(p, q, r))
end

function Random.rand!(rng::Random.AbstractRNG, a::AbstractVector, spl::Random.SamplerSimple{<:Triangulation})
    n = length(a)
    tri = spl.data.tri
    triangles = spl.data.triangles
    areas = spl.data.areas
    multi = Multinomial(1, areas)
    samples = rand(rng, multi, n)
    ntri = sum(samples; dims=2)
    cur_idx = 0
    for (idx, T) in enumerate(triangles)
        ntri[idx] == 0 && continue
        i, j, k = triangle_vertices(T)
        p, q, r = get_point(tri, i, j, k)
        _spl = TriangleSampler(p, q, r)
        m = ntri[idx]
        rand!(rng, view(a, (cur_idx+1):(cur_idx+m)), _spl)
        cur_idx += m
    end
    return a
end
