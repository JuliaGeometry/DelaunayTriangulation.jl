# Taken and modified from https://github.com/JuliaGeometry/DistMesh.jl/blob/4f0a64a44bd807e2aef7bdfdadcf3a33386619ee/src/hilbertsort.jl#L1C3-L1C4 with permission from sjkelly.
# Copyright (c) 2019: Steve Kelly (DistMesh.jl)
# Copyright (c) 2014: Ariel Keselman (Original from GeometricalPredicates.jl)

using DelaunayTriangulation, CairoMakie, EnumX

############################################
## Drawing Hilbert curves
function apply_transform(p, transform)
    # A transform is specified as a 6-tuple (a, b, c, d, e, f), representing the transformation
    # p -> Mp + t, where M = [a b; c d] and t = [e f].
    a, b, c, d, e, f = transform
    x, y = getxy(p)
    x′ = x * a + y * b + e 
    y′ = x * c + y * d + f
    return (x′, y′)
end
function apply_transform!(points, transform)
    for i in DelaunayTriangulation.each_point_index(points)
        p = get_point(points, i)
        DelaunayTriangulation.set_point!(points, i, apply_transform(p, transform))
    end
    return points
end
const SW_TRANSFORM = (0//1, 1//2, 1//2, 0//1, -1//4, -1//4) # Scale -> Translate -> Rotate -> Mirror
const SE_TRANSFORM = (0//1, -1//2, -1//2, 0//1, 1//4, -1//4) # Scale -> Translate -> Rotate -> Mirror
const NW_TRANSFORM = (1//2, 0//1, 0//1, 1//2, -1//4, 1//4) # Scale -> Translate
const NE_TRANSFORM = (1//2, 0//1, 0//1, 1//2, 1//4, 1//4) # Scale -> Translate
function hilbert_step(points)
    sw, se, nw, ne = copy(points), copy(points), copy(points), copy(points)
    apply_transform!(sw, SW_TRANSFORM)
    apply_transform!(se, SE_TRANSFORM)
    apply_transform!(nw, NW_TRANSFORM)
    apply_transform!(ne, NE_TRANSFORM)
    append!(sw, nw, ne, se)
    return sw
end
function hilbert(order; points = [(0.0, 0.0)]) # use points kwarg if you want to look at the transformed points
    # This function assumes that the Hilbert curve lives in [-1/2, 1/2]².
    # Since we centre at (0.0, 0.0), we don't have to keep track of the width,
    # allowing the transforms above to be constant.
    for _ in 1:order
        points = hilbert_step(points)
    end
    return points
end
function hilbert(order, lo, hi; points = [(0.0, 0.0)]) # kwarg is in [-1/2, 1/2]² space
    points = hilbert(order; points = points)
    width = hi - lo
    for i in DelaunayTriangulation.each_point_index(points)
        p = get_point(points, i)
        x, y = getxy(p)
        x′ = lo + width * (x + 1/2)
        y′ = lo + width * (y + 1/2)
        DelaunayTriangulation.set_point!(points, i, (x′, y′))
    end
    return points
end

function philbert(p::Int)
    # This curve will span [0, 2^p]².
    hilb = hilbert(p, -1//2, -1//2 + (1 << p); points = [(0//1, 0//1)])
    return convert(Vector{NTuple{2, Int}}, hilb)
end

function plot_hilbert(order, lo = -1/2, hi = 1/2)
    points = hilbert(order, lo, hi)
    width = hi - lo
    steps = LinRange(lo, hi, order + 2)
    fig, ax, sc = lines(points, color = :blue, linewidth = 2)
    vlines!(ax, steps, color = :black, linewidth = 1)
    hlines!(ax, steps, color = :black, linewidth = 1)
    xlims!(ax, lo, hi)
    ylims!(ax, lo, hi)
    display(fig)
    return fig, ax, sc 
end
function plot_philbert(p)
    hilb = philbert(p)
    fig, ax, sc = lines(hilb, color = :blue, linewidth = 2)
    xlims!(ax, -1e-2, (1 << p) + 1e-2)
    ylims!(ax, -1e-2, (1 << p) + 1e-2)
    display(fig)
    return fig, ax, sc
end

############################################

const forward = true
const backward = false
next2d(c) = c % 2 + 1
function hilbertsort!(a)
    return hilbertsort!(backward, backward, 2, a, 1, length(a), 1)
end
function hilbertsort!(directionx, directiony, coordinate, a, lo, hi, lim)
    hi - lo ≤ lim && return a

    i2 = (lo + hi) >>> 1
    i1 = (lo + i2) >>> 1
    i3 = (i2 + hi) >>> 1

    select!(directionx, coordinate, a, i2, lo, hi)
    select!(directiony, next2d(coordinate), a, i1, lo, i2)
    select!(!directiony, next2d(coordinate), a, i3, i2, hi)

    hilbertsort!(directionx, directiony, next2d(coordinate), a, lo, i1, lim)
    hilbertsort!(directionx, directiony, coordinate, a, i1, i2, lim)
    hilbertsort!(directionx, directiony, coordinate, a, i2, i3, lim)
    hilbertsort!(!directiony, !directionx, next2d(coordinate), a, i3, hi, lim)

    return a
end

function select!(direction, coord, v, k, lo, hi)
    if direction === forward
        while lo < hi
            if isone(hi - lo)
                if v[hi][coord] < a[lo][coord]
                    v[lo], v[hi] = v[hi], v[lo]
                end
                return
            end
            pivot = v[(lo+hi)>>>1]
            i, j = lo, hi
            while true
                pivot_elt = pivot[coord]
                while v[i][coord] < pivot_elt
                    i += 1
                end
                while pivot_elt < v[j][coord]
                    j -= 1
                end
                i ≤ j || break
                v[i], v[j] = v[j], v[i]
                i += 1
                j -= 1
            end
            if k ≤ j
                hi = j
            elseif i ≤ k
                lo = i
            else
                return
            end
        end
    else
        while lo < hi
            if isone(hi - lo)
                if v[hi][coord] > v[lo][coord]
                    v[lo], v[hi] = v[hi], v[lo]
                end
                return
            end
            pivot = v[(lo+hi)>>>1]
            i, j = lo, hi
            while true
                pivot_elt = pivot[coord]
                while v[i][coord] > pivot_elt
                    i += 1
                end
                while pivot_elt > v[j][coord]
                    j -= 1
                end
                i ≤ j || break
                v[i], v[j] = v[j], v[i]
                i += 1
                j -= 1
            end
            if k ≤ j
                hi = j
            elseif i ≤ k
                lo = i
            else
                return
            end
        end
    end
    return
end