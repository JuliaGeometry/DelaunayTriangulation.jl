## The tests here are for testing a use of triangulations where each point contains some label 
abstract type AbstractLabel end
struct Physical <: AbstractLabel end
struct Ghost <: AbstractLabel end
struct Boundary <: AbstractLabel end
struct Cell{L<:AbstractLabel} <: DelaunayTriangulation.AbstractPoint{Float64,Vector{Float64}}
    x::Float64 
    y::Float64
    label::Type{L}
    Cell(x::Float64, y::Float64, label::Type{L}) where {L} = new{L}(x, y, label)
    Cell(x, y, label) = Cell(Float64(x), Float64(y), label)
    Cell{L}(x, y) where {L} = Cell(x, y, L)
end
DelaunayTriangulation.getx(p::Cell) = p.x
DelaunayTriangulation.gety(p::Cell) = p.y
getlabel(::Cell{L}) where {L} = L

function DelaunayTriangulation.Points(points::Vector{Cell{L} where L})
    xmin, xmax, ymin, ymax = typemax(Float64), typemin(Float64), typemax(Float64), typemin(Float64)
    for pt in points
        getx(pt) < xmin && (xmin = getx(pt))
        getx(pt) > xmax && (xmax = getx(pt))
        gety(pt) < ymin && (ymin = gety(pt))
        gety(pt) > ymax && (ymax = gety(pt))
    end
    width, height = xmax - xmin, ymax - ymin
    xc, yc = (xmax + xmin) / 2, (ymax + ymin) / 2
    max_width_height = max(width, height)
    lower_left = Cell{Ghost}(xc - DelaunayTriangulation.BoundingTriangleShift * max_width_height, yc - max_width_height)
    lower_right = Cell{Ghost}(xc + DelaunayTriangulation.BoundingTriangleShift * max_width_height, yc - max_width_height)
    upper_bounding = Cell{Ghost}(xc, yc + DelaunayTriangulation.BoundingTriangleShift * max_width_height)
    return Points{Float64,Vector{Float64},Cell{L} where L}(points, xc, yc, xmin, xmax, ymin,
        ymax, width, height, max_width_height, lower_left,
        lower_right, upper_bounding)
end

n = 100
labels = [Physical, Ghost, Boundary]
Random.seed!(29291)
cell_x, cell_y, cell_k = 10randn(n), 10randn(n), rand(eachindex(labels), n)
cells = convert(Vector{Cell{L} where L}, Cell.(cell_x, cell_y, getindex.(Ref(labels), cell_k)))
pts = Points(cells)
DTri = triangulate(pts; shuffle_pts=false)
@test DT.is_delaunay(DTri)

import CairoMakie: poly!, Figure, Axis, lines!, scatter!, save
fig = Figure()
ax = Axis(fig[1, 1])
Tmat = zeros(Int64, num_triangles(DTri), 3)
pmat = zeros(num_points(DTri), 2)
for (i, T) in enumerate(triangles(DTri))
    Tmat[i, :] = [geti(T), getj(T), getk(T)]
end
for (i, p) in enumerate(points(DTri))
    pmat[i, :] = [getx(p), gety(p)]
end
poly!(ax, pmat, Tmat, strokewidth=2, color=(:white, 0))
for (i, j) in adjacent2vertex(DTri, DelaunayTriangulation.BoundaryIndex)
    p = DT.get_point(DTri, i)
    q = DT.get_point(DTri, j)
    lines!(ax, [getx(p), getx(q)], [gety(p), gety(q)], color=:red, linewidth=3)
end
colors = Dict(labels .=> [:green, :blue, :magenta])
scatter!(ax, getx.(points(DTri)),
    gety.(points(DTri)),
    color=getindex.(Ref(colors), getlabel.(points(DTri))),
    markersize=11)
save("figures/cell_application.png", fig)

n = 1000
labels = [Physical, Ghost, Boundary]
cell_x, cell_y, cell_k = 10randn(n), 10randn(n), rand(eachindex(labels), n)
cells = convert(Vector{Cell{L} where L}, Cell.(cell_x, cell_y, getindex.(Ref(labels), cell_k)))
pts = Points(cells)
function cell_study(pts)
    DTri = DelaunayTriangulation.triangulate(pts; shuffle_pts=false)
    return DTri
end
standard_points = Point.(cell_x, cell_y)
standard_pts = Points(standard_points)
function standard_study(standard_pts)
    DTri = DelaunayTriangulation.triangulate(standard_pts; shuffle_pts=false)
    return DTri
end
using BenchmarkTools
cell_timings = @benchmark $cell_study($pts)
standard_timings = @benchmark $standard_study($standard_pts)

cell_timings
standard_timings