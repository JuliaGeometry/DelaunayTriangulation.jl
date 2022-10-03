# DelaunayTriangulation

- [DelaunayTriangulation](#delaunaytriangulation)
- [Triangulation Structure](#triangulation-structure)
- [Visualisation](#visualisation)
- [Custom Points, Edges, and Triangles](#custom-points-edges-and-triangles)
  - [Application: Cells with labels](#application-points-with-labels)
  
Package for computing the Delaunay triangulation for a given set of points, using randomised incremental insertion. The method used is that given in the second edition of the book "Computational Geometry: Algorithms and Applications" by de Berg et al. (1999). All the necessary geometrical predicates are computed rigorously using the great package [ExactPredicates.jl](https://github.com/lairez/ExactPredicates.jl), and support is provided for custom points, triangles, edges, and data types other than `Float64` and `Int64` (with a slight caveat for different float types; see the end of this README).

The main function to be used is `triangulate`, which takes as input a set of points. These points can be provided as a vector of 2-vectors, or a vector of `Point`s, where `Point <: AbstractPoint` (defined in this package). For example,

```julia
using DelaunayTriangulation
p1 = Point(-6.88, 3.61) # Point comes from the package
p2 = Point(-6.08, -0.43)
p3 = Point(-0.3, 2.01)
p4 = Point(5.1, -1.27)
p5 = Point(6.18, 1.87)
p6 = Point(3.08, 4.43)
p7 = Point(-1.34, 4.83)
p8 = Point(-1.68, -0.77)
pts = Points(p1, p2, p3, p4, p5, p6, p7, p8) # A collection of Points - also comes from the package
DT = triangulate(pts) # You can use DelaunayTriangulation.is_delaunay(DT) to verify this is Delaunay 
```

The following would also work:

```julia
p1 = [-6.88, 3.61]
p2 = [-6.08, -0.43]
p3 = [-0.3, 2.01]
p4 = [5.1, -1.27]
p5 = [6.18, 1.87]
p6 = [3.08, 4.43]
p7 = [-1.34, 4.83]
p8 = [-1.68, -0.77]
pts = [p1, p2, p3, p4, p5, p6, p7, p8]
DT = triangulate(pts)
```

This function `triangulate` will shuffle the points randomly in-place, but you can disable this by using the keyword argument `shuffle_pts=false`.

If a new point is to be added from outside the original point set, then the original triangulation needs to use `trim = false` so that the bounding triangle used in the construction of the algorithm is not removed. For example:
```julia
pts = Points(p1, p2, p3, p4, p5, p6, p7, p8)
DTri = triangulate(pts; shuffle_pts=false, trim=false)
DTri2 = triangulate(Points(p1, p2, p3, p4, p5, p6, p7); shuffle_pts=false, trim=false)
add_point!(DTri2, [-1.68, -0.77])
# triangles(DTri2).triangles == triangles(DTri).triangles
```
The bounding triangle can be removed separately from `triangulate` by using `DelaunayTriangulation.remove_bounding_triangle!`:
```julia
DTri = triangulate(pts; shuffle_pts=false, trim=true)
DTri2 = triangulate(Points(p1, p2, p3, p4, p5, p6, p7); shuffle_pts=false, trim=false)
add_point!(DTri2, [-1.68, -0.77])
DT.remove_bounding_triangle!(DTri2)
# triangles(DTri2).triangles == triangles(DTri).triangles
```

See the tests for more examples of what can be done with the triangulation. 

# Triangulation Structure 

The result `DT` above is of type `Triangulation`, which is defined by:
```julia
struct Triangulation{A,A2V,DG,H,T,P,R}
    adjacent::A
    adjacent2vertex::A2V
    graph::DG
    history::H
    triangles::T
    points::P
    root::R
end
```
The definition of each field is as follows, letting `DT` be of type `Triangulation`:
- `adjacent`: This is a map that takes an edge `(u, v)` to a vertex `w` such that `(u, v, w)` is a positively oriented triangle. You can access this via `adjacent(DT)` or `adjacent(DT, u, v)` for the edge `(u, v)`. If `adjacent(DT, u, v) = 0`, then `(u, v)` is a boundary segment.
- `adjacent2vertex`: This is a map that takes a vertex `w` and maps it to the set of all edges `(u, v)` such that `(u, v, w)` is a positively oriented triangle. You can access this via `adjacent2vertex(DT)` or `adjacent2vertex(DT, w)` for the vertex `w`. You can use `adjacent2vertex(DT, DelaunayTriangulation.BoundaryIndex)` (currently `DelaunayTriangulation.BoundaryIndex = 0`) to access all edges on the boundary of the triangulation (in particular, the convex hull), although these edges are not given in any particular order.
- `graph`: This is the graph representation of the triangulation, defined as an undirected graph between vertices. You can access this via `graph(DT)`, and the neighbours of a point `u` can be accessed using `neighbours(DT, u)`.
- `history`: This is the history structure of the triangulation, and is useful for point location during the triangulation. You can access it using `history(DT)`. See the function `locate_triangle` to see how it's used.
- `triangles`: This is the set of triangles in the triangulation, all given as positively oriented `Triangle`s. The triangles are stored in a `Triangles` struct which collects the triangles into a `Set`. You can access the triangles using `triangles(DT)`.
- `points`: This is the point set of the triangulation, stored in a `Points` struct which puts all the points into a vector. Note that if you provided a vector of 2-vectors into `triangulate`, these points will be transformed into `Point` types. You can access the points using `points(DT)`, and a specific point `i` can be accessed via its index using `get_point(DT, i)`. Note that if `shuffle_pts=true` in `triangulate`, the ordering of the points here will be different.
- `root`: This is the root of `history(DT)`, and can be accessed using `DelaunayTriangulation.root(DT)`. You'll only really need this if you happen to use `locate_triangle` outside of the incremental algorithm.

# Visualisation

The package does not currently provide any visualisation features, but note that it is very easy to visualise a triangulation using the above structures. For example, the following code will plot a triangulation along with the convex hull, using very basic ideas. Note that this could of course be refined (e.g. the convex hull could be obtained and plotted as one line, rather than plotting many boundary segments separately).

```julia
using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1])
x = rand(100)
y = rand(100)
pts = Points([Point(x, y) for (x, y) in zip(x, y)])
DTri = triangulate(pts)
Tmat = zeros(Int64, num_triangles(DTri), 3)
for (i, T) in enumerate(triangles(DTri))
  Tmat[i, :] = [geti(T), getj(T), getk(T)]
end
pmat = zeros(num_points(DTri), 2)
for (i, p) in enumerate(points(DTri))
  pmat[i, :] = [getx(p), gety(p)]
end
poly!(ax, pmat, Tmat, strokewidth=2)
for (i, j) in adjacent2vertex(DTri, DelaunayTriangulation.BoundaryIndex)
  p = get_point(DTri, i)
  q = get_point(DTri, j)
  lines!(ax, [getx(p), getx(q)], [gety(p), gety(q)], color=:red, linewidth=5)
end
```

![A triangulation with a convex hull](https://github.com/DanielVandH/DelaunayTriangulation.jl/blob/main/test/figures/test_triangulation.png)

# Custom Points, Edges, and Triangles

The package is designed to support general definitions of points, edges, triangles, and also support integer types other than `Int64`. The definitions just have to be subtypes of `AbstractPoint`, `AbstractEdge`, and `AbstractTriangle`, respectively, and the methods required to be defined are outlined in their respective docstrings. See also the abstract primitives part of the [primitives.jl](https://github.com/DanielVandH/DelaunayTriangulation.jl/blob/main/src/primitives.jl) file. We define specific types `Point`, `Edge`, and `Triangle` in this same file that are used, but we only use these in tests.

As an example, let us implement an interface which treats points as `NTuple{2, Float32}`, edges as vectors of `Int16`, and triangles as `SVector{3, Int16}`. We will then triangulate a small set of points.
```julia
using StaticArrays
struct CustomEdge <: DelaunayTriangulation.AbstractEdge{Int16}
    uv::Vector{Int16}
    CustomEdge(i::Int16, j::Int16) = new([i, j])
    CustomEdge(i, j) = CustomEdge(Int16(i), Int16(j))
end
DelaunayTriangulation.initial(e::CustomEdge) = first(e.uv)
DelaunayTriangulation.terminal(e::CustomEdge) = last(e.uv)

struct CustomTriangle <: DelaunayTriangulation.AbstractTriangle{Int16}
    ijk::SVector{3,Int16}
    CustomTriangle(i, j, k) = new(SVector{3,Int16}(i, j, k))
    CustomTriangle(ijk::NTuple{3,Int16}) = CustomTriangle(ijk[1], ijk[2], ijk[3])
end
DelaunayTriangulation.geti(T::CustomTriangle) = T.ijk[1]
DelaunayTriangulation.getj(T::CustomTriangle) = T.ijk[2]
DelaunayTriangulation.getk(T::CustomTriangle) = T.ijk[3]

struct CustomPoint <: DelaunayTriangulation.AbstractPoint{Float32,NTuple{2,Float32}}
    coords::NTuple{2,Float32}
    CustomPoint(x::Float32, y::Float32) = new((x, y))
    CustomPoint(x, y) = CustomPoint(Float32(x), Float32(y))
    CustomPoint(xy::NTuple{2,T}) where {T} = CustomPoint(xy[1], xy[2])
end
DelaunayTriangulation.getx(p::CustomPoint) = p.coords[1]
DelaunayTriangulation.gety(p::CustomPoint) = p.coords[2]

p1 = CustomPoint(2.0, 0.0)
p2 = CustomPoint(2.8, -1.0)
p3 = CustomPoint(-1.37, 3.5)
p4 = CustomPoint(5.81, 3.0)
p5 = CustomPoint(17.1, 10.2)
pts = Points(p1, p2, p3, p4, p5)
DTri = triangulate(pts;
    IntegerType=Int16, EdgeType=CustomEdge, TriangleType=CustomTriangle)
```

It is important to note that the use of `Float32` is not currently being used, since the predicates used all require a conversion to `Float64`. If you do want to compute the predicates using a type other than `Float64`, you will have to sacrifice the guarantee that the predicates are exact and define new methods for `ExactPredicates.orient` and `ExactPredicates.incircle`. (Note that converting to `Float64` already loses some of this guarantee, anyway.) The following retriangulates using new definitions of the predicates.

```julia
function ExactPredicates.orient(p::CustomPoint, q::CustomPoint, r::CustomPoint)
    detval = (getx(p) - getx(r)) * (gety(q) - gety(r)) - (gety(p) - gety(r)) * (getx(q) - getx(r))
    return sign(detval)
end
function ExactPredicates.incircle(p::CustomPoint, q::CustomPoint, r::CustomPoint, s::CustomPoint)
    val11 = getx(p) - getx(s)
    val21 = getx(q) - getx(s)
    val31 = getx(r) - getx(s)
    val12 = gety(p) - gety(s)
    val22 = gety(q) - gety(s)
    val32 = gety(r) - gety(s)
    val13 = (getx(p) - getx(s))^2 + (gety(p) - gety(s))^2
    val23 = (getx(q) - getx(s))^2 + (gety(q) - gety(s))^2
    val33 = (getx(r) - getx(s))^2 + (gety(r) - gety(s))^2
    return sign(val11 * val22 * val33 -
           val11 * val23 * val32 -
           val12 * val21 * val33 +
           val12 * val23 * val31 +
           val13 * val21 * val32 -
           val13 * val22 * val31)
end

pts = Points(p1, p2, p3, p4, p5)
DTri = triangulate(pts;
    IntegerType=Int16, EdgeType=CustomEdge, TriangleType=CustomTriangle, shuffle_pts=false)
```

## Application: Points with labels 

Let us give an example that takes the above ideas with custom points and applies it to a problem where the points we want to triangulate have labels. In particular, suppose each point represents a cell that is either physical, non-physical (a 'ghost' cell), or a boundary cell. One possible way to implement such a point is:
```julia
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
```
Here we are defining separate structs to define the labels, and a struct `Cell` that stores the point's coordinates and its label.

To now compute the triangulations, we unfortunately need to redefine a new method for `Points` since we currently use coordinates for the bounding triangle's vertices (see [#5](https://github.com/DanielVandH/DelaunayTriangulation.jl/issues/5)). The following function defines the new method, with the only change being that we now define the bounding triangle's coordinates as being ghost cells.
```julia
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
```

This is all that we need to do. Let's generate a set of cells and triangulate, colouring each point by its label.
```julia

n = 100
labels = [Physical, Ghost, Boundary]
Random.seed!(29291)
cell_x, cell_y, cell_k = 10randn(n), 10randn(n), rand(eachindex(labels), n)
cells = convert(Vector{Cell{L} where L}, Cell.(cell_x, cell_y, getindex.(Ref(labels), cell_k)))
pts = Points(cells)
DTri = triangulate(pts; shuffle_pts=false)

using CairoMakie
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
```


![A triangulation of cells coloured by label](https://github.com/DanielVandH/DelaunayTriangulation.jl/blob/main/test/figures/cell_application.png)

There appear to be some slight overhead to this approach (about 1.4x slower, see the benchmark below), for reasons that are as of yet unclear, but it is probably still faster than working with standard points and then mapping back onto the correct labels after the triangulation.

```julia
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
```
```julia
julia> cell_timings
BenchmarkTools.Trial: 119 samples with 1 evaluation.
 Range (min … max):  33.847 ms … 69.224 ms  ┊ GC (min … max): 0.00% … 38.12%
 Time  (median):     40.107 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   42.061 ms ±  7.878 ms  ┊ GC (mean ± σ):  5.18% ± 10.41%

    █ ▄▄ ▅ ▂ ▄ ▇ 
  ▆▃█▆██▇███▇█▇█▇▆▅▅▅▅▆▃▅▁▁▃▁▁▁▁▁▁▁▁▁▁▃▁▁▁▃▃▁▅▁▅▃▁▁▁▃▁▁▁▁▃▁▃▃ ▃
  33.8 ms         Histogram: frequency by time        68.5 ms <

 Memory estimate: 19.39 MiB, allocs estimate: 190911.

julia> standard_timings
BenchmarkTools.Trial: 169 samples with 1 evaluation.
 Range (min … max):  21.595 ms … 58.689 ms  ┊ GC (min … max): 0.00% … 44.66%
 Time  (median):     28.298 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   29.674 ms ±  7.227 ms  ┊ GC (mean ± σ):  6.13% ± 11.89%

    ▁▁▂▁  ▇ ▃▄▄▅█
  ▆▃████▆██▇█████▇▆▄▃▄▃▄▄▁▁▁▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃▄▁▃▁▃▅▃▁▁▁▁▃▁▁▃▁▃ ▃
  21.6 ms         Histogram: frequency by time        56.4 ms <

 Memory estimate: 16.24 MiB, allocs estimate: 89960.
 ```

