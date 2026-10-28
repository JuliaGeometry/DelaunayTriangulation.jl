# # Voronoi Treemaps
#

# In this tutorial, we show how power diagrams can be used to build _Voronoi treemaps_, which display hierarchical
# data by dividing a region into nested polygons. We first describe what treemaps are for and how the algorithm works,
# and then implement it in the form that it would have as a part of this package: a function `voronoi_treemap` that
# computes a treemap, and a plotting function `voronoitreemap`, written as a Makie recipe in the same way as Makie's
# `triplot` and `voronoiplot`. At the end, we recreate the examples of the D3 plugin
# [d3-voronoi-treemap](https://github.com/Kcnarf/d3-voronoi-treemap).
#
# ## Treemaps
#
# Hierarchical data assigns values to items that are grouped, where the groups may themselves be grouped, and so on: a
# budget is split into departments that are split into programmes, the world's economy into continents and countries,
# and a file system into folders and files. A _treemap_, introduced by
# [Shneiderman (1992)](https://doi.org/10.1145/102377.115768), shows such data by area. A region is divided into one cell
# for each node at the top of the hierarchy, with areas proportional to their values, and the cell of each group is then
# divided in the same way for its children, and so on. The leaves of the hierarchy end up with areas proportional to
# their values, and each group is the union of the cells of its leaves.
#
# The first treemaps used rectangles, which often gives long and thin cells, and the levels of the hierarchy are hard to
# tell apart since the borders of cells at different levels line up. _Voronoi treemaps_, introduced by
# [Balzer and Deussen (2005)](https://doi.org/10.1109/INFVIS.2005.1532128), instead use the cells of power diagrams
# (see [this tutorial](power.md)). These are convex polygons of any shape, which can divide any convex region, e.g. a
# circle, and which can be kept compact, so that they are easier to compare by size. The price is that the cells are
# found by an iteration, rather than directly.
#
# ## The Algorithm
#
# ### Power Diagrams
#
# Let $\Omega$ be a convex polygon with area $|\Omega|$. Given sites $s_1, \ldots, s_n$ in $\Omega$ with weights
# $w_1, \ldots, w_n$, the _power distance_ from a point $x$ to the site $s_i$ is
# ```math
# \pi_i(x) = \|x - s_i\|^2 - w_i,
# ```
# and the cell of $s_i$ in the power diagram, clipped to $\Omega$, is
# ```math
# \mathcal V_i = \{x \in \Omega : \pi_i(x) \leq \pi_j(x) \text{ for all } j\}.
# ```
# We write $A_i = |\mathcal V_i|$ and $d_{ij} = \|s_i - s_j\|$. We use four properties of these cells.
#
# 1. Each $\mathcal V_i$ is convex, since it is the intersection of $\Omega$ with the half-planes $\pi_i(x) \leq \pi_j(x)$.
# 2. The cells do not change if the same number is added to all of the weights, since only the differences $w_i - w_j$
#    appear in the inequalities $\pi_i(x) \leq \pi_j(x)$.
# 3. The border between $\mathcal V_i$ and $\mathcal V_j$ lies on the line $\pi_i(x) = \pi_j(x)$, which is perpendicular
#    to the segment from $s_i$ to $s_j$, at a distance $(d_{ij}^2 + w_i - w_j)/(2d_{ij})$ from $s_i$. Increasing $w_i$ by
#    $\delta$ thus moves the border away from $s_i$ by $\delta/(2d_{ij})$.
# 4. The site $s_j$ is in its own cell if and only if $\pi_j(s_j) \leq \pi_i(s_j)$ for all $i$, i.e. if
#    $w_i - w_j \leq d_{ij}^2$ for all $i$. If this holds for every site, then no cell is empty. Otherwise, a cell can
#    be empty, in which case its site is said to be _submerged_.
#
# ### Dividing a Polygon
#
# The main step of a treemap is to divide $\Omega$ into cells with prescribed areas, which we call a _Voronoi map_. Given
# values $v_1, \ldots, v_n > 0$, the target areas are $A_i^* = |\Omega| v_i / \sum_j v_j$, and we look for sites and
# weights for which $A_i = A_i^*$ for each $i$. We measure how far the areas are from their targets by
# ```math
# E = \frac{1}{|\Omega|}\sum_{i=1}^n |A_i - A_i^*|,
# ```
# which is twice the fraction of $\Omega$ that is in the wrong cell, as in d3-voronoi-treemap. The algorithm starts from
# random sites in $\Omega$ with zero weights, and repeats the following steps:
#
# 1. Compute the power diagram of the sites, clipped to $\Omega$, and the areas $A_i$ of its cells. Stop if $E$ is less
#    than a tolerance $\varepsilon$, or if the maximum number of iterations has been reached.
# 2. Move each site to the centroid of its cell.
# 3. Add $2(A_i^* - A_i)/\pi$ to each weight $w_i$, and then subtract the smallest weight from all of the weights.
# 4. While $w_i - w_j > d_{ij}^2$ for some $i$ and $j$, set $w_i = w_j + d_{ij}^2$.
#
# Step 2 is Lloyd's algorithm, as for [centroidal Voronoi tessellations](centroidal.md): moving the sites to the
# centroids of their cells makes the cells compact, so that the sites control the shapes of the cells and the weights
# control their areas.
#
# Step 3 adjusts the areas. By property 3, increasing $w_i$ by $\delta$ increases $A_i$ by about
# $\delta \sum_j \ell_{ij}/(2d_{ij})$, where $\ell_{ij}$ is the length of the border between $\mathcal V_i$ and
# $\mathcal V_j$. For a roughly round cell of radius $r$, whose neighbouring sites are at a distance of about $2r$, this
# is about $2\pi r \delta/(4r) = \pi\delta/2$, so adding $2(A_i^* - A_i)/\pi$ to $w_i$ would bring $A_i$ to $A_i^*$ if
# the other weights did not change. This is one step of the Jacobi method for solving the equations $A_i = A_i^*$ with
# Newton's method. Borders on the boundary of $\Omega$ do not move, so for cells on the boundary the step is too small
# rather than too large. By property 2, subtracting the smallest weight does not change the cells, and it keeps the
# weights non-negative, which we need for the anchors below.
#
# Step 4 makes sure that every site is in its own cell, by property 4, so that no cell is empty and step 2 is well
# defined. It only decreases weights: it finds the largest weights, below the current ones, that satisfy property 4, in
# the same way as the Bellman-Ford algorithm finds shortest paths, and so it takes at most $n$ passes over the pairs of
# sites.
#
# A triangulation needs three points that are not collinear, but $n$ may be 1 or 2, and the sites can be collinear,
# e.g. for a row of cells in a long and thin region. We therefore add three sites far outside of $\Omega$, which we call
# _anchors_, with zero weights, at a distance $R = 10m$ from the centre of the bounding box of $\Omega$, where $m$ is the
# larger side of the box. The power distance from a point of $\Omega$ to an anchor is then at least
# $(R - m/\sqrt 2)^2 > 85m^2$, and to one of the sites at most $2m^2$, since the weights are non-negative, so the anchors
# never get any part of $\Omega$.
#
# ### From Maps to Treemaps
#
# A Voronoi treemap divides the region for the nodes at the top of the hierarchy, where the value of a group is the sum
# of the values of the leaves below it, then divides the cell of each group for its children, and so on. By property 1,
# the cells are convex, so they can be divided in the same way as the region. Since each division is only accurate up to
# the tolerance, the areas of the leaves are only approximately proportional to their values, with errors that add up
# over the levels of the hierarchy.
#
# ### Comparison with Nocaj and Brandes (2012)
#
# The algorithm above follows the scheme of [Nocaj and Brandes (2012)](https://doi.org/10.1111/j.1467-8659.2012.03078.x),
# who alternate between moving the sites and adapting the weights, except for the weights. They treat $\sqrt{w_i}$ as
# the radius of a circle around $s_i$, multiply the radius by $A_i^*/A_i$, and keep each circle inside the cell of its
# site, by limiting the radius to the distance to the nearest site and, after moving the site, to the distance to the
# border of its cell. This also guarantees that no cell is empty, but it limits how quickly the cells can change, and
# they note that many iterations are needed when the target areas vary a lot. In our tests, a direct implementation of
# their algorithm (without the heuristic of their Section 5.2, which displaces the sites near growing cells) reached
# $E < 0.01$ within 500 iterations in fewer than half of 780 divisions, of the data in this tutorial and of four other
# data sets in six regions, while the update above reached it in all but one of them, typically within a few dozen
# iterations. d3-voronoi-treemap uses a variant of their scheme that multiplies the weights themselves by $A_i^*/A_i$.
#
# ## Implementation
#
# We now implement the algorithm in the form that it would have in this package, e.g. in a file
# `src/voronoi/treemap.jl`. As for the `clip_polygon` of [`voronoi`](@ref), a polygon is given as a `Tuple`
# `(points, boundary_nodes)`, where `boundary_nodes` lists the vertices in counter-clockwise order, with the first
# vertex repeated at the end, and cells are returned as their coordinates, as from
# [`get_polygon_coordinates`](@ref). We load the packages and the internal functions that we need.
using DelaunayTriangulation
using DelaunayTriangulation: AbstractPredicateKernel, dist_sqr, distance_to_polygon, get_centroid, get_children,
    get_parent, get_point, is_negatively_oriented, number_type, pole_of_inaccessibility, polygon_bounds, polygon_features,
    triangle_orientation
using CairoMakie
using Random
using StableRNGs
using ReferenceTests #src
using Test #src
fig_path = joinpath(@__DIR__, "..", "..", "..", "..", "assets", "docs", "figures") #src

# ### Polygons
#
# We start with a few small functions for polygons.

"""
    unit_square() -> (points, boundary_nodes)

Returns the unit square `[0, 1]²`, as a polygon for the `clip_polygon` of [`voronoi_map`](@ref) or [`voronoi_treemap`](@ref).
"""
unit_square() = ([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)], [1, 2, 3, 4, 1])

"""
    polygon_from_coordinates(coordinates) -> (points, boundary_nodes)

Returns the polygon with the given `coordinates`, which are in counter-clockwise order with the first coordinate
repeated at the end, as returned by [`get_polygon_coordinates`](@ref), in the form `(points, boundary_nodes)`.
"""
polygon_from_coordinates(coordinates) = (coordinates, [1:(length(coordinates) - 1); 1])

"""
    is_convex_polygon(points, boundary_nodes) -> Bool

Returns `true` if the polygon given by `points` and `boundary_nodes` is convex, with its vertices in counter-clockwise
order, and `false` otherwise.
"""
function is_convex_polygon(points, boundary_nodes)
    n = length(boundary_nodes) - 1
    n < 3 && return false
    for i in 1:n
        p, q, r = get_point(points, boundary_nodes[i], boundary_nodes[mod1(i + 1, n)], boundary_nodes[mod1(i + 2, n)])
        is_negatively_oriented(triangle_orientation(p, q, r)) && return false
    end
    A, _ = polygon_features(points, boundary_nodes)
    return A > 0
end

"""
    random_point_in_polygon(points, boundary_nodes, rng::AbstractRNG) -> NTuple{2, Number}

Returns a random point inside of the polygon given by `points` and `boundary_nodes`, using rejection sampling from its
bounding box.
"""
function random_point_in_polygon(points, boundary_nodes, rng::AbstractRNG)
    xmin, xmax, ymin, ymax = polygon_bounds(points, boundary_nodes)
    T = number_type(points)
    while true
        p = (xmin + (xmax - xmin) * rand(rng, T), ymin + (ymax - ymin) * rand(rng, T))
        distance_to_polygon(p, points, boundary_nodes) > 0 && return p
    end
end

"""
    get_anchors(points, boundary_nodes) -> Vector{NTuple{2, Number}}

Returns the three anchors for the polygon given by `points` and `boundary_nodes`: points at a distance of `10m` from the
centre of its bounding box, where `m` is the larger side of the box. With zero weights, and non-negative weights for
the other sites, the cells of the anchors never meet the polygon.
"""
function get_anchors(points, boundary_nodes)
    xmin, xmax, ymin, ymax = polygon_bounds(points, boundary_nodes)
    cx, cy = (xmin + xmax) / 2, (ymin + ymax) / 2
    R = 10max(xmax - xmin, ymax - ymin)
    return [(cx + R * cos(θ), cy + R * sin(θ)) for θ in (π / 2, 7π / 6, 11π / 6)]
end

# ### Voronoi Maps
#
# A Voronoi map is stored in a `VoronoiMap`, with getters in the style of the package; where the package already has a
# function for the same information about a [`VoronoiTessellation`](@ref), we add a method to it.

"""
    VoronoiMap{T}

A division of a convex polygon into cells with prescribed areas, computed by [`voronoi_map`](@ref). The cells are the
cells of a power diagram, clipped to the polygon.

## Fields
- `polygons::Vector{Vector{NTuple{2, T}}}`: The coordinates of the cells, in counter-clockwise order with the first coordinate repeated at the end. The cell of a submerged site is empty.
- `sites::Vector{NTuple{2, T}}`: The sites of the power diagram.
- `weights::Vector{T}`: The weights of the sites.
- `areas::Vector{T}`: The areas of the cells.
- `target_areas::Vector{T}`: The areas that the cells should have.
- `error::T`: The sum of the errors in the areas of the cells, divided by the area of the polygon.
- `iterations::Int`: The number of iterations that were used.
"""
struct VoronoiMap{T}
    polygons::Vector{Vector{NTuple{2, T}}}
    sites::Vector{NTuple{2, T}}
    weights::Vector{T}
    areas::Vector{T}
    target_areas::Vector{T}
    error::T
    iterations::Int
end

DelaunayTriangulation.num_polygons(vmap::VoronoiMap) = length(vmap.polygons)
DelaunayTriangulation.each_polygon_index(vmap::VoronoiMap) = eachindex(vmap.polygons)
DelaunayTriangulation.get_polygon_coordinates(vmap::VoronoiMap, i) = vmap.polygons[i]
DelaunayTriangulation.get_area(vmap::VoronoiMap, i) = vmap.areas[i]
DelaunayTriangulation.get_weight(vmap::VoronoiMap, i) = vmap.weights[i]
get_site(vmap::VoronoiMap, i) = vmap.sites[i]
get_target_area(vmap::VoronoiMap, i) = vmap.target_areas[i]
get_error(vmap::VoronoiMap) = vmap.error
num_iterations(vmap::VoronoiMap) = vmap.iterations

# Each step of the algorithm gets its own function, so that the main function reads like the algorithm.

"""
    clipped_power_diagram!(polygons, areas, sites, weights, anchors, clip_polygon, rng, predicates) -> VoronoiTessellation

Computes the power diagram of the `sites` with the `weights`, together with the `anchors` with zero weights, clipped to
`clip_polygon`, and stores the coordinates and the areas of the cells of the sites in `polygons` and `areas` (step 1).
A submerged site gets an empty polygon and a zero area.
"""
function clipped_power_diagram!(polygons, areas, sites, weights, anchors, clip_polygon, rng, predicates)
    T = eltype(weights)
    tri = triangulate([sites; anchors]; weights = [weights; zeros(T, length(anchors))], rng, predicates)
    vorn = voronoi(tri; clip = true, clip_polygon, rng, predicates)
    for i in eachindex(sites)
        if i ∈ each_polygon_index(vorn)
            polygons[i] = get_polygon_coordinates(vorn, i)
            areas[i] = get_area(vorn, i)
        else
            polygons[i] = NTuple{2, T}[]
            areas[i] = zero(T)
        end
    end
    return vorn
end

"""
    move_sites_to_centroids!(sites, vorn, areas) -> sites

Moves each site whose cell is not empty to the centroid of its cell in `vorn` (step 2).
"""
function move_sites_to_centroids!(sites, vorn, areas)
    for i in eachindex(sites)
        areas[i] > 0 && (sites[i] = get_centroid(vorn, i))
    end
    return sites
end

"""
    adjust_weights!(weights, areas, target_areas) -> weights

Adds `2(target_areas[i] - areas[i]) / π` to each weight, and then subtracts the smallest weight from all of them (step 3).
"""
function adjust_weights!(weights, areas, target_areas)
    for i in eachindex(weights)
        weights[i] += 2(target_areas[i] - areas[i]) / π
    end
    weights .-= minimum(weights)
    return weights
end

"""
    cap_weights!(weights, sites) -> weights

Decreases the weights until `weights[i] - weights[j] ≤ dist_sqr(sites[i], sites[j])` for all `i` and `j`, so that each
site is in its own cell (step 4).
"""
function cap_weights!(weights, sites)
    changed = true
    while changed # at most length(sites) + 1 passes, as for the Bellman-Ford algorithm
        changed = false
        for i in eachindex(sites), j in eachindex(sites)
            i == j && continue
            wmax = weights[j] + dist_sqr(sites[i], sites[j])
            if weights[i] > wmax
                weights[i] = wmax
                changed = true
            end
        end
    end
    return weights
end

function check_voronoi_map_args(values, clip_polygon)
    isempty(values) && throw(ArgumentError("There must be at least one value."))
    all(v -> isfinite(v) && v > 0, values) || throw(ArgumentError("The values must be positive and finite, but got $(values)."))
    is_convex_polygon(clip_polygon...) || throw(ArgumentError("The clip_polygon must be convex, with its vertices in counter-clockwise order."))
    return true
end

# The main function puts the steps together.

"""
    voronoi_map(values; clip_polygon = unit_square(), rng = Random.default_rng(), predicates = AdaptiveKernel(), tol = 1.0e-3, maxiters = 500) -> VoronoiMap

Divides a convex polygon into cells whose areas are proportional to `values`, which are the cells of a power diagram
clipped to the polygon. See the description of the algorithm above.

## Arguments
- `values`: The values of the cells, which must be positive and finite.

## Keyword Arguments
- `clip_polygon = unit_square()`: The polygon to divide, as a `Tuple` `(points, boundary_nodes)`, as for [`voronoi`](@ref). The polygon must be convex, with its vertices in counter-clockwise order.
- `rng::AbstractRNG = Random.default_rng()`: The random number generator, used for the starting sites and for the triangulations.
- `predicates::AbstractPredicateKernel = AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref).
- `tol = 1.0e-3`: The iteration stops once the sum of the errors in the areas of the cells, divided by the area of the polygon, is less than `tol`.
- `maxiters = 500`: The maximum number of iterations.

## Output
- `vmap`: The [`VoronoiMap`](@ref). If the error did not drop below `tol` within `maxiters` iterations, this is the last iterate, and `get_error(vmap) ≥ tol`.
"""
function voronoi_map(
        values; clip_polygon = unit_square(), rng::AbstractRNG = Random.default_rng(),
        predicates::AbstractPredicateKernel = AdaptiveKernel(), tol = 1.0e-3, maxiters = 500,
    )
    check_voronoi_map_args(values, clip_polygon)
    points, boundary_nodes = clip_polygon
    T = number_type(points)
    n = length(values)
    A, _ = polygon_features(points, boundary_nodes)
    V = sum(values)
    target_areas = [T(A * v / V) for v in values]
    sites = [random_point_in_polygon(points, boundary_nodes, rng) for _ in 1:n]
    weights = zeros(T, n)
    anchors = get_anchors(points, boundary_nodes)
    polygons = [NTuple{2, T}[] for _ in 1:n]
    areas = zeros(T, n)
    error, iterations = T(Inf), 0
    for iteration in 1:maxiters
        iterations = iteration
        vorn = clipped_power_diagram!(polygons, areas, sites, weights, anchors, clip_polygon, rng, predicates)
        error = sum(abs(areas[i] - target_areas[i]) for i in 1:n) / A
        (error < tol || iteration == maxiters) && break
        move_sites_to_centroids!(sites, vorn, areas)
        adjust_weights!(weights, areas, target_areas)
        cap_weights!(weights, sites)
    end
    return VoronoiMap(polygons, sites, weights, areas, target_areas, error, iterations)
end

# ### Hierarchies
#
# We represent a hierarchy as a vector of nodes, where a node is either a leaf, `label => value`, whose value is a
# positive number, or a group, `label => children`, whose children are a non-empty vector of nodes. The value of a
# group is the sum of the values of its children, and the labels can be anything.

is_group(node) = last(node) isa AbstractVector
node_value(node) = is_group(node) ? sum(node_value, last(node)) : last(node)

function check_hierarchy(nodes)
    (nodes isa AbstractVector && !isempty(nodes)) || throw(ArgumentError("A hierarchy must be a non-empty vector of nodes, but got $(repr(nodes))."))
    for node in nodes
        node isa Pair || throw(ArgumentError("A node must be a Pair, label => value or label => children, but got $(repr(node))."))
        if is_group(node)
            check_hierarchy(last(node))
        else
            value = last(node)
            (value isa Real && isfinite(value) && value > 0) || throw(ArgumentError("The value of $(repr(first(node))) must be a positive and finite number, but got $(repr(value))."))
        end
    end
    return true
end

# Hierarchical data is often stored in a table instead, with a row for each node that gives its label, the label of its
# parent, and its value, as for e.g. Plotly's treemaps. The function below converts such a table into nested pairs.

"""
    hierarchy_from_table(labels, parents, values; root = "") -> Vector{Pair}

Converts hierarchical data given as a table, with a row for each node that gives its label, the label of its parent,
and its value, into the nested pairs used by [`voronoi_treemap`](@ref). The labels must be unique, the parent of a node
at the top of the hierarchy is `root`, and the values of groups are ignored, since they are the sums of the values of
their children.
"""
function hierarchy_from_table(labels, parents, values; root = "")
    allunique(labels) || throw(ArgumentError("The labels must be unique."))
    children = Dict(label => Int[] for label in labels)
    top = Int[]
    for (i, parent) in enumerate(parents)
        if parent == root
            push!(top, i)
        elseif parent == labels[i] || !haskey(children, parent)
            throw(ArgumentError("The parent $(repr(parent)) of $(repr(labels[i])) must be another label or the root, $(repr(root))."))
        else
            push!(children[parent], i)
        end
    end
    node(i) = isempty(children[labels[i]]) ? labels[i] => values[i] : labels[i] => [node(j) for j in children[labels[i]]]
    return Pair{eltype(labels), Any}[node(i) for i in top]
end

# ### Voronoi Treemaps
#
# A treemap is stored in a `VoronoiTreemap`, with the nodes in breadth-first order, so that the nodes at the top of the
# hierarchy come first. Its getters again extend those of the package where they exist, e.g. `get_parent` and
# `get_children`, which the package uses for its polygon hierarchies.

"""
    VoronoiTreemap{L, T, P}

A Voronoi treemap, computed by [`voronoi_treemap`](@ref). The nodes of the hierarchy are numbered in breadth-first
order, so that the nodes at the top of the hierarchy come first, in their given order.

## Fields
- `labels::Vector{L}`: The labels of the nodes.
- `values::Vector{T}`: The values of the nodes; for a group, this is the sum of the values of its children.
- `parents::Vector{Int}`: The parent of each node, or `0` for the nodes at the top of the hierarchy.
- `children::Vector{Vector{Int}}`: The children of each node, which are empty for a leaf.
- `polygons::Vector{Vector{NTuple{2, T}}}`: The coordinates of the cell of each node, in counter-clockwise order with the first coordinate repeated at the end.
- `error::T`: The largest error of the Voronoi maps that were used to divide the cells, as for a [`VoronoiMap`](@ref).
- `clip_polygon::P`: The polygon that was divided.
"""
struct VoronoiTreemap{L, T, P}
    labels::Vector{L}
    values::Vector{T}
    parents::Vector{Int}
    children::Vector{Vector{Int}}
    polygons::Vector{Vector{NTuple{2, T}}}
    error::T
    clip_polygon::P
end

DelaunayTriangulation.num_polygons(treemap::VoronoiTreemap) = length(treemap.polygons)
DelaunayTriangulation.each_polygon_index(treemap::VoronoiTreemap) = eachindex(treemap.polygons)
DelaunayTriangulation.get_polygon_coordinates(treemap::VoronoiTreemap, i) = treemap.polygons[i]
DelaunayTriangulation.get_parent(treemap::VoronoiTreemap, i) = treemap.parents[i]
DelaunayTriangulation.get_children(treemap::VoronoiTreemap, i) = treemap.children[i]
function DelaunayTriangulation.get_area(treemap::VoronoiTreemap, i)
    coordinates = get_polygon_coordinates(treemap, i)
    isempty(coordinates) && return zero(eltype(treemap.values))
    A, _ = polygon_features(polygon_from_coordinates(coordinates)...)
    return A
end
get_label(treemap::VoronoiTreemap, i) = treemap.labels[i]
get_value(treemap::VoronoiTreemap, i) = treemap.values[i]
get_error(treemap::VoronoiTreemap) = treemap.error
is_leaf(treemap::VoronoiTreemap, i) = isempty(get_children(treemap, i))
each_leaf(treemap::VoronoiTreemap) = Iterators.filter(i -> is_leaf(treemap, i), each_polygon_index(treemap))
function get_depth(treemap::VoronoiTreemap, i) # 1 for the nodes at the top of the hierarchy
    depth = 1
    while get_parent(treemap, i) ≠ 0
        i = get_parent(treemap, i)
        depth += 1
    end
    return depth
end
function get_top_level_node(treemap::VoronoiTreemap, i) # the node at the top of the hierarchy that i is below
    while get_parent(treemap, i) ≠ 0
        i = get_parent(treemap, i)
    end
    return i
end

# The treemap itself is computed level by level, with a queue of the groups whose cells are still to be divided, rather
# than recursively.

"""
    voronoi_treemap(hierarchy; clip_polygon = unit_square(), rng = Random.default_rng(), predicates = AdaptiveKernel(), tol = 1.0e-3, maxiters = 500) -> VoronoiTreemap
    voronoi_treemap(labels, parents, values; root = "", kwargs...) -> VoronoiTreemap

Computes a Voronoi treemap of `hierarchy`, by dividing `clip_polygon` into a cell for each node at the top of the
hierarchy, and then the cell of each group into a cell for each of its children, using [`voronoi_map`](@ref). The second
method takes the hierarchy as a table; see [`hierarchy_from_table`](@ref).

## Arguments
- `hierarchy`: A vector of nodes, where a node is either a leaf, `label => value`, whose value is a positive number, or a group, `label => children`, whose children are a non-empty vector of nodes. The labels can be anything.

## Keyword Arguments
- `clip_polygon`, `rng`, `predicates`, `tol`, `maxiters`: As for [`voronoi_map`](@ref), which uses the same `tol` and `maxiters` for every group.

## Output
- `treemap`: The [`VoronoiTreemap`](@ref).
"""
function voronoi_treemap(
        hierarchy; clip_polygon = unit_square(), rng::AbstractRNG = Random.default_rng(),
        predicates::AbstractPredicateKernel = AdaptiveKernel(), tol = 1.0e-3, maxiters = 500,
    )
    check_hierarchy(hierarchy)
    T = number_type(first(clip_polygon))
    labels, values, parents, polygons = [], T[], Int[], Vector{NTuple{2, T}}[]
    error = zero(T)
    queue = Tuple{Any, Any, Int}[(hierarchy, clip_polygon, 0)] # the nodes to divide a polygon between, the polygon, and the parent of the nodes
    while !isempty(queue)
        nodes, polygon, parent = popfirst!(queue)
        vmap = voronoi_map(node_value.(nodes); clip_polygon = polygon, rng, predicates, tol, maxiters)
        error = max(error, get_error(vmap))
        for (k, node) in enumerate(nodes)
            push!(labels, first(node))
            push!(values, node_value(node))
            push!(parents, parent)
            push!(polygons, get_polygon_coordinates(vmap, k))
            is_group(node) && push!(queue, (last(node), polygon_from_coordinates(polygons[end]), length(labels)))
        end
    end
    children = [Int[] for _ in parents]
    for (i, parent) in enumerate(parents)
        parent ≠ 0 && push!(children[parent], i)
    end
    return VoronoiTreemap(identity.(labels), values, parents, children, polygons, error, clip_polygon)
end
voronoi_treemap(labels, parents, values; root = "", kwargs...) = voronoi_treemap(hierarchy_from_table(labels, parents, values; root); kwargs...)

# ### A Makie Recipe
#
# Finally, a plotting function. In Makie, this would be written next to `voronoiplot`, which also accepts either the
# data or a precomputed `VoronoiTessellation`. A [Makie recipe](https://docs.makie.org/stable/explanations/recipes)
# declares the attributes of the plot with their defaults and documentation, converts the arguments with
# `convert_arguments`, and builds the plot in `plot!` from Makie's compute graph, where each `map!` computes new
# attributes that are updated whenever their inputs change. Here, the treemap is computed from the data with
# `voronoi_treemap` unless it is given directly, the leaves are filled with the colours of their groups at the top of
# the hierarchy, the groups get thicker borders, and the leaves are labelled at their
# [poles of inaccessibility](pole_of_inaccessibility.md), the points inside of them furthest from their borders. The region to divide can be given as a `Rect2`, a `Circle`, or a vector of points.

"""
    voronoitreemap(hierarchy; kwargs...)
    voronoitreemap(labels, parents, values; kwargs...)
    voronoitreemap(treemap::VoronoiTreemap; kwargs...)

Plots a Voronoi treemap of hierarchical data, given as nested pairs (see `voronoi_treemap`), as a table of labels,
parents and values (see `hierarchy_from_table`), or as a `VoronoiTreemap` that has already been computed.
"""
@recipe Voronoitreemap (treemap,) begin
    "The region to divide, if the treemap is computed: a `Rect2`, a `Circle`, or a vector of the vertices of a convex polygon in counter-clockwise order."
    clip = Rect2d(0, 0, 1, 1)
    "The number of vertices of the polygon that approximates a `Circle` `clip`."
    clip_resolution = 128
    "The seed of the random number generator used for computing the treemap."
    seed = 1
    "The tolerance for the areas of the cells, as for `voronoi_treemap`."
    tol = 1.0e-3
    "The maximum number of iterations for dividing each cell, as for `voronoi_treemap`."
    maxiters = 500
    "The colors of the groups at the top of the hierarchy, as a vector that is reused if there are more groups than colors, or a single color."
    color = Makie.wong_colors(0.6)
    "Sets the color of the borders."
    strokecolor = :white
    "Sets the width of the borders of the leaves."
    strokewidth = 1.0
    "Sets the width of the borders of the groups at the top of the hierarchy. A group at depth `d` gets `group_strokewidth / d`, but at least `strokewidth`."
    group_strokewidth = 4.0
    "Determines whether to label the leaves."
    show_labels = true
    "A function `(label, value) -> text` that gives the text of the label of a leaf."
    label_formatter = (label, value) -> string(label)
    "Sets the size of the labels, as a number or as a function of the value of a leaf."
    label_size = @inherit fontsize
    "Sets the color of the labels."
    label_color = @inherit textcolor
    "Sets the font of the labels."
    label_font = @inherit font
    Makie.mixin_generic_plot_attributes()...
end

Makie.preferred_axis_type(::Voronoitreemap) = Axis

function Makie.convert_arguments(::Type{<:Voronoitreemap}, labels::AbstractVector, parents::AbstractVector, values::AbstractVector)
    return (hierarchy_from_table(labels, parents, values),)
end

"""
    treemap_clip_polygon(clip, resolution) -> (points, boundary_nodes)

Converts the `clip` attribute of [`voronoitreemap`](@ref) into a polygon for [`voronoi_treemap`](@ref), approximating a
`Circle` by a polygon with `resolution` vertices.
"""
function treemap_clip_polygon(rect::Rect2, _)
    (x, y), (w, h) = Float64.(rect.origin), Float64.(rect.widths)
    return [(x, y), (x + w, y), (x + w, y + h), (x, y + h)], [1, 2, 3, 4, 1]
end
function treemap_clip_polygon(circle::Circle, resolution)
    (cx, cy), r = Float64.(circle.center), Float64(circle.r)
    points = [(cx + r * cos(θ), cy + r * sin(θ)) for θ in LinRange(0, 2π, resolution + 1)[1:(end - 1)]]
    return points, [1:resolution; 1]
end
function treemap_clip_polygon(vertices::AbstractVector, _)
    points = [(Float64(p[1]), Float64(p[2])) for p in vertices]
    first(points) == last(points) && pop!(points)
    return points, [eachindex(points); 1]
end

function Makie.plot!(p::Voronoitreemap)
    map!(p, [:treemap, :clip, :clip_resolution, :seed, :tol, :maxiters], :computed_treemap) do treemap, clip, resolution, seed, tol, maxiters
        treemap isa VoronoiTreemap && return treemap
        return voronoi_treemap(treemap; clip_polygon = treemap_clip_polygon(clip, resolution), rng = StableRNG(seed), tol, maxiters)
    end
    map!(p, [:computed_treemap, :color], [:leaf_polygons, :leaf_colors]) do treemap, color
        leaves = collect(each_leaf(treemap))
        polygons = [Point2f.(get_polygon_coordinates(treemap, i)) for i in leaves]
        colors = color isa AbstractVector ? [color[mod1(get_top_level_node(treemap, i), length(color))] for i in leaves] : fill(color, length(leaves))
        return polygons, colors
    end
    map!(p, [:computed_treemap, :strokewidth, :group_strokewidth], [:group_polygons, :group_strokewidths]) do treemap, strokewidth, group_strokewidth
        groups = [i for i in each_polygon_index(treemap) if !is_leaf(treemap, i)]
        polygons = [Point2f.(get_polygon_coordinates(treemap, i)) for i in groups]
        widths = [max(strokewidth, group_strokewidth / get_depth(treemap, i)) for i in groups]
        return polygons, widths
    end
    map!(p, [:computed_treemap, :label_formatter, :label_size], [:label_positions, :label_texts, :label_sizes]) do treemap, formatter, label_size
        leaves = collect(each_leaf(treemap))
        positions = [Point2f(pole_of_inaccessibility(polygon_from_coordinates(get_polygon_coordinates(treemap, i))...)) for i in leaves]
        texts = [formatter(get_label(treemap, i), get_value(treemap, i)) for i in leaves]
        sizes = [label_size isa Function ? label_size(get_value(treemap, i)) : label_size for i in leaves]
        return positions, texts, sizes
    end
    poly!(p, p.leaf_polygons; color = p.leaf_colors, strokecolor = p.strokecolor, strokewidth = p.strokewidth)
    poly!(p, p.group_polygons; color = :transparent, strokecolor = p.strokecolor, strokewidth = p.group_strokewidths)
    text!(
        p, p.label_positions; text = p.label_texts, fontsize = p.label_sizes, color = p.label_color, font = p.label_font,
        align = (:center, :center), visible = p.show_labels,
    )
    return p
end

# ## Examples
#
# ### Dividing a Circle
#
# Let's start with a single Voronoi map, dividing a circle, approximated by a polygon with 128 sides, into five cells.
# The table compares the target areas with the areas of the cells.
circle = ([(cos(θ), sin(θ)) for θ in LinRange(0, 2π, 129)[1:(end - 1)]], [1:128; 1])
vmap = voronoi_map([30.0, 20.0, 10.0, 5.0, 2.0]; clip_polygon = circle, rng = StableRNG(123))
[[get_target_area(vmap, i) for i in each_polygon_index(vmap)] [get_area(vmap, i) for i in each_polygon_index(vmap)]]

#-
(get_error(vmap), num_iterations(vmap))

# Here are the cells, with their sites.
fig, ax, sc = poly(get_polygon_coordinates(vmap, 1), color = (:blue, 0.5), strokewidth = 2, axis = (aspect = DataAspect(),))
for (i, color) in zip(2:5, (:red, :green, :orange, :purple))
    poly!(ax, get_polygon_coordinates(vmap, i), color = (color, 0.5), strokewidth = 2)
end
scatter!(ax, [get_site(vmap, i) for i in each_polygon_index(vmap)], color = :black)
fig
@test get_error(vmap) < 1.0e-3 #src
@test sum(get_area(vmap, i) for i in each_polygon_index(vmap)) ≈ polygon_features(circle...)[1] #src
@test all(i -> distance_to_polygon(get_site(vmap, i), polygon_from_coordinates(get_polygon_coordinates(vmap, i))...) > 0, each_polygon_index(vmap)) #src
@test_throws ArgumentError voronoi_map([1.0, -1.0]) #src
@test_throws ArgumentError voronoi_map(Float64[]) #src
@test_throws ArgumentError voronoi_map([1.0, 2.0]; clip_polygon = ([(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)], [1, 2, 3, 4, 1])) #src
@test_throws ArgumentError voronoi_map([1.0, 2.0]; clip_polygon = ([(0.0, 0.0), (2.0, 0.0), (1.0, 0.5), (2.0, 2.0), (0.0, 2.0)], [1, 2, 3, 4, 5, 1])) #src

# ### A Budget
#
# Here is a fictional budget, with five categories of items, one of which has a group of its own.
budget = [
    "Health" => ["Hospitals" => 30.0, "Primary care" => ["Doctors" => 12.0, "Dentists" => 4.0, "Pharmacies" => 2.0], "Medicines" => 12.0, "Public health" => 5.0],
    "Education" => ["Schools" => 25.0, "Universities" => 14.0, "Early childhood" => 6.0],
    "Transport" => ["Roads" => 12.0, "Rail" => 9.0, "Buses" => 4.0, "Cycling" => 1.5],
    "Housing" => ["Social housing" => 10.0, "Grants" => 3.0],
    "Environment" => ["Parks" => 3.0, "Waste" => 4.0, "Water" => 5.0, "Energy" => 6.0, "Research" => 2.0],
]
treemap = voronoi_treemap(budget; clip_polygon = circle, rng = StableRNG(123))
[(get_label(treemap, i), get_depth(treemap, i), get_value(treemap, i)) for i in each_polygon_index(treemap)]
@test sum(i -> get_area(treemap, i), each_leaf(treemap)) ≈ polygon_features(circle...)[1] #src
@test all(i -> isapprox(get_area(treemap, i), polygon_features(circle...)[1] * get_value(treemap, i) / sum(j -> get_value(treemap, j), 1:5); rtol = 0.02), 1:5) #src
@test get_error(treemap) < 1.0e-3 #src
@test_throws ArgumentError voronoi_treemap(["a" => 1.0, "b" => -2.0]) #src
@test_throws ArgumentError voronoi_treemap(["a" => 1.0, "b" => String[]]) #src

# With the recipe, the whole treemap is a single call.
fig, ax, plt = voronoitreemap(budget; clip = Circle(Point2f(0, 0), 1.0f0), label_size = 10, figure = (size = (700, 700),), axis = (aspect = DataAspect(),))
hidedecorations!(ax)
hidespines!(ax)
fig
@test_reference reference_path(fig_path, "voronoi_treemap_ex_1.png") fig #src

# ### Data from a Table
#
# Here is some data in the form of a table. The recipe converts it with `hierarchy_from_table`, and divides the unit
# square, which is the default region.
labels = ["Fruit", "Apples", "Bananas", "Cherries", "Vegetables", "Carrots", "Peas", "Leeks", "Grains", "Rice", "Oats"]
parents = ["", "Fruit", "Fruit", "Fruit", "", "Vegetables", "Vegetables", "Vegetables", "", "Grains", "Grains"]
amounts = [0, 40, 25, 10, 0, 30, 15, 12, 0, 35, 20] # the values of the groups are ignored
fig, ax, plt = voronoitreemap(labels, parents, amounts; label_size = 14, axis = (aspect = DataAspect(),))
hidedecorations!(ax)
hidespines!(ax)
fig
@test_reference reference_path(fig_path, "voronoi_treemap_ex_2.png") fig #src
@test hierarchy_from_table(labels, parents, amounts)[1] == ("Fruit" => ["Apples" => 40, "Bananas" => 25, "Cherries" => 10]) #src
@test_throws ArgumentError hierarchy_from_table(["a", "b"], ["", "c"], [1, 2]) #src

# ### The Global Economy by GDP
#
# The example [The Global Economy by GDP](https://bl.ocks.org/Kcnarf/fa95aa7b076f537c00aed614c29bb568) of
# d3-voronoi-treemap, a remake of [an article by HowMuch.net](https://howmuch.net/articles/the-global-economy-by-GDP),
# shows each country's share of the world's GDP, in percent, as of January 2017, grouped by continent. The labels of
# the leaves can be anything, so here we label each country with its name and its code, and the `label_formatter` shows
# the code instead of the name for the smaller countries, as in the original.
gdp = [
    "Asia" => [
        (name = "China", code = "CN") => 14.84, (name = "Japan", code = "JP") => 5.91, (name = "India", code = "IN") => 2.83,
        (name = "South Korea", code = "KR") => 1.86, (name = "Russia", code = "RU") => 1.8, (name = "Indonesia", code = "ID") => 1.16,
        (name = "Turkey", code = "TR") => 0.97, (name = "Saudi Arabia", code = "SA") => 0.87, (name = "Iran", code = "IR") => 0.57,
        (name = "Thailand", code = "TH") => 0.53, (name = "United Arab Emirates", code = "AE") => 0.5, (name = "Hong Kong", code = "HK") => 0.42,
        (name = "Israel", code = "IL") => 0.4, (name = "Malaysia", code = "MY") => 0.4, (name = "Singapore", code = "SG") => 0.39,
        (name = "Philippines", code = "PH") => 0.39,
    ],
    "North America" => [(name = "United States", code = "US") => 24.32, (name = "Canada", code = "CA") => 2.09, (name = "Mexico", code = "MX") => 1.54],
    "Europe" => [
        (name = "Germany", code = "DE") => 4.54, (name = "United Kingdom", code = "UK") => 3.85, (name = "France", code = "FR") => 3.26,
        (name = "Italy", code = "IT") => 2.46, (name = "Spain", code = "ES") => 1.62, (name = "Netherlands", code = "NL") => 1.01,
        (name = "Switzerland", code = "CH") => 0.9, (name = "Sweden", code = "SE") => 0.67, (name = "Poland", code = "PL") => 0.64,
        (name = "Belgium", code = "BE") => 0.61, (name = "Norway", code = "NO") => 0.52, (name = "Austria", code = "AT") => 0.51,
        (name = "Denmark", code = "DK") => 0.4, (name = "Ireland", code = "IE") => 0.38,
    ],
    "South America" => [(name = "Brazil", code = "BR") => 2.39, (name = "Argentina", code = "AR") => 0.79, (name = "Venezuela", code = "VE") => 0.5, (name = "Colombia", code = "CO") => 0.39],
    "Australia" => [(name = "Australia", code = "AU") => 1.81],
    "Africa" => [(name = "Nigeria", code = "NG") => 0.65, (name = "Egypt", code = "EG") => 0.45, (name = "South Africa", code = "ZA") => 0.42],
    "Rest of the World" => [(name = "Rest of the World", code = "RotW") => 9.41],
]
continent_colors = ["#f58321", "#ef1621", "#77bc45", "#4aaaea", "#00acad", "#f575a3", "#592c94"]
gdp_label(label, value) = string(value < 1 ? label.code : label.name, "\n", value, "%")
gdp_label_size(value) = clamp(8 + 12(value - 3) / 17, 8, 20) # from 8 for values up to 3 to 20 for values from 20

# The original divides a polygon with 60 sides, representing the Earth, and adds a title and a legend. This time, we
# compute the treemap first and then plot it, which lets us check its areas.
earth = ([(cos(θ), sin(θ)) for θ in LinRange(0, 2π, 61)[1:(end - 1)]], [1:60; 1])
treemap = voronoi_treemap(gdp; clip_polygon = earth, rng = StableRNG(1))
fig = Figure(size = (960, 500), backgroundcolor = RGBf(0.98, 0.98, 0.98))
ax = Axis(fig[1, 1], aspect = DataAspect(), title = "The Global Economy by GDP (as of 01/2017)", titlesize = 20, backgroundcolor = :transparent)
hidedecorations!(ax)
hidespines!(ax)
voronoitreemap!(
    ax, treemap; color = continent_colors, strokewidth = 1, group_strokewidth = 1,
    label_formatter = gdp_label, label_size = gdp_label_size, label_color = :white,
)
poly!(ax, Point2f.(first(earth)), color = :transparent, strokecolor = :lightgrey, strokewidth = 4)
Legend(
    fig[1, 1], [PolyElement(color = color) for color in continent_colors], first.(gdp), "Continents",
    halign = :left, valign = :bottom, tellwidth = false, tellheight = false, framevisible = false, backgroundcolor = :transparent,
)
fig
@test_reference reference_path(fig_path, "voronoi_treemap_ex_3.png") fig #src
@test all(i -> isapprox(get_area(treemap, i), polygon_features(earth...)[1] * get_value(treemap, i) / sum(j -> get_value(treemap, j), 1:7); rtol = 0.05), 1:7) #src
@test count(i -> is_leaf(treemap, i), each_polygon_index(treemap)) == 42 #src

# ### Other Shapes
#
# The README of d3-voronoi-treemap shows the same treemap in a square, a hexagon, a diamond and a circle.
clips = [
    Rect2d(-1, -1, 2, 2),
    [Point2f(cos(θ), sin(θ)) for θ in (0:5) .* (π / 3)],
    [Point2f(0, -1), Point2f(1, 0), Point2f(0, 1), Point2f(-1, 0)],
    Circle(Point2f(0, 0), 1.0f0),
]
fig = Figure(size = (800, 800))
for (k, clip) in enumerate(clips)
    shape_ax = Axis(fig[fldmod1(k, 2)...], aspect = DataAspect())
    hidedecorations!(shape_ax)
    hidespines!(shape_ax)
    voronoitreemap!(
        shape_ax, gdp; clip, clip_resolution = 60, color = continent_colors, strokewidth = 1, group_strokewidth = 1,
        label_formatter = gdp_label, label_size = value -> gdp_label_size(value) - 2, label_color = :white,
    )
end
fig
@test_reference reference_path(fig_path, "voronoi_treemap_ex_4.png") fig #src
