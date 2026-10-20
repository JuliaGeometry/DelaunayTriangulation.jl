# Weighted triangulations computed independently in MATLAB, and slow reference computations for weighted triangulations.

export get_weighted_example, get_nearest_power_point, NUM_WEGT, NUM_CWEGT, get_convex_polygon_weighted_example

#=
These weighted Delaunay examples come from the following MATLAB code, which exploits 
the relationship between the convex hull of the lifting map and the Delaunay triangulation.
1:78 are random triangulations, and 79:155 are triangulations of convex polygons.

For computing these triangulations, we first define the MATLAB functions:

    ```matlab
    function [cells, lifted] = LiftingMap(points,weights)
    num = size(points, 1);
    lifted = zeros(num, 3);
    for i = 1:num
        p = points(i, :);
        lifted(i, :) = [p sum(p.^2) - weights(i)];
    end
    ch = convhull(lifted);
    sep_ch = reshape(lifted(ch, :), [size(ch) 3]);
    vecs = diff(sep_ch, 1, 2);
    n = cross(vecs(:, 1, :), vecs(:, 2, :), 3);
    downward_faces = n(:, :, 3) <= 0;
    cells = ch(downward_faces, :);
    end
    function save_to_file(points, weights, n)
    [lift_cells, ~] = LiftingMap(points,weights);
    res = [(1+size(points, 1)) * ones(1, 3); [points weights]; lift_cells];
    % so the points are in res(2:res(1), 1:2), the weights are in res(2:res(1), 3),
    % and the triangles are in res((res(1)+1):end, :).
    writematrix(res, fullfile(pwd, 'tempmats', strcat(string(n), '.txt')));
    end
    ```
Then, for the random triangulations, we use:
    ```matlab
    digits(16)
    rng(123)
    ctr = 1;
    for pass = 1:3
        if pass == 1
            range = 4:20;
        elseif pass == 2
            range = 21:25:1500;
        else 
            range = 15000;
        end
        for n = range
            points = randn(n, 2);
            weights = randn(n, 1);
            s = 10*rand;
            u = rand;
            if u < 1/250
                weights = 0 * weights;
            elseif (1/250 <= u) && (u < 1/5)
                weights = weights.^2;
            elseif (1/5 <= u) && (u < 0.6)
                weights = s * weights;
            elseif (u <= 0.6) && (u < 0.65)
                sparse_prob = rand(n);
                weights(sparse_prob >= 0.95) = 0;
            end
            save_to_file(points, weights, ctr)
            ctr = ctr + 1;
        end
    end
    ```
The triangulations of convex polygons are computed with the code:
    ```matlab
    digits(16)
    rng(123)
    ctr = 79;
    for pass = 1:3
        if pass == 1
            range = 6:20;
        elseif pass == 2
            range = 21:25:1500;
        else
            range = [100, 15000];
        end
        for n = range
            if pass < 3
                n0 = 0;
                while n0 <= 3
                    points = randn(n, 2);
                    convex_poly = points(convhull(points), :);
                    points = convex_poly(1:end-1, :);
                    n0 = size(points, 1);
                end
                weights = randn(n0, 1);
            else 
                theta = linspace(0, 2*pi, n);
                x = cos(theta);
                y = sin(theta);
                points = [x(:), y(:)];
                points = points(1:end-1, :);
                n0 = size(points, 1);
                weights = randn(n0, 1);
            end
            s = 10*rand;
            u = rand;
            if u < 1/250
                weights = 0 * weights;
            elseif (1/250 <= u) && (u < 1/5)
                weights = weights.^2;
            elseif (1/5 <= u) && (u < 0.6)
                weights = s * weights;
            elseif (u <= 0.6) && (u < 0.65)
                sparse_prob = rand(n0);
                weights(sparse_prob >= 0.95) = 0;
            end
            save_to_file(points, weights, ctr)
            ctr = ctr + 1;
        end
    end
    ```
=#
const WGT_DIR = joinpath(dirname(dirname(dirname(pathof(DelaunayTriangulation)))), "assets", "test", "triangulation", "weighted_delaunay_mats")

const NUM_WEGT = length(readdir(WGT_DIR))

const NUM_CWEGT = NUM_WEGT - 78 # number of convex polygon weighted examples

function get_weighted_example(i)
    # get the points, weights, and triangles
    file = "$i.txt"
    mat = readdlm(joinpath(WGT_DIR, file), ',')
    num_pts = Int(mat[1, 1])
    points = mat[2:num_pts, 1:2]
    weights = mat[2:num_pts, 3]
    triangles = Int.(mat[(num_pts+1):end, :])
    triangles = Set([Tuple(tri) for tri in eachrow(triangles)])
    _triangles = empty(triangles)
    for T in triangles
        i, j, k = T
        p, q, r = eachrow(points[[i, j, k], :])
        or = DT.triangle_orientation(p, q, r)
        if !DT.is_positively_oriented(or)
            push!(_triangles, (j, i, k))
        else
            push!(_triangles, (i, j, k))
        end
    end
    triangles = _triangles
    points = ElasticMatrix(points')

    # since some of the vertices might be submerged, we can't just get the 
    # convex hull of points to find the boundary (in cases of collinear points). instead, let's find 
    # all the edges which have only one adjoining triangle, and then sort them.
    d = Dict{NTuple{2,Int},Int}()
    for T in triangles
        i, j, k = T
        d[(i, j)] = k
        d[(j, k)] = i
        d[(k, i)] = j
    end
    boundary_nodes = Set{Int}()
    for (i, j) in keys(d)
        if !haskey(d, (j, i))
            push!(boundary_nodes, i, j)
        end
    end
    boundary_nodes = collect(boundary_nodes)
    DT.sort_convex_polygon!(boundary_nodes, points)
    push!(boundary_nodes, boundary_nodes[1])

    # now get all the submerged vertices 
    submerged_vertices = Set{Int}()
    nonsubmerged_vertices = Set{Int}()
    for T in triangles
        push!(nonsubmerged_vertices, T...)
    end
    for i in axes(points, 2)
        if i ∉ nonsubmerged_vertices
            push!(submerged_vertices, i)
        end
    end

    # return 
    tri = Triangulation(points, triangles, boundary_nodes; weights)
    unlock_convex_hull!(tri)
    return (tri=tri, submerged_vertices=submerged_vertices, nonsubmerged_vertices=nonsubmerged_vertices, weights=weights)
end

get_unstructured_weighted_example(i) = get_weighted_example(i)

function get_convex_polygon_weighted_example(i)
    tri, submerged_vertices, nonsubmerged_vertices, weights = get_weighted_example(78 + i)
    S = eachindex(weights)
    representative_point_list = DT.get_representative_point_list(tri)
    cx, cy = DT.mean_points(get_points(tri), S)
    representative_point_list[1] = DT.RepresentativeCoordinates(cx, cy, length(S))
    return (tri=tri, submerged_vertices=submerged_vertices, nonsubmerged_vertices=nonsubmerged_vertices, weights=weights, S=S)
end

function get_all_distances_to_witness_planes(tri, i)
    distances = Dict{NTuple{3,Int},Float64}()
    for T in each_solid_triangle(tri)
        T = DT.sort_triangle(T)
        δ = DT.get_distance_to_witness_plane(tri, i, T)
        distances[T] = δ
    end
    return distances
end

function get_nearest_power_point(tri, i)
    all_dists = Dict{Int,Float64}()
    for j in each_solid_vertex(tri)
        p = DT.get_point(tri, i)
        q = DT.get_point(tri, j)
        all_dists[j] = norm(p .- q)^2 - get_weight(tri, i) - get_weight(tri, j)
    end
    return argmin(all_dists)
end
