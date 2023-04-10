```@meta
CurrentModule = DelaunayTriangulation
```

# Pole of Inaccessibility and Polygons 

We provide a function for computing the pole of inaccessibility of a given polygon, namely the point inside the polygon that is furthest from the boundary. Our method is primarily based on [this blogpost](https://blog.mapbox.com/a-new-algorithm-for-finding-a-visual-center-of-a-polygon-7c77e6492fbc), recursively subdividing the polygon using quadtree partitioning. The function for this is `pole_of_inaccessibility`:

```@docs 
pole_of_inaccessibility
```

We needed this method since the point we need to associate ghost vertices with must be inside the domain, and so other representative points like centroids or arithmetic averages would not be sufficient if the domain is non-convex.

Below we also list some other relevant docstrings.

```@docs 
Cell 
CellQueue 
polygon_features 
squared_distance_to_segment 
distance_to_polygon 
polygon_bounds 
```

If you need to compute this for multiple boundaries, meaning multiple poles, use `compute_representative_points!`.

```@docs 
compute_representative_points!
```

## Example 

Below is a simple example of computing this pole of inaccessibility.

```julia 
pts = [0.0 8.0
      2.0 5.0
      3.0 7.0
      1.81907 8.13422
      3.22963 8.865
      4.24931 7.74335
      4.50423 5.87393
      3.67149 4.3784
      2.73678 2.62795
      5.50691 1.38734
      8.43 2.74691
      9.7046 5.53404
      8.56595 7.79433
      6.71353 9.03494
      4.13034 9.66375
      2.75378 10.3775
      1.0883 10.4965
      -1.138 9.83369
      -2.25965 8.45712
      -2.78649 5.94191
      -1.39292 3.64763
      0.323538 4.97322
      -0.900078 6.6217
      0.98633 9.68074
      0.153591 9.54478
      0.272554 8.66106
      2.90673 8.18521
      2.12497 9.42582
      7.27436 2.7979
      3.0 4.0
      5.33697 1.88019]'
boundary_nodes = [
      [[1, 4, 3, 2], [2, 9, 10, 11, 8, 7, 12], [12, 6, 13, 5, 14, 15, 16, 17, 16], [16, 17, 18, 19, 20, 21, 22, 23, 1]],
      [[26, 25, 24], [24, 28, 27, 26]],
      [[29, 30, 31, 29]]
]
x, y = DT.pole_of_inaccessibility(pts, boundary_nodes)

fig = Figure()
ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y")
bn1 = pts[:, unique(reduce(vcat, boundary_nodes[1]))] |> x -> hcat(x, x[:, begin])
bn2 = pts[:, unique(reduce(vcat, boundary_nodes[2]))] |> x -> hcat(x, x[:, begin])
bn3 = pts[:, unique(reduce(vcat, boundary_nodes[3]))] |> x -> hcat(x, x[:, begin])
lines!(ax, bn1, color=:red, linewidth=4)
lines!(ax, bn2, color=:red, linewidth=4)
lines!(ax, bn3, color=:red, linewidth=4)
scatter!(ax, [x], [y], color=:blue, markersize=23)
```

```@raw html
<figure>
    <img src='../figs/pole_of_inaccessibility.png', alt='Pole of inaccessibility'><br>
</figure>
```