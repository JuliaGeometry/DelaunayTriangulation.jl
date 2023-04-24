"""
    voronoi(!)(vor::VoronoiTessellation)

Plot a Voronoi tessellation.

# Attributes 
- `markersize=11`

Size of the generators.
- `show_generators=true`

Whether to show the generators.
- `generator_color=:black`

Colour of the generators.
- `strokecolor=:black`

Colour of the edges.
- `polygon_color=(:white, 0)`

Colour of the polygons.
- `strokewidth=1`

Width of the edges.
- `unbounded_edge_extension_factor=2.0`

Factor that scales the bounding box for clipping unbounded edges.
- `colormap=th.colormap`

Colormap for the polygons.
- `colorrange=get(th.attributes, :colorrange, MakieCore.automatic)`

Range of the colormap.
- `cycle=[:color]`

Cycle of the colormap.
"""
MakieCore.@recipe(Voronoiplot, vorn) do scene
    th = MakieCore.default_theme(scene, Mesh)
    return MakieCore.Attributes(;
        markersize=11,
        show_generators=true,
        generator_color=:black,
        strokecolor=:black,
        polygon_color=(:white, 0),
        strokewidth=1,
        unbounded_edge_extension_factor=2.0,
        colormap=th.colormap,
        colorrange=get(th.attributes, :colorrange, MakieCore.automatic),
        cycle=[:color]
    )
end

function MakieCore.plot!(p::Voronoiplot)
    ## Extract 
    vorn = p[:vorn]
    markersize = p[:markersize]
    show_generators = p[:show_generators]
    generator_color = p[:generator_color]
    strokecolor = p[:strokecolor]
    polygon_color = p[:polygon_color]
    strokewidth = p[:strokewidth]
    unbounded_edge_extension_factor = p[:unbounded_edge_extension_factor]
    colormap = p[:colormap]
    colorrange = p[:colorrange]
    cycle = p[:cycle]
    bbox = polygon_bounds(vorn[], unbounded_edge_extension_factor[])
    xmin, xmax, ymin, ymax = bbox
    bbox = [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)]
    bbox_order = [1, 2, 3, 4, 1]

    if !(typeof(polygon_color[]) <: AbstractVector)
        n = num_polygons(vorn[])
        polygon_color[] = [polygon_color[] for _ in 1:n]
    end

    ## Define all the necessary observables 
    generators_2f = MakieCore.Observable(NTuple{2,Float64}[])
    polygons = MakieCore.Observable([NTuple{2,Float64}[] for _ in each_polygon(vorn[])])

    ## Define the plotting function 
    function update_plot(vorn)
        empty!(generators_2f[])
        empty!(polygons[])
        resize!(generators_2f[], num_generators(vorn))
        resize!(polygons[], num_generators(vorn))
        for i in each_generator(vorn)
            generators_2f[][i] = get_generator(vorn, i)
            polygons[][i] = get_polygon_coordinates(vorn, i, bbox, bbox_order)
        end
    end

    ## Connect the plot so that it updates whenever we change a value 
    MakieCore.Observables.onany(vorn)

    ## Call it once to prepopulate with current values 
    update_plot(vorn[])

    ## Now plot 
    for i in eachindex(polygons[])
        poly!(p, polygons[][i], color=polygon_color[][i], strokecolor=strokecolor[],
            strokewidth=strokewidth[], colormap=colormap[],
            colorrange=colorrange[], cycle=cycle[])
    end
    if show_generators[]
        scatter!(p, generators_2f[], markersize=markersize[], color=generator_color[])
    end
    return p
end

function get_polygon_colors(vorn::VoronoiTessellation, cmap)
    F = number_type(vorn)
    gtr = [get_generator(vorn, i) for i in each_generator(vorn)]
    gtr_mat = reinterpret(reshape, F, gtr)
    colors = get(cmap, gtr_mat, :extrema)
    return [(a + b) / 2 for (a, b) in eachcol(colors)]
end