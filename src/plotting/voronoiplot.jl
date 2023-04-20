MakieCore.@recipe(Voronoiplot, vorn) do scene
    th = MakieCore.default_theme(scene, Mesh)
    return MakieCore.Attributes(;
        markersize=11,
        show_generators=true,
        show_cell_points=false,
        generator_color=:black,
        cell_point_color=:red,
        strokecolor=:black,
        cell_color=(:white, 0),
        strokewidth=1,
        show_obstacles=true,
        obstacle_color=:magenta,
        obstacle_linewidth=2,
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
    show_cell_points = p[:show_cell_points]
    generator_color = p[:generator_color]
    cell_point_color = p[:cell_point_color]
    strokecolor = p[:strokecolor]
    cell_color = p[:cell_color]
    strokewidth = p[:strokewidth]
    show_obstacles = p[:show_obstacles]
    obstacle_color = p[:obstacle_color]
    obstacle_linewidth = p[:obstacle_linewidth]
    unbounded_edge_extension_factor = p[:unbounded_edge_extension_factor]
    colormap = p[:colormap]
    colorrange = p[:colorrange]
    cycle = p[:cycle]
    bbox = polygon_bounds(vorn[], unbounded_edge_extension_factor[])
    xmin, xmax, ymin, ymax = bbox
    bbox = [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)]
    bbox_order = [1, 2, 3, 4, 1]

    if !(typeof(cell_color[]) <: AbstractVector)
        n = num_cells(vorn[])
        cell_color[] = [cell_color[] for _ in 1:n]
    end

    ## Define all the necessary observables 
    generators_2f = MakieCore.Observable(NTuple{2,Float64}[])
    cells = MakieCore.Observable([NTuple{2,Float64}[] for _ in each_cell(vorn[])])
    obstacles_2f = MakieCore.Observable(NTuple{2,Float64}[])
    cell_points_2f = MakieCore.Observable(NTuple{2,Float64}[])

    ## Define the plotting function 
    function update_plot(vorn)
        empty!(generators_2f[])
        empty!(cells[])
        empty!(obstacles_2f[])

        for i in each_cell(vorn)
            push!(generators_2f[], get_generator(vorn, i))
            push!(cells[], get_cell_coordinates(vorn, i, bbox, bbox_order))
        end

        for e in each_obstacle(vorn)
            u, v = edge_indices(e)
            _p, q = get_cell_point(vorn, u, v)
            push!(obstacles_2f[], _p, q)
        end

        for i in each_cell_point(vorn)
            push!(cell_points_2f[], get_cell_point(vorn, i))
        end
    end

    ## Connect the plot so that it updates whenever we change a value 
    MakieCore.Observables.onany(vorn)

    ## Call it once to prepopulate with current values 
    update_plot(vorn[])

    ## Now plot 
    for i in eachindex(cells[])
        poly!(p, cells[][i], color=cell_color[][i], strokecolor=strokecolor[],
            strokewidth=strokewidth[], colormap=colormap[],
            colorrange=colorrange[], cycle=cycle[])
    end
    if show_generators[]
        scatter!(p, generators_2f[], markersize=markersize[], color=generator_color[])
    end
    if show_cell_points[]
        scatter!(p, cell_points_2f[], markersize=markersize[], color=cell_point_color[])
    end
    if show_obstacles[]
        linesegments!(p, obstacles_2f[], color=obstacle_color[], linewidth=obstacle_linewidth[])
    end
    return p
end