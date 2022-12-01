include("functions.jl")

function prepare_lattice(nx, ny, Δx, Δy)
    cells = zeros(SVector{2,Float64}, nx * ny)
    cell_ij_idx = CartesianIndices((1:nx, 1:ny))
    cell_k_idx = LinearIndices(cell_ij_idx)
    for (ij, k) in zip(cell_ij_idx, cell_k_idx)
        i, j = Tuple(ij)
        # Compute the cells
        cell_x = (i - 1) * Δx
        cell_y = (j - 1) * Δy
        cells[k] = @SVector[cell_x, cell_y]
    end
    return cells
end

nx = 6
ny = 6
Δx = 1.0
Δy = 1.0
pts = prepare_lattice(nx, ny, Δx, Δy)

T, adj, adj2v, DG = DT.lazy_triangulate_bowyer(pts)
fig = Figure(fontsize=55)
ax = Axis(fig[1, 1], aspect=1)
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0.0), plot_ghost_edges=true, DG=DG)
hidedecorations!(ax)
xlims!(minimum(getx.(pts))-1, maximum(getx.(pts))+1)
ylims!(minimum(gety.(pts))-1, maximum(gety.(pts))+1)

fig