x = [2.0, -4.0, -5.0, 5.0, 10.0]
y = [3.0, 3.0, -5.0, -5.0, 2.0]
const GMSH_PATH = "./gmsh-4.9.4-Windows64/gmsh.exe"
for r in 0.05:10.0
   ( T, adj, adj2v, DG, nodes), BN = DT.generate_mesh(x, y, 0.2; gmsh_path=GMSH_PATH)

    function DT._get_point(pts::AbstractMatrix, i)
        return @view pts[:, i]
    end
    function DT._eachindex(pts::AbstractMatrix)
        return axes(pts, 2)
    end

    _T, _adj, _adj2v, _DG = DT.triangulate_bowyer(nodes; trim=true)
    @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)
end
