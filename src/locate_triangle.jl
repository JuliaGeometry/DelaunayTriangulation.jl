function locate_triangle(HG::HistoryGraph, pts, r::I, init=find_root(HG; method=:rng)) where {I}
    if out_deg(HG, init) == 0
        return init, isintriangle(init, pts, r)
    end
    out = out_neighbors(HG, init)
    for T in out
        isin = isintriangle(T, pts, r)
        (isin == I(1) || isin == I(0)) && return locate_triangle(HG, pts, r, T)
    end
end