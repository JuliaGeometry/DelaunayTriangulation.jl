function quadrant_sort!(sorter::HilbertSorter, xsplit, ysplit)
    points = sorter.points
    directions = sorter.directions
    op1 = HilbertOp(directions.xdir, directions.coord, xsplit)
    op2 = HilbertOp(directions.ydir, next(directions.coord), ysplit)
    op3 = HilbertOp(reverse(directions.ydir), next(directions.coord), ysplit)
    i₁, i₅ = sorter.first, sorter.last
    i₃ = separate!(view(points, i₁:i₅), op1; carry=view(sorter.perm, i₁:i₅)) + i₁
    i₂ = separate!(view(points, i₁:i₃-1), op2; carry=view(sorter.perm, i₁:i₃-1)) + i₁
    i₄ = separate!(view(points, i₃:i₅), op3; carry=view(sorter.perm, i₃:i₅)) + i₃
    return i₁, i₂, i₃, i₄, i₅
end

function to_quadrant1(sorter::HilbertSorter, xsplit, ysplit, i₁, i₂, _, _, _)
    xmin, _, ymin, _ = sorter.bbox
    rotated_bbox = (ymin, ysplit, xmin, xsplit) # scaled, rotated, and reflected
    directions = sorter.directions
    new_directions = HilbertDirections(next(directions.coord), directions.ydir, directions.xdir)
    return HilbertSorter(sorter.points, sorter.perm, rotated_bbox, i₁, i₂ - 1, new_directions, sorter.splitter, sorter.capacity, sorter.rng)
end
function to_quadrant2(sorter::HilbertSorter, xsplit, ysplit, _, i₂, i₃, _, _)
    xmin, _, _, ymax = sorter.bbox
    rotated_bbox = (xmin, xsplit, ysplit, ymax) # scaled
    directions = sorter.directions
    new_directions = HilbertDirections(directions.coord, directions.ydir, directions.xdir)
    return HilbertSorter(sorter.points, sorter.perm, rotated_bbox, i₂, i₃ - 1, new_directions, sorter.splitter, sorter.capacity, sorter.rng)
end
function to_quadrant3(sorter::HilbertSorter, xsplit, ysplit, _, _, i₃, i₄, _)
    _, xmax, _, ymax = sorter.bbox
    rotated_bbox = (xsplit, xmax, ysplit, ymax) # scaled
    return HilbertSorter(sorter.points, sorter.perm, rotated_bbox, i₃, i₄ - 1, sorter.directions, sorter.splitter, sorter.capacity, sorter.rng)
end
function to_quadrant4(sorter::HilbertSorter, xsplit, ysplit, _, _, _, i₄, i₅)
    _, xmax, ymin, _ = sorter.bbox
    rotated_bbox = (ysplit, ymin, xmax, xsplit) # scaled, rotated, and reflected
    directions = sorter.directions
    new_directions = HilbertDirections(next(directions.coord), reverse(directions.ydir), reverse(directions.xdir))
    return HilbertSorter(sorter.points, sorter.perm, rotated_bbox, i₄, i₅, new_directions, sorter.splitter, sorter.capacity, sorter.rng)
end

function hilbert_sort(points; capacity=8, splitter=Middle(), rng=Random.default_rng())
    if splitter === Median()
        throw(ArgumentError("Median splitting is not yet implemented."))
        #=
        The reason it's not yet implemented is because 
            E = (12.8728090546081, 1.057828011921)
            F = (2.9914997728536, 4.9162440171775)
            G = (4.9677616292045, 10.8920834399529)
            H = (2.7091766505178, 21.1498235514886)
            I = (11.0847626131478, 18.9853462802471)
            J = (17.0263810513976, 26.8690055903223)
            K = (23.2289344966548, 24.9012989801027)
            L = (22.9722771127131, 22.7197112165985)
            M = (25.0683124149035, 22.9763686005402)
            N = (25.1110886455604, 16.6027102326552)
            O = (24.9399837229326, 15.1055421596621)
            P = (31.014208476219, 9.3307510209744)
            Q = (25.0683124149035, 6.8497296428715)
            R = (28.6938091798674, 0.903833581556)
            S = (22.9722771127131, 0.7755048895852)
            T = (30.6719986309634, 31.1038524253599)
            points = [E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T]
            spoints, perm, order = DelaunayTriangulation.hilbert_sort(points; capacity=1, splitter=DelaunayTriangulation.Median())
        leads to a stack overflow for some reason. Maybe I am not understanding the actual definition of median splitting.
        =#
    end
    sorter = HilbertSorter(points; capacity, splitter, rng)
    order = hilbert_sort!(sorter)
    return unwrap(sorter.points), sorter.perm, order
end

function hilbert_sort!(sorter::HilbertSorter, depth=0)
    stop_sorter(sorter) && return depth
    depth += 1
    xsplit, ysplit = get_split(sorter)
    i₁, i₂, i₃, i₄, i₅ = quadrant_sort!(sorter, xsplit, ysplit)
    depth1 = hilbert_sort!(to_quadrant1(sorter, xsplit, ysplit, i₁, i₂, i₃, i₄, i₅), depth)
    depth2 = hilbert_sort!(to_quadrant2(sorter, xsplit, ysplit, i₁, i₂, i₃, i₄, i₅), depth)
    depth3 = hilbert_sort!(to_quadrant3(sorter, xsplit, ysplit, i₁, i₂, i₃, i₄, i₅), depth)
    depth4 = hilbert_sort!(to_quadrant4(sorter, xsplit, ysplit, i₁, i₂, i₃, i₄, i₅), depth)
    return max(depth1, depth2, depth3, depth4)
end
