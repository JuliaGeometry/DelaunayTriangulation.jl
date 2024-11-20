struct LiangBarskyClipper{T}
    p::T 
    q::T
    t1::T 
    t2::T
end 

@inline get_p(clipper::LiangBarskyClipper) = clipper.p 
@inline get_q(clipper::LiangBarskyClipper) = clipper.q
@inline get_t1(clipper::LiangBarskyClipper) = clipper.t1 
@inline get_t2(clipper::LiangBarskyClipper) = clipper.t2