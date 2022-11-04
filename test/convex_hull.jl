using Random 
Random.seed!(2939818188)
p1 = [-3.0,6.0]
p2=[-4.0,4.0]
p3=[-3.0,3.0]
p4=[-3.0,-2.0]
p5=[-2.0,2.0]
p6=[3.08,6.87]
p7=[5.0,5.0]
p8=[3.0,0.0]
p9=[4.0,-4.0]
p10=[0.0,-1.0]
p11=[-1.16,4.37]
p12=[1.94,1.55]
pts=[p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12]
T,adj,adj2v,DG=DT.triangulate_bowyer(pts)
idx = DT.convex_hull(DG, pts)
@test idx == [4, 9, 7, 6, 1, 2]

p11=[0.0,8.0]
pts=[p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12]
T,adj,adj2v,DG=DT.triangulate_bowyer(pts)
idx = DT.convex_hull(DG, pts)
@test idx == [4,9,7,6,11,1,2]
