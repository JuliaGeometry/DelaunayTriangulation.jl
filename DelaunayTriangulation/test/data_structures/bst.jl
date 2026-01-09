using ..DelaunayTriangulation
using Test
using DataStructures
const DT = DelaunayTriangulation

tree = DT.BalancedBST{Int}()
tree_avl = AVLTree{Int}()
for _ in 1:1000000
    i = rand(1:10000)
    push!(tree, i)
    push!(tree_avl, i)
end
@test tree.count == tree_avl.count
@test DT.get_height(tree.root) == DataStructures.get_height(tree_avl.root)
@test DT.inorder(tree) == __inorder(tree_avl)
@test length(DT.inorder(tree)) == tree.count
inord = DT.inorder(tree)
for i in [inord; 0; -1;-2;-3]
    delete!(tree, i)
    delete!(tree_avl, i)
end
@test tree.count == tree_avl.count
@test DT.get_height(tree.root) == DataStructures.get_height(tree_avl.root)
@test DT.inorder(tree) == __inorder(tree_avl)
@test length(DT.inorder(tree)) == tree.count

tree = DT.BalancedBST{NTuple{2, Int}}()
tree_avl = AVLTree{NTuple{2, Int}}()
for _ in 1:100000
    i = (rand(1:10000), rand(1:10000))
    push!(tree, i)
    push!(tree_avl, i)
end
@test tree.count == tree_avl.count
@test DT.get_height(tree.root) == DataStructures.get_height(tree_avl.root)
@test DT.inorder(tree) == __inorder(tree_avl)
@test length(DT.inorder(tree)) == tree.count
inord = DT.inorder(tree)
for i in [inord; (0, 0); (-1, -1);(-2, -2);(-3, -3)]
    delete!(tree, i)
    delete!(tree_avl, i)
end
@test tree.count == tree_avl.count
@test DT.get_height(tree.root) == DataStructures.get_height(tree_avl.root)
@test DT.inorder(tree) == __inorder(tree_avl)
@test length(DT.inorder(tree)) == tree.count
