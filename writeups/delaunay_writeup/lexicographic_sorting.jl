using Test
x = [1, 3, 4, 1, 2, 5, 10, 3]
y = [5, 1, 2, 4, 10, 2, 3, 3]
p = [[x, y] for (x, y) in zip(x, y)]
p_ls = [[1, 4], [1, 5], [2, 10], [3, 1],
    [3, 3], [4, 2], [5, 2], [10, 3]] # True lexicographic sorting
sort!(p)
@test p == p_ls
