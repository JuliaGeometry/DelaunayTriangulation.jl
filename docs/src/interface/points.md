```@meta
CurrentModule = DelaunayTriangulation
```

## Individual Points 

Points are assumed to take the form `(x, y)`, but we allow for customisation in how we represent these points. The following methods are used, where we first list the methods that must be defined and then methods that extend these former methods. 

### Necessary Methods 

```@docs 
getx 
gety 
```

### Generic Methods 

```@docs 
getxy 
```

## Collection of Points 

A collection of points simply store many points. It does not need to be mutable (unless you want to add points into the triangulation not already in `tri.points`, or if you want mesh refinement). The following methods are used, where we first list the methods that must be defined and then methods that extend these former methods. 

### Necessary Methods 

```@docs 
getpoint 
each_point_index 
each_point
num_points 
number_type 
push_point!
pop_point!
```

### Generic Methods 

```@docs 
get_point 
points_are_unique 
lexicographic_order 
```