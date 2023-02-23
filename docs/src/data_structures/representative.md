```@meta
CurrentModule = DelaunayTriangulation
```

# Representative Coordinates

When we are building the triangulations, we need to have a point that represents the center of the domain at each stage so that ghost triangles can be worked with for point location. The main type we use for this is `RepresentativeCoordinates`:

```@docs 
RepresentativeCoordinates
```

More accurately, we maintain a list of such coordinates, one for each boundary curve. This list is a `const`, given below.

```@docs 
RepresentativePointList
```

This list is updated at the end of a triangulation to be given by the poles of inaccessibilities (see the sidebar), but when the triangulation is being built it is given by the arithmetic mean of all points in the triangulation. Some other docstrings are below.

```@docs 
get_representative_point_coordinates 
compute_representative_points!(::Any, ::Any)
```