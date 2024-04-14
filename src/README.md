This folder contains the source code for DelaunayTriangulation.jl. Since there are quite a large number of files, it may be a bit daunting to navigate. Here is a brief description of the organisation of the source code.

# `src/interfaces`

This folder defines the basic interfaces for representing points, triangles, edges, and boundaries. These interfaces are what `Triangulation`s rely on to do their work.

# `src/structs`

This folder contains all the definitions for the structs used in this package. The name of each file is the name of the struct defined in it. 