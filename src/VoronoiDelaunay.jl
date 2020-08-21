# Fast, robust 2D Voronoi/Delaunay tesselation
# Implementation follows algorithms described in http://arxiv.org/abs/0901.4107
# and used (for e.g.) in the Illustris Simulation
# http://www.illustris-project.org/
#
# Author: Ariel Keselman (skariel@gmail.com)
# License: MIT
# Bug reports welcome!

module VoronoiDelaunay
    using GeometricalPredicates
    const GP = GeometricalPredicates
    import GeometricalPredicates: geta, getb, getc

    import Base: push!, iterate, copy, sizehint!
    import Colors: RGB, RGBA
    using Random: shuffle!

    #TODO: What is up with adding ±eps(Float64)?
    const min_coord = GeometricalPredicates.min_coord + eps(Float64)
    const max_coord = GeometricalPredicates.max_coord - eps(Float64)

    include("exports.jl")
    include("scaling.jl")
    include("datastructures.jl")
    include("tesselation.jl")
    include("images.jl")
    include("plotting.jl")
    include("delaunay2.jl")

end 
