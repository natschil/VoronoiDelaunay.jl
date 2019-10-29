module VoronoiDelaunay

# Fast, robust 2D Voronoi/Delaunay tesselation
# Implementation follows algorithms described in http://arxiv.org/abs/0901.4107
# and used (for e.g.) in the Illustris Simulation
# http://www.illustris-project.org/
#
# Author: Ariel Keselman (skariel@gmail.com)
# License: MIT
# Bug reports welcome!

export
DelaunayTessellation,
delaunayTriangles,
delaunayEdges,
Vertices,
DelaunayTessellation2D, 
sizehint!, 
isexternal,
min_coord, 
max_coord, 
locate, 
movea, 
moveb, 
movec,
delaunayedges, 
voronoiedges, 
voronoiedgeswithoutgenerators,
iterate, 
findindex, 
push!,
Point, 
Point2D, 
AbstractPoint2D, 
getx, 
gety, 
geta, 
getb, 
getc,
getgena, 
getgenb, 
getplotxy,
scaleShiftPoints,
expand,
Triple,
DelaunayTriangle



using GeometricalPredicates
import GeometricalPredicates: geta, getb, getc

import Base: push!, iterate, copy, sizehint!
import Colors: RGB, RGBA
using Random: shuffle!

using OffsetArrays

include("VoronoiDelaunayExtensions.jl")

const min_coord = GeometricalPredicates.min_coord + eps(Float64)
const max_coord = GeometricalPredicates.max_coord - eps(Float64)

const Triple = NTuple{3,Int64}   



mutable struct DelaunayTriangle{T<:AbstractPoint2D} <: AbstractNegativelyOrientedTriangle
    _a::T; _b::T; _c::T
    _bx::Float64; _by::Float64
    _cx::Float64; _cy::Float64
    _px::Float64; _py::Float64
    _pr2::Float64
    n::Triple                   # neighbours
    v::Triple                   # vertices

    function DelaunayTriangle{T}(pa::T, pb::T, pc::T,
                                 n::Triple, v::Triple) where T
        t = new(pa, pb, pc, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, n, v)
        clean!(t)
        t
    end

    # this constructor is good for making copies
    function DelaunayTriangle{T}(pa::T, pb::T, pc::T,
                                 bx::Float64, by::Float64,
                                 cx::Float64, cy::Float64,
                                 px::Float64, py::Float64,
                                 pr2::Float64,
                                 n::Triple, v::Triple) where T
        new(pa, pb, pc, bx, by, cx, cy, px, py, pr2, n, v)
    end
end

function DelaunayTriangle(pa::T, pb::T, pc::T,
                          n::Triple, v::Triple) where T<:AbstractPoint2D
    DelaunayTriangle{T}(pa, pb, pc, n, v)
end

function DelaunayTriangle(pa::T, pb::T, pc::T,
                          bx::Float64, by::Float64,
                          cx::Float64, cy::Float64,
                          px::Float64, py::Float64,
                          pr2::Float64,
                          n::Triple, v::Triple) where T<:AbstractPoint2D
    DelaunayTriangle{T}(pa, pb, pc, bx, by, cx, cy, px, py, pr2, n, v)
end

function copy(t::DelaunayTriangle{T}) where T<:AbstractPoint2D
    DelaunayTriangle( t._a, t._b, t._c,
                      t._bx, t._by,
                      t._cx, t._cy,
                      t._px, t._py,
                      t._pr2,
                      t.n, t.v )
end



mutable struct DelaunayTessellation2D{T<:AbstractPoint2D}
    _trigs::Vector{DelaunayTriangle{T}}
    _workpoints::OffsetVector{T}
    _last_trig_index::Int64
    _edges_to_check::Vector{Tuple{Int64,Int64}}
    _total_points_added::Int64

    function DelaunayTessellation2D{T}(n::Int64 = 100) where T
        a = T(GeometricalPredicates.min_coord, GeometricalPredicates.min_coord)
        b = T(GeometricalPredicates.min_coord, GeometricalPredicates.max_coord)
        c = T(GeometricalPredicates.max_coord, GeometricalPredicates.min_coord)
        d = T(GeometricalPredicates.max_coord, GeometricalPredicates.max_coord)
        t1 = DelaunayTriangle{T}( d,c,b, (2,1,1), (0,-1,-2) )
        t2 = DelaunayTriangle{T}( a,b,c, (3,1,1), (-3,-2,-1) )
        t3 = DelaunayTriangle{T}( d,c,b, (2,1,1), (0,-1,-2) )
        _trigs = DelaunayTriangle{T}[t1, t2, t3]
        _workpoints = OffsetVector([a,b,c,d],-3:0)
        t = new(_trigs, _workpoints, 3, Int64[], 0)
        sizehint!(t._edges_to_check, 1000)
        sizehint!(t, n)
    end
end

DelaunayTessellation2D(n::Int64) = DelaunayTessellation2D{Point2D}(n)

DelaunayTessellation2D(n::Int64, ::T) where {T<:AbstractPoint2D} = DelaunayTessellation2D{T}(n)

DelaunayTessellation(n::Int64=100) = DelaunayTessellation2D(n)

# tessellation function that can deal with points outside [1,2]x[1,2]
# and that always returns a convex tessellation
function DelaunayTessellation( points::Array{Point2D,1} )
    scaledPoints, ranges = scaleShiftPoints( points )
    scaledTess = DelaunayTessellation( length( points ) )
    push!( scaledTess, scaledPoints )  
    tess = expand( scaledTess, ranges )
    return tess
end



# filter out triangles whose vertices are all from the given point set
delaunayTriangles( trigs::Vector{DelaunayTriangle{T}} ) where T<:AbstractPoint2D = filter( x -> x.v[1]>0 && x.v[2]>0 && x.v[3]>0, trigs ) 

# edges given as a vector of sets of two point indices
function delaunayEdges( trigs::Vector{DelaunayTriangle{T}} ) where T<:AbstractPoint2D
    edges = Vector{Set{Int64}}()
    for t in trigs
        e1 = Set(t.v[1:2]); e2 = Set(t.v[2:3]); e3 = Set(( t.v[3], t.v[1] ));
        union!( edges, [ e1, e2, e3 ] )
    end
    return edges
end

Vertices( tess::DelaunayTessellation2D{T} ) where T<:AbstractPoint2D = tess._workpoints[1:end]

function getplotxy( Edges::Vector{Set{Int64}}, Points::Vector{T} ) where T<:AbstractPoint2D
    x = Float64[]
    y = Float64[]
    for e in Edges
        ab = Points[collect(e)]
        push!(x, ab[1]._x, ab[2]._x, NaN)
        push!(y, ab[1]._y, ab[2]._y, NaN)
    end
    (x, y)
end



function sizehint!(t::DelaunayTessellation2D{T}, n::Int64) where T<:AbstractPoint2D
    required_total_size = 2n + 10
    required_total_size <= length(t._trigs) && return
    sizehint!(t._trigs, required_total_size)
    while length(t._trigs) < required_total_size
        push!(t._trigs, copy(t._trigs[end]))
    end
    t
end

# convert the tessellation back to original scale after the tessellation
function expand( tess::DelaunayTessellation2D{T}, ranges::NTuple{4,Float64} ) where T<:AbstractPoint2D
    xmin = ranges[1]
    ymin = ranges[3]
    scale = max( ranges[4] - ranges[3], ranges[2] - ranges[1] ) / 0.98
    offset = 1.01  
    tess._workpoints = [ Point2D( ( p._x - offset ) * scale + xmin, ( p._y - offset ) * scale + ymin ) for p in tess._workpoints ]
    for t in tess._trigs
        t._a = Point2D( ( t._a._x - offset ) * scale + xmin, ( t._a._y - offset ) * scale + ymin )
        t._b = Point2D( ( t._b._x - offset ) * scale + xmin, ( t._b._y - offset ) * scale + ymin )
        t._c = Point2D( ( t._c._x - offset ) * scale + xmin, ( t._c._y - offset ) * scale + ymin )
        t._bx *= scale
        t._by *= scale
        t._cx *= scale
        t._cy *= scale
        t._px *= scale^3
        t._py *= scale^3
        t._pr2 *= scale^2
    end
    return tess
end



# push an array but sort it first for better performance
function push!(tess::DelaunayTessellation2D{T}, a::Array{T, 1}) where T<:AbstractPoint2D
    shuffle!(a)
    mssort!(a)
    _pushunsorted!(tess, a)
end

# push an array in given order
function _pushunsorted!(tess::DelaunayTessellation2D{T}, a::Array{T, 1}) where T<:AbstractPoint2D
    sizehint!(tess, length(a))
    for p in a
        push!( tess._workpoints, p )
        push!( tess, p, last(axes(tess._workpoints,1)) )
    end
end

# push a single point. Grows tessellation as needed
function push!( tess::DelaunayTessellation2D{T}, p::T, v::Int64 ) where T<:AbstractPoint2D
    tess._total_points_added += 1
    sizefit_at_least(tess, tess._total_points_added)
    i = _pushunfixed!( tess, p, v )
    _restoredelaunayhood!(tess, i)
end

function _pushunfixed!( tess::DelaunayTessellation2D{T}, p::T, v::Int64 ) where T<:AbstractPoint2D

    # indices of new triangles
    i = findindex(tess, p)
    ltrigs1::Int64 = tess._last_trig_index + 1
    ltrigs2::Int64 = tess._last_trig_index + 2

    # new triangles
    @inbounds t1 = tess._trigs[i]
    @inbounds t2 = tess._trigs[ltrigs1]
    @inbounds t3 = tess._trigs[ltrigs2]

    # update point indices
    t3.v = ( t1.v[1], t1.v[2], v )
    t2.v = ( t1.v[1], v, t1.v[3] )
    t1.v = ( v, t1.v[2], t1.v[3] )

    # update neighbours
    t3.n = ( i, ltrigs1, t1.n[3] )
    t2.n = ( i, t1.n[2], ltrigs2 )
    t1.n = ( t1.n[1], ltrigs1, ltrigs2 )

    # update point values
    t1._a, t1._b, t1._c = tess._workpoints[t1.v[1]], tess._workpoints[t1.v[2]], tess._workpoints[t1.v[3]]
    t2._a, t2._b, t2._c = tess._workpoints[t2.v[1]], tess._workpoints[t2.v[2]], tess._workpoints[t2.v[3]]
    t3._a, t3._b, t3._c = tess._workpoints[t3.v[1]], tess._workpoints[t3.v[2]], tess._workpoints[t3.v[3]]
    
    # update precomputed determinants etc.
    clean!(t1)
    clean!(t2)
    clean!(t3)

    # update neighbours of neighbours
    @inbounds nt2 = tess._trigs[t2.n[2]]
    index = findfirst( x -> x == i , nt2.n )
    if index != nothing
        nt2.n = Base.setindex( nt2.n, ltrigs1, index )
    end
    @inbounds nt3 = tess._trigs[t3.n[3]]
    index = findfirst( x -> x == i , nt3.n )
    if index != nothing
        nt3.n = Base.setindex( nt3.n, ltrigs2, index )
    end

    tess._last_trig_index += 2

    i
end

function _restoredelaunayhood!(tess::DelaunayTessellation2D{T},
                               ix_trig::Int64) where T<:AbstractPoint2D

    ltrig1 = tess._last_trig_index - 1
    ltrig2 = tess._last_trig_index
    # A edge
    push!( tess._edges_to_check, ( ix_trig, tess._trigs[ix_trig].n[1] ) )
    # B edge
    push!( tess._edges_to_check, ( ltrig1, tess._trigs[ltrig1].n[2] ) )
    # C edge
    push!( tess._edges_to_check, ( ltrig2, tess._trigs[ltrig2].n[3] ) )

    while !isempty(tess._edges_to_check)
        t_pair = pop!( tess._edges_to_check )
        if t_pair[2] > 1
            @inbounds t1 = tess._trigs[t_pair[1]]
            @inbounds t2 = tess._trigs[t_pair[2]]
            p = tess._workpoints[ setdiff( t1.v, t2.v )[1] ]
            if incircle( t2, p ) > 0
                i1, i2 = flip!( tess, t_pair[1], t_pair[2] )
                push!( tess._edges_to_check, ( t_pair[1], i1 ) )
                push!( tess._edges_to_check, ( t_pair[2], i2 ) )
            end
        end
    end

end

# generalised flip function
function flip!( tess::DelaunayTessellation2D, i1::Int64, i2::Int64 )

    # triangles with a common edge to be flipped
    t1 = tess._trigs[i1]
    t2 = tess._trigs[i2]

    # to generalise the flipping functions, circshift the neighbours and vertices lists
    # such that the triangles lists begin with the triangle opposite the edge to be flipped
    cs1 = 1 - findfirst( x -> x == i2, t1.n )
    cs2 = 1 - findfirst( x -> x == i1, t2.n )
    t1.n = circshiftTuple( t1.n, cs1 )
    t2.n = circshiftTuple( t2.n, cs2 )
    t1.v = circshiftTuple( t1.v, cs1 )
    t2.v = circshiftTuple( t2.v, cs2 )

    # new neighbours
    t_union = ( t1.n[2:3]..., t2.n[2:3]... )
    t_union = circshiftTuple( t_union, 1 )
    t1.n = ( t_union[1:2]..., i2 )
    t2.n = ( t_union[3:4]..., i1 )

    # new vertices
    v_union = ( t1.v[1:2]..., t2.v[1:2]... ) 
    t1.v = ( t1.v[1], v_union[3:4]... )
    t2.v = ( t2.v[1], v_union[1:2]... )

    # new neighbours of neighbours
    n1 = tess._trigs[t1.n[1]]
    n2 = tess._trigs[t2.n[1]]
    index = findfirst( x -> x == i2, n1.n )
    if index != nothing
        n1.n = Base.setindex( n1.n, i1, index )
    end
    index = findfirst( x -> x == i1, n2.n )
    if index != nothing
        n2.n = Base.setindex( n2.n, i2, index )
    end

    # new point values
    t1._a, t1._b, t1._c = tess._workpoints[t1.v[1]], tess._workpoints[t1.v[2]], tess._workpoints[t1.v[3]]
    t2._a, t2._b, t2._c = tess._workpoints[t2.v[1]], tess._workpoints[t2.v[2]], tess._workpoints[t2.v[3]]

    # update precomputed determinants etc.
    clean!(t1)
    clean!(t2)

    # new edges to check: ( i1, t1.n[1] ) and ( i2, t2.n[2] )
    return t1.n[1], t2.n[2]

end

function circshiftTuple( src::Tuple{Vararg{Int64}}, shift::Int64 )
    len = length(src)
    split = len - mod( shift, len )
    r1, r2 = 1:split, split+1:len
    dest = ( src[r2]..., src[r1]... )
end

function findindex(tess::DelaunayTessellation2D{T}, p::T) where T<:AbstractPoint2D
    i::Int64 = tess._last_trig_index
    while true
        @inbounds w = intriangle(tess._trigs[i], p)
        w > 0 && return i
        @inbounds tr = tess._trigs[i]
        if w == -1
            i = tr.n[1]
        elseif w == -2
            i = tr.n[2]
        else
            i = tr.n[3]
        end
    end
end

locate(t::DelaunayTessellation2D{T}, p::T) where {T<:AbstractPoint2D} = t._trigs[findindex(t, p)]

movea(tess::DelaunayTessellation2D{T},
      trig::DelaunayTriangle{T}) where {T<:AbstractPoint2D} = tess._trigs[trig.n[1]]
moveb(tess::DelaunayTessellation2D{T},
      trig::DelaunayTriangle{T}) where {T<:AbstractPoint2D} = tess._trigs[trig.n[2]]
movec(tess::DelaunayTessellation2D{T},
      trig::DelaunayTriangle{T}) where {T<:AbstractPoint2D} = tess._trigs[trig.n[3]]

# growing strategy
function sizefit_at_least(t::DelaunayTessellation2D{T}, n::Int64) where T<:AbstractPoint2D
    minimal_acceptable_actual_size = 2*n+10
    minimal_acceptable_actual_size <= length(t._trigs) && return
    required_total_size = length(t._trigs)
    while required_total_size < minimal_acceptable_actual_size
        required_total_size += required_total_size >>> 1
    end
    sizehint!(t._trigs, required_total_size)
    while length(t._trigs) < required_total_size
        push!(t._trigs, copy(t._trigs[end]))
    end
    t
end



mutable struct TrigIter
    ix::Int64
end

function iterate( t::DelaunayTessellation2D, it::TrigIter=TrigIter(2), ranges::NTuple{4,Float64}=(min_coord,max_coord,min_coord,max_coord) ) # default it resembles old start
    # resembles old done
    while isexternal(t._trigs[it.ix], ranges) && it.ix <= t._last_trig_index
        it.ix += 1
    end
    if it.ix > t._last_trig_index
        return nothing
    end
    # resembles old next
    @inbounds trig = t._trigs[it.ix]
    it.ix += 1
    return (trig, it)
end

struct DelaunayEdge{T<:AbstractPoint2D}
    _a::T
    _b::T
end
geta(e::DelaunayEdge{T}) where {T<:AbstractPoint2D} = e._a
getb(e::DelaunayEdge{T}) where {T<:AbstractPoint2D} = e._b

struct VoronoiEdge{T<:AbstractPoint2D}
    _a::Point2D
    _b::Point2D
    _generator_a::T
    _generator_b::T
end
geta(e::VoronoiEdge{T}) where {T<:AbstractPoint2D} = e._a
getb(e::VoronoiEdge{T}) where {T<:AbstractPoint2D} = e._b
getgena(e::VoronoiEdge{T}) where {T<:AbstractPoint2D} = e._generator_a
getgenb(e::VoronoiEdge{T}) where {T<:AbstractPoint2D} = e._generator_b

struct VoronoiEdgeWithoutGenerators
    _a::Point2D
    _b::Point2D
end
geta(e::VoronoiEdgeWithoutGenerators) = e._a
getb(e::VoronoiEdgeWithoutGenerators) = e._b

# TODO: is an iterator faster?
function delaunayedges( t::DelaunayTessellation2D, ranges::NTuple{4,Float64}=(min_coord,max_coord,min_coord,max_coord) )
    visited = zeros(Bool, t._last_trig_index)
    function delaunayiterator(c::Channel)
        @inbounds for ix in 2:t._last_trig_index
            tr = t._trigs[ix]
            isexternal( tr, ranges ) && continue
            visited[ix] && continue
            visited[ix] = true
            if !visited[tr.n[1]]
                put!(c, DelaunayEdge(getb(tr), getc(tr)))
            end
            if !visited[tr.n[2]]
                put!(c, DelaunayEdge(geta(tr), getc(tr)))
            end
            if !visited[tr.n[3]]
                put!(c, DelaunayEdge(geta(tr), getb(tr)))
            end
        end
    end
    Channel(delaunayiterator)
end

# TODO: is an iterator faster?
function voronoiedges(t::DelaunayTessellation2D)
    visited = zeros(Bool, t._last_trig_index)
    visited[1] = true
    function voronoiiterator(c::Channel)
        for ix in 2:t._last_trig_index
            visited[ix] && continue
            tr = t._trigs[ix]
            visited[ix] = true
            #isexternal(tr) && continue
            cc = circumcenter(tr)

            if !visited[tr.n[1]] #&& !isexternal(t._trigs[ix_na])
                nb = t._trigs[tr.n[1]]
                put!(c, VoronoiEdge(cc, circumcenter(nb), getb(tr), getc(tr)))
            end
            if !visited[tr.n[2]] #&& !isexternal(t._trigs[ix_nb])
                nb = t._trigs[tr.n[2]]
                put!(c, VoronoiEdge(cc, circumcenter(nb), geta(tr), getc(tr)))
            end
            if !visited[tr.n[3]] #&& !isexternal(t._trigs[ix_nb])
                nb = t._trigs[tr.n[3]]
                put!(c, VoronoiEdge(cc, circumcenter(nb), geta(tr), getb(tr)))
            end
        end
    end
    Channel(voronoiiterator)
end

function voronoiedgeswithoutgenerators(t::DelaunayTessellation2D)
    visited = zeros(Bool, t._last_trig_index)
    visited[1] = true
    function voronoiiterator(c::Channel)
        for ix in 2:t._last_trig_index
            visited[ix] && continue
            tr = t._trigs[ix]
            visited[ix] = true
            #isexternal(tr) && continue
            cc = circumcenter(tr)

            if !visited[tr.n[1]] #&& !isexternal(t._trigs[ix_na])
                nb = t._trigs[tr.n[1]]
                put!(c, VoronoiEdgeWithoutGenerators(cc, circumcenter(nb)))
            end
            if !visited[tr.n[2]] #&& !isexternal(t._trigs[ix_nb])
                nb = t._trigs[tr.n[2]]
                put!(c, VoronoiEdgeWithoutGenerators(cc, circumcenter(nb)))
            end
            if !visited[tr.n[3]] #&& !isexternal(t._trigs[ix_nc])
                nb = t._trigs[tr.n[3]]
                put!(c, VoronoiEdgeWithoutGenerators(cc, circumcenter(nb)))
            end
        end
    end
    Channel(voronoiiterator)
end

function isexternal( t::DelaunayTriangle{T}, ranges::NTuple{4,Float64} ) where T<:AbstractPoint2D
    getx(geta(t)) < ranges[1] || getx(geta(t)) > ranges[2] ||
    getx(getb(t)) < ranges[1] || getx(getb(t)) > ranges[2] ||
    getx(getc(t)) < ranges[1] || getx(getc(t)) > ranges[2]
end



intensity(c::RGB)  = c.b
intensity(c::RGBA) = c.b
intensity(c) 	     = getfield(c, 1) # Workaround. Gray needs to be imported from images, which would take to long.

# Create DelaunayTessellation with npts points from an image
function from_image(img, npts)

    # placing points in places that represent the image
    pts = Point2D[]
    for i in 1:npts
        x = rand()
        y = rand()
        if intensity(img[Int64(floor(x * size(img)[1])) + 1, Int64(floor(y * size(img)[2])) + 1]) > 0.5
            if rand() < 0.100
                push!(pts, Point2D(1.0 + rand(), 1.0 + rand()))
            end
            continue
        end
        x /= 2
        y /= 2
        push!(pts, Point2D(1.0 + x + 0.3, 2.0 - y*2/3 - 0.3))
    end

    tess = DelaunayTessellation(npts)
    push!(tess, pts)

    tess
end

function getplotxy(edges)
    edges = collect(edges)
    x = Float64[]
    y = Float64[]
    for e in edges
        push!(x, getx(geta(e)))
        push!(x, getx(getb(e)))
        push!(x, NaN)
        push!(y, gety(geta(e)))
        push!(y, gety(getb(e)))
        push!(y, NaN)
    end
    (x, y)
end

end # module VoronoiDelaunay