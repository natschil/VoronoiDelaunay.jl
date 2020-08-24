#Struct for representing a single triangle

mutable struct DelaunayTriangle{T<:AbstractPoint2D} <: AbstractNegativelyOrientedTriangle
    #The three points of the triangle
    #points::Array{T}
    point_a::T
    point_b::T
    point_c::T

    #Other constants required by GeometricalPredicates.jl
    _bx::Float64; _by::Float64
    _cx::Float64; _cy::Float64

    _px::Float64; _py::Float64
    _pr2::Float64

    #Indexes of neighboring points
    #neighbors::Array{Int64}
    neighbor_a::Int64
    neighbor_b::Int64
    neighbor_c::Int64

    function DelaunayTriangle{T}(pa::T, pb::T, pc::T,
                                 na::Int64, nb::Int64, nc::Int64) where {T<:AbstractPoint2D}
        t = new(pa, pb, pc, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, na, nb, nc)
        clean!(t)
        return t
    end

    # this constructor is good for making copies
    function DelaunayTriangle{T}(pa::T, pb::T, pc::T,
                                 bx::Float64, by::Float64,
                                 cx::Float64, cy::Float64,
                                 px::Float64, py::Float64,
                                 pr2::Float64,
                                 na::Int64, nb::Int64, nc::Int64) where {T<:AbstractPoint2D}
        new(pa, pb, pc, bx, by, cx, cy, px, py, pr2, na, nb, nc)
    end
end

#Get/set methods
GP.geta(tri::DelaunayTriangle{T}) where T<:AbstractPoint2D = tri.point_a
GP.getb(tri::DelaunayTriangle{T}) where T<:AbstractPoint2D = tri.point_b
GP.getc(tri::DelaunayTriangle{T}) where T<:AbstractPoint2D = tri.point_c
GP.seta(tri::DelaunayTriangle{T},val::T) where T<:AbstractPoint2D = (tri.point_a = val; clean!(tri))
GP.setb(tri::DelaunayTriangle{T},val::T) where T<:AbstractPoint2D = (tri.point_b = val;clean!(tri))
GP.setc(tri::DelaunayTriangle{T},val::T) where T<:AbstractPoint2D = (tri.point_c = val;clean!(tri))
function GP.setabc(tri::DelaunayTriangle{T},vala::T,valb::T,valc::T) where T<:AbstractPoint2D
    tri.point_a = vala
    tri.point_b = valb
    tri.point_c = valc
    clean!(tri)
    return tri
end
function getpoint(tri::DelaunayTriangle{T}, index::Int64) where T<:AbstractPoint2D
    if index == 1
        return tri.point_a
    elseif index == 2
        return tri.point_b
    else
        return tri.point_c
    end
end
function setpoint!(tri::DelaunayTriangle{T}, index::Int64,val::T) where T<:AbstractPoint2D
    if index == 1
        tri.point_a = val
    elseif index == 2
        tri.point_b = val
    else
        tri.point_c = val
    end
    clean!(tri)
end

function getneighbor(tri::DelaunayTriangle{T}, index::Int64) where T<:AbstractPoint2D 
    if index == 1
        return tri.neighbor_a
    elseif index == 2
        return tri.neighbor_b
    else
        return tri.neighbor_c
    end
end
function setneighbor!(tri::DelaunayTriangle{T}, index::Int64,val::Int64) where T<:AbstractPoint2D
    if index == 1
        tri.neighbor_a = val
    elseif index == 2
        tri.neighbor_b = val
    else
        tri.neighbor_c = val
    end
end

#=TODO: Figure out if the two functions below are needed or not...
function DelaunayTriangle(pa::T, pb::T, pc::T,
                          na::Int64, nb::Int64, nc::Int64) where T<:AbstractPoint2D
DelaunayTriangle{T}(pa, pb, pc, na, nb, nc)
end

function DelaunayTriangle(pa::T, pb::T, pc::T,
                          bx::Float64, by::Float64,
                          cx::Float64, cy::Float64,
                          px::Float64, py::Float64,
                          pr2::Float64,
                          na::Int64, nb::Int64, nc::Int64) where T<:AbstractPoint2D
    DelaunayTriangle{T}(pa, pb, pc, bx, by, cx, cy, px, py, pr2, na, nb, nc)
end
=#

function copy(t::DelaunayTriangle{T}) where T<:AbstractPoint2D
    DelaunayTriangle{T}(
                     getpoint(t,1),getpoint(t,2),getpoint(t,3),
                     t._bx, t._by,
                     t._cx, t._cy,
                     t._px, t._py,
                     t._pr2,
                     getneighbor(t,1),getneighbor(t,2),getneighbor(t,3)
                     )
end

function isexternal(t::T) where T<:AbstractPoint2D
    return gety(t) < min_coord || gety(t) > max_coord
end

function isexternal( t::DelaunayTriangle{T}) where T<:AbstractPoint2D
    return isexternal(geta(t)) || isexternal(getb(t)) || isexternal(getc(t))
end


mutable struct DelaunayTessellation2D{T<:AbstractPoint2D}
    _trigs::Vector{DelaunayTriangle{T}}
    _nonexternal_indexes::Vector{Int64} #_nonexternal_indexes[i] gives index of ith non external triangle
    _nonexternal_indexes_inv::Vector{Int64} # inverse of the above
    _last_trig_index::Int64 #How many of the elements in _trigs are actually being used
    _prev_trig_index::Int64 #The triangle that the previous search ended in
    _edges_to_check::Vector{Int64} #A buffer used by _restoredelaunayhood
    _total_points_added::Int64 # Number of points added to the tesselation
    _scaling::NTuple{4,Float64} # xmin,xwidth,ymin,yw of the domain used (relative to (min_coord, max_coord-min_coord) in each direction)
    _invscaling::NTuple{4,Float64}
    _convex::Bool #Triangulation has 3 points "at infinity" to enforce convexity

    function DelaunayTessellation2D{T}(n::Int64 = 100;convex=false) where T
        a = T(GeometricalPredicates.min_coord, GeometricalPredicates.min_coord)
        b = T(GeometricalPredicates.min_coord, GeometricalPredicates.max_coord)
        c = T(GeometricalPredicates.max_coord, GeometricalPredicates.min_coord)
        d = T(GeometricalPredicates.max_coord, GeometricalPredicates.max_coord)
        e = T(0.5*(GeometricalPredicates.min_coord + GeometricalPredicates.max_coord), GeometricalPredicates.max_coord)
        _scaling=(0.0,1.0,0.0,1.0)
        if !convex
            t1 = DelaunayTriangle{T}(d,c,b, 2,1,1) #Here t1 is a "dummy" triangle that is not part of the domain
            t2 = DelaunayTriangle{T}(a,b,c, 3,1,1)
            t3 = DelaunayTriangle{T}(d,c,b, 2,1,1)
            _trigs = DelaunayTriangle{T}[t1, t2, t3]
            _nonexternal_indexes = Int64[0,0,0]
            t = new(_trigs,_nonexternal_indexes,copy(_nonexternal_indexes), 3,2, Int64[], 0,_scaling,_scaling,false)
        else
            t1 = DelaunayTriangle{T}(a,e,c,1,1,1) #Dummy triangle
            t2 = DelaunayTriangle{T}(a,e,c,1,1,1)
            _trigs = DelaunayTriangle{T}[t1, t2]
            _nonexternal_indexes = Int64[0,0,0]
            t = new(_trigs,_nonexternal_indexes, copy(_nonexternal_indexes), 2,2, Int64[], 0,_scaling,_scaling,true)
        end
        sizehint!(t._edges_to_check, 1000)
        sizehint!(t, n)
        return t
    end
end
DelaunayTessellation2D(n::Int64;convex=false) = DelaunayTessellation2D{Point2D}(n;convex=convex)
DelaunayTessellation2D(n::Int64, ::T;convex=false) where {T<:AbstractPoint2D} = DelaunayTessellation2D{T}(n; convex=convex)
DelaunayTessellation(n::Int64=100; convex=false) = DelaunayTessellation2D(n;convex=convex)

function _clean_nonexternal_indexes!(tess)
    j = 1
    for (i,t) in enumerate(tess._trigs)
        if i > tess._last_trig_index
            tess._nonexternal_indexes_inv[i] = 0
            continue
        end
        if isexternal(t)
            tess._nonexternal_indexes_inv[i] = 0
            continue
        else
            tess._nonexternal_indexes[j] = i
            tess._nonexternal_indexes_inv[i] = j
            j += 1
        end
    end
    for i in j:length(tess._nonexternal_indexes)
        tess._nonexternal_indexes[i] = 0
    end
    return 0
end

function sizehint!(t::DelaunayTessellation2D{T}, n::Int64) where T<:AbstractPoint2D
    required_total_size = 2n + 10
    required_total_size <= length(t._trigs) && return
    sizehint!(t._trigs, required_total_size)
    sizehint!(t._nonexternal_indexes, required_total_size)
    while length(t._trigs) < required_total_size
        push!(t._trigs, copy(t._trigs[end]))
        push!(t._nonexternal_indexes, 0)
        push!(t._nonexternal_indexes_inv, 0)
    end
    return t
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
function delaunayedges( t::DelaunayTessellation2D,include_external=false)
    sc = t._invscaling
    visited = zeros(Bool, t._last_trig_index)
    function delaunayiterator(c::Channel)
        @inbounds for ix in 2:t._last_trig_index
            tr = t._trigs[ix]
            !include_external && isexternal( tr ) && continue
            visited[ix] && continue
            visited[ix] = true
            ix_na = getneighbor(tr,1)
            if !visited[ix_na]
                put!(c, DelaunayEdge(rescale(getb(tr),sc), rescale(getc(tr),sc)))
            end
            ix_nb = getneighbor(tr,2)
            if !visited[ix_nb]
                put!(c, DelaunayEdge(rescale(geta(tr),sc), rescale(getc(tr),sc)))
            end
            ix_nc = getneighbor(tr,3)
            if !visited[ix_nc]
                put!(c, DelaunayEdge(rescale(geta(tr),sc), rescale(getb(tr),sc)))
            end
        end
    end
    Channel(delaunayiterator)
end

# TODO: is an iterator faster?
function voronoiedges(t::DelaunayTessellation2D)
    if tess._convex
        error("voronoiedges has not been implemented for tess._convex == true")
    end
    visited = zeros(Bool, t._last_trig_index)
    visited[1] = true
    function voronoiiterator(c::Channel)
        for ix in 2:t._last_trig_index
            visited[ix] && continue
            tr = t._trigs[ix]
            visited[ix] = true
            #isexternal(tr) && continue
            cc = circumcenter(tr)

            ix_na = getneighbor(tr,1)
            if !visited[ix_na] #&& !isexternal(t._trigs[ix_na])
                nb = t._trigs[ix_na]
                put!(c, VoronoiEdge(cc, circumcenter(nb), getb(tr), getc(tr)))
            end
            ix_nb = getneighbor(tr,2)
            if !visited[ix_nb] #&& !isexternal(t._trigs[ix_nb])
                nb = t._trigs[ix_nb]
                put!(c, VoronoiEdge(cc, circumcenter(nb), geta(tr), getc(tr)))
            end
            ix_nc = getneighbor(tr,3)
            if !visited[ix_nc] #&& !isexternal(t._trigs[ix_nb])
                nb = t._trigs[ix_nc]
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

            ix_na = getneighbor(tr,1)
            if !visited[ix_na] #&& !isexternal(t._trigs[ix_na])
                nb = t._trigs[ix_na]
                put!(c, VoronoiEdgeWithoutGenerators(cc, circumcenter(nb)))
            end
            ix_nb = getneighbor(tr,2)
            if !visited[ix_nb] #&& !isexternal(t._trigs[ix_nb])
                nb = t._trigs[ix_nb]
                put!(c, VoronoiEdgeWithoutGenerators(cc, circumcenter(nb)))
            end
            ix_nc = getneighbor(tr,3)
            if !visited[ix_nc] #&& !isexternal(t._trigs[ix_nc])
                nb = t._trigs[ix_nc]
                put!(c, VoronoiEdgeWithoutGenerators(cc, circumcenter(nb)))
            end
        end
    end
    Channel(voronoiiterator)
end


mutable struct TrigIter
    ix::Int64
end


function iterate( t::DelaunayTessellation2D, it::TrigIter=TrigIter(2) ) # default it resembles old start
    # resembles old done
    while isexternal(t._trigs[it.ix]) && it.ix <= t._last_trig_index
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

function _findindex(tess::DelaunayTessellation2D{T}, p::S) where {T<:AbstractPoint2D, S<:AbstractPoint2D}
    i::Int64 = tess._prev_trig_index
    #i::Int64 = tess._last_trig_index
    counter::Int64 = 0
    ntriangles = tess._last_trig_index
    while true
        counter += 1
        @inbounds w = intriangle(tess._trigs[i], p)
        if w > 0
            tess._prev_trig_index = i
            return i
        end
        @inbounds tr = tess._trigs[i]
        if w == -1
            i = getneighbor(tr,1)
        elseif w == -2
            i = getneighbor(tr,2)
        else
            i = getneighbor(tr,3)
        end
        if counter > ntriangles
            throw(DomainError("Point was not found in any triangle"))
        end
    end
end

function triangles(tess::DelaunayTessellation2D{T}) where T<:AbstractPoint2D
    res = NTuple{3,T}[]
    sc = tess._invscaling
    for (i,t) in enumerate(tess._trigs)
        if i > tess._last_trig_index
            break
        end
        isexternal(t) && continue
        push!(res,(rescale(geta(t),sc), rescale(getb(t),sc), rescale(getc(t),sc)))
    end
    return res
end

function Base.getindex(tess,i)
        if i > tess._last_trig_index
            error("Out of bounds")
        end
        newindex = tess._nonexternal_indexes[i]
        if newindex == 0
            error("Out of bounds")
        end
        t = tess._trigs[newindex]
        sc = tess._invscaling
        return (rescale(geta(t),sc), rescale(getb(t),sc), rescale(getc(t),sc))
end

function Base.show(io::IO,tess::DelaunayTessellation2D{T}) where {T<:AbstractPoint2D}
    println(io,"DelaunayTesselation2D object with triangles:")
    Base.show(triangles(tess))
end

_locate(t::DelaunayTessellation2D{T}, p::S) where {T<:AbstractPoint2D, S<:AbstractPoint2D} = t._trigs[_findindex(t, p)]

function locate(t::DelaunayTessellation2D{T},p::S) where {T<:AbstractPoint2D, S<:AbstractPoint2D}
    tri = _locate(t,rescale(p,t._scaling))
    if isexternal(tri)
        throw(DomainError("Point is outside of triangulation"))
    end
    return (rescale(geta(tri),t._invscaling),rescale(getb(tri),t._invscaling),rescale(getc(tri),t._invscaling))
end


movea(tess::DelaunayTessellation2D{T},
      trig::DelaunayTriangle{T}) where {T<:AbstractPoint2D} = tess._trigs[getneighbor(trig,1)]
moveb(tess::DelaunayTessellation2D{T},
      trig::DelaunayTriangle{T}) where {T<:AbstractPoint2D} = tess._trigs[getneighbor(trig,2)]
movec(tess::DelaunayTessellation2D{T},
      trig::DelaunayTriangle{T}) where {T<:AbstractPoint2D} = tess._trigs[getneighbor(trig,3)]
