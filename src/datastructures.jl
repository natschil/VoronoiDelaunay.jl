#Struct for representing a single triangle

mutable struct DelaunayTriangle{T<:AbstractPoint2D} <: AbstractNegativelyOrientedTriangle
    #The three points of the triangle
    points::Array{T} 

    #TODO: Figure out what these are....
    _bx::Float64; _by::Float64
    _cx::Float64; _cy::Float64

    _px::Float64; _py::Float64
    _pr2::Float64

    #Indexes of neighboring points
    neighbors::Array{Int64}

    function DelaunayTriangle{T}(pa::T, pb::T, pc::T,
                                 na::Int64, nb::Int64, nc::Int64) where {T<:AbstractPoint2D}
        t = new([pa, pb, pc], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, [na, nb, nc])
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
        new([pa, pb, pc], bx, by, cx, cy, px, py, pr2, [na, nb, nc])
    end
end

GP.geta(tri::DelaunayTriangle{T}) where T<:AbstractPoint2D = tri.points[1]
GP.getb(tri::DelaunayTriangle{T}) where T<:AbstractPoint2D = tri.points[2]
GP.getc(tri::DelaunayTriangle{T}) where T<:AbstractPoint2D = tri.points[3]
GP.seta(tri::DelaunayTriangle{T},val::T) where T<:AbstractPoint2D = (tri.points[1] = val; clean!(tri))
GP.setb(tri::DelaunayTriangle{T},val::T) where T<:AbstractPoint2D = (tri.points[2] = val;clean!(tri))
GP.setc(tri::DelaunayTriangle{T},val::T) where T<:AbstractPoint2D = (tri.points[3] = val;clean!(tri))
function GP.setabc(tri::DelaunayTriangle{T},vala::T,valb::T,valc::T) where T<:AbstractPoint2D
    tri.points[1] = vala
    tri.points[2] = valb
    tri.points[3] = valc
    clean!(tri)
    return tri
end


getpoint(tri::DelaunayTriangle{T}, index::Int64) where T<:AbstractPoint2D = tri.points[index]
setpoint!(tri::DelaunayTriangle{T}, index::Int64,val::T) where T<:AbstractPoint2D = (tri.points[index] = val; clean!(tri))
getneighbor(tri::DelaunayTriangle{T}, index::Int64) where T<:AbstractPoint2D = tri.neighbors[index]
setneighbor!(tri::DelaunayTriangle{T}, index::Int64,val::Int64) where T<:AbstractPoint2D = (tri.neighbors[index] = val)

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

function isexternal( t::DelaunayTriangle{T}, ranges::NTuple{4,Float64} ) where T<:AbstractPoint2D
    getx(geta(t)) < ranges[1] || getx(geta(t)) > ranges[2] ||
    getx(getb(t)) < ranges[1] || getx(getb(t)) > ranges[2] ||
    getx(getc(t)) < ranges[1] || getx(getc(t)) > ranges[2]
end

mutable struct DelaunayTessellation2D{T<:AbstractPoint2D}
    _trigs::Vector{DelaunayTriangle{T}}
    _last_trig_index::Int64 #How many of the elements in _trigs are actually being used
    _prev_trig_index::Int64 #The triangle that the previous search ended in
    _edges_to_check::Vector{Int64}
    _total_points_added::Int64
    _ranges::NTuple{4,Float64}

    function DelaunayTessellation2D{T}(n::Int64 = 100) where T
        a = T(GeometricalPredicates.min_coord, GeometricalPredicates.min_coord)
        b = T(GeometricalPredicates.min_coord, GeometricalPredicates.max_coord)
        c = T(GeometricalPredicates.max_coord, GeometricalPredicates.min_coord)
        d = T(GeometricalPredicates.max_coord, GeometricalPredicates.max_coord)
        t1 = DelaunayTriangle{T}(d,c,b, 2,1,1)
        t2 = DelaunayTriangle{T}(a,b,c, 3,1,1)
        t3 = DelaunayTriangle{T}(d,c,b, 2,1,1)
        _trigs = DelaunayTriangle{T}[t1, t2, t3]
        _ranges=(min_coord,max_coord,min_coord,max_coord)
        t = new(_trigs, 3,1, Int64[], 0,_ranges)
        sizehint!(t._edges_to_check, 1000)
        sizehint!(t, n)
        return t
    end
end
DelaunayTessellation2D(n::Int64) = DelaunayTessellation2D{Point2D}(n)
DelaunayTessellation2D(n::Int64, ::T) where {T<:AbstractPoint2D} = DelaunayTessellation2D{T}(n)
DelaunayTessellation(n::Int64=100) = DelaunayTessellation2D(n)

function sizehint!(t::DelaunayTessellation2D{T}, n::Int64) where T<:AbstractPoint2D
    required_total_size = 2n + 10
    required_total_size <= length(t._trigs) && return
    sizehint!(t._trigs, required_total_size)
    while length(t._trigs) < required_total_size
        push!(t._trigs, copy(t._trigs[end]))
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
    ranges = t._ranges
    visited = zeros(Bool, t._last_trig_index)
    function delaunayiterator(c::Channel)
        @inbounds for ix in 2:t._last_trig_index
            tr = t._trigs[ix]
            !include_external && isexternal( tr, ranges ) && continue
            visited[ix] && continue
            visited[ix] = true
            ix_na = getneighbor(tr,1)
            if !visited[ix_na]
                put!(c, DelaunayEdge(getb(tr), getc(tr)))
            end
            ix_nb = getneighbor(tr,2)
            if !visited[ix_nb]
                put!(c, DelaunayEdge(geta(tr), getc(tr)))
            end
            ix_nc = getneighbor(tr,3)
            if !visited[ix_nc]
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
    ranges = t._ranges
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

function findindex(tess::DelaunayTessellation2D{T}, p::T) where T<:AbstractPoint2D
    #i::Int64 = tess._prev_trig_index
    i::Int64 = tess._last_trig_index
    counter::Int64 = 0
    ntriangles = length(tess._trigs)
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

locate(t::DelaunayTessellation2D{T}, p::T) where {T<:AbstractPoint2D} = t._trigs[findindex(t, p)]


movea(tess::DelaunayTessellation2D{T},
      trig::DelaunayTriangle{T}) where {T<:AbstractPoint2D} = tess._trigs[getneighbor(trig,1)]
moveb(tess::DelaunayTessellation2D{T},
      trig::DelaunayTriangle{T}) where {T<:AbstractPoint2D} = tess._trigs[getneighbor(trig,2)]
movec(tess::DelaunayTessellation2D{T},
      trig::DelaunayTriangle{T}) where {T<:AbstractPoint2D} = tess._trigs[getneighbor(trig,3)]




