"""
    DelaunayTesselation(points)

Returns a Delaunay Tesselation of the points in `points`. 
The result is convex. The points in `points` are not required 
to be in [1,2]x[1,2]
"""
function DelaunayTessellation( points::Array{T,1} ) where T<:AbstractPoint2D
  scaledPoints, ranges = scaleShiftPoints( points )
  scaledTess = DelaunayTessellation2D{T}(length( points ))
  push!( scaledTess, scaledPoints )
  tess = expand( scaledTess, ranges )
  return tess
end

# convert the tessellation back to original scale after the tessellation
function expand( tess::DelaunayTessellation2D{T},ranges) where T<:AbstractPoint2D
  xmin = ranges[1]
  ymin = ranges[3]
  scale = max( ranges[4] - ranges[3], ranges[2] - ranges[1] ) / 0.98
  offset = 1.01
  for i in 1:length(tess._trigs)
    #tess._trigs[i]._a = Point2D( ( tess._trigs[i]._a._x - offset ) * scale + xmin, ( tess._trigs[i]._a._y - offset ) * scale + ymin )
    tess._trigs[i]._a = expand([tess._trigs[i]._a],ranges)[1]
    tess._trigs[i]._b = expand([tess._trigs[i]._b],ranges)[1]
    tess._trigs[i]._c = expand([tess._trigs[i]._c],ranges)[1]
    #tess._trigs[i]._b = Point2D( ( tess._trigs[i]._b._x - offset ) * scale + xmin, ( tess._trigs[i]._b._y - offset ) * scale + ymin )
    #tess._trigs[i]._c = Point2D( ( tess._trigs[i]._c._x - offset ) * scale + xmin, ( tess._trigs[i]._c._y - offset ) * scale + ymin )
    tess._trigs[i]._bx   = tess._trigs[i]._bx  * scale
    tess._trigs[i]._by   = tess._trigs[i]._by  * scale
    tess._trigs[i]._cx   = tess._trigs[i]._cx  * scale
    tess._trigs[i]._cy   = tess._trigs[i]._cy  * scale
    tess._trigs[i]._px   = tess._trigs[i]._px  * scale ^ 3
    tess._trigs[i]._py   = tess._trigs[i]._py  * scale ^ 3
    tess._trigs[i]._pr2  = tess._trigs[i]._pr2 * scale ^ 2
  end
  tess._ranges = ranges
  return tess
end


# in order to reduce computation, if you are interested in the edges only, you can also
# apply the expand function to the points of the edges directly.
# if this is what you want, you can (see VoronoiDelaunayExtensions.jl)
# 1. apply scaleShiftPoints to your point set
# 2. do the tessellation and push the scaled and shifted point set
# 3. get the edges with the delaunayedges function and
# 4. expand the end points of the edges with the expand function



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


movea(tess::DelaunayTessellation2D{T},
      trig::DelaunayTriangle{T}) where {T<:AbstractPoint2D} = tess._trigs[trig._neighbour_a]
moveb(tess::DelaunayTessellation2D{T},
      trig::DelaunayTriangle{T}) where {T<:AbstractPoint2D} = tess._trigs[trig._neighbour_b]
movec(tess::DelaunayTessellation2D{T},
      trig::DelaunayTriangle{T}) where {T<:AbstractPoint2D} = tess._trigs[trig._neighbour_c]

#Finds out in which triangle the point was, and then "splits" the triangle in question into three sub-triangles
#This function does *not* enforce the delaunay property!
function _pushunfixed!(tess::DelaunayTessellation2D{T}, p::T) where T<:AbstractPoint2D
    #Find out which triangle the point `p` is in
    i = findindex(tess, p)
    
    #Add two new triangles at the end of the list
    ltrigs1::Int64 = tess._last_trig_index+1
    ltrigs2::Int64 = tess._last_trig_index+2
    tess._last_trig_index += 2

    #Modify the current triangle to be one of the three new triangles 
    #The point being added is vertex "a"
    @inbounds t1 = tess._trigs[i]
    old_t1_a = geta(t1)
    seta(t1, p)
    old_t1_b = t1._neighbour_b
    old_t1_c = t1._neighbour_c
    t1._neighbour_b = ltrigs1
    t1._neighbour_c = ltrigs2

    #Modify one of the new triangles to be the second of the new triangles
    #The point being added is vertex "b"
    @inbounds t2 = tess._trigs[ltrigs1]
    setabc(t2, old_t1_a, p, getc(t1))
    t2._neighbour_a = i
    t2._neighbour_b = old_t1_b
    t2._neighbour_c = ltrigs2

    #Modify one of the new triangles to be the third of the new triangles
    #The point being added is vertex "c"
    @inbounds t3 = tess._trigs[ltrigs2]
    setabc(t3, old_t1_a, getb(t1), p)
    t3._neighbour_a = i
    t3._neighbour_b = ltrigs1
    t3._neighbour_c = old_t1_c

    #Fix the neighbourhood relations of the previous neighbors
    @inbounds nt2 = tess._trigs[t2._neighbour_b]
    if nt2._neighbour_a==i
        nt2._neighbour_a = ltrigs1
    elseif nt2._neighbour_b==i
        nt2._neighbour_b = ltrigs1
    else
        nt2._neighbour_c = ltrigs1
    end

    @inbounds nt3 = tess._trigs[t3._neighbour_c]
    if nt3._neighbour_a==i
        nt3._neighbour_a = ltrigs2
    elseif nt3._neighbour_b==i
        nt3._neighbour_b = ltrigs2
    else
        nt3._neighbour_c = ltrigs2
    end


    return i
end

function _flipa!(tess::DelaunayTessellation2D, ix1::Int64, ix2::Int64)
    @inbounds ot1 = tess._trigs[ix1]
    @inbounds ot2 = tess._trigs[ix2]
    if ot2._neighbour_a == ix1
        _flipaa!(tess, ix1, ix2, ot1, ot2)
    elseif ot2._neighbour_b == ix1
        _flipab!(tess, ix1, ix2, ot1, ot2)
    else
        _flipac!(tess, ix1, ix2, ot1, ot2)
    end
end

function _endflipa!(tess::DelaunayTessellation2D,
                    ix1::Int64, ix2::Int64,
                    ot1::DelaunayTriangle, ot2::DelaunayTriangle)
    @inbounds n1 = tess._trigs[ot1._neighbour_a]
    if n1._neighbour_a==ix2
        n1._neighbour_a = ix1
    elseif n1._neighbour_b==ix2
        n1._neighbour_b = ix1
    else
        n1._neighbour_c = ix1
    end
    @inbounds n2 = tess._trigs[ot2._neighbour_c]
    if n2._neighbour_a==ix1
        n2._neighbour_a = ix2
    elseif n2._neighbour_b==ix1
        n2._neighbour_b = ix2
    else
        n2._neighbour_c = ix2
    end
end

function _flipaa!(tess::DelaunayTessellation2D{T},
                  ix1::Int64, ix2::Int64,
                  ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where T<:AbstractPoint2D
    old_ot1_geom_b = getb(ot1)
    setb(ot1, geta(ot2))
    ot1._neighbour_a = ot2._neighbour_c
    old_ot1_neighbour_c = ot1._neighbour_c
    ot1._neighbour_c = ix2

    newc = geta(ot2)
    setabc(ot2, geta(ot1), old_ot1_geom_b, newc)
    ot2._neighbour_a = ot2._neighbour_b
    ot2._neighbour_b = ix1
    ot2._neighbour_c = old_ot1_neighbour_c

    _endflipa!(tess,ix1,ix2,ot1,ot2)
end

function _flipab!(tess::DelaunayTessellation2D{T},
                  ix1::Int64, ix2::Int64,
                  ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where T<:AbstractPoint2D
    old_ot1_geom_b = getb(ot1)
    setb(ot1, getb(ot2))
    ot1._neighbour_a = ot2._neighbour_a
    old_ot1_neighbour_c = ot1._neighbour_c
    ot1._neighbour_c = ix2

    newc = getb(ot2)
    setabc(ot2, geta(ot1), old_ot1_geom_b, newc)
    ot2._neighbour_a = ot2._neighbour_c
    ot2._neighbour_b = ix1
    ot2._neighbour_c = old_ot1_neighbour_c

    _endflipa!(tess,ix1,ix2,ot1,ot2)
end

function _flipac!(tess::DelaunayTessellation2D{T},
                  ix1::Int64, ix2::Int64,
                  ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where T<:AbstractPoint2D
    old_ot1_geom_b = getb(ot1)
    setb(ot1, getc(ot2))
    ot1._neighbour_a = ot2._neighbour_b
    old_ot1_neighbour_c = ot1._neighbour_c
    ot1._neighbour_c = ix2

    setab(ot2, geta(ot1), old_ot1_geom_b)
    ot2._neighbour_b = ix1
    ot2._neighbour_c = old_ot1_neighbour_c

    _endflipa!(tess,ix1,ix2,ot1,ot2)
end

########################

function _flipb!(tess::DelaunayTessellation2D{T},
                 ix1::Int64, ix2::Int64) where T<:AbstractPoint2D
    @inbounds ot1 = tess._trigs[ix1]
    @inbounds ot2 = tess._trigs[ix2]
    if ot2._neighbour_a == ix1
        _flipba!(tess, ix1, ix2, ot1, ot2)
    elseif ot2._neighbour_b == ix1
        _flipbb!(tess, ix1, ix2, ot1, ot2)
    else
        _flipbc!(tess, ix1, ix2, ot1, ot2)
    end
end

function _endflipb!(tess::DelaunayTessellation2D{T},
                    ix1::Int64, ix2::Int64,
                    ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where T<:AbstractPoint2D
    @inbounds n1 = tess._trigs[ot1._neighbour_b]
    if n1._neighbour_a==ix2
        n1._neighbour_a = ix1
    elseif n1._neighbour_b==ix2
        n1._neighbour_b = ix1
    else
        n1._neighbour_c = ix1
    end
    @inbounds n2 = tess._trigs[ot2._neighbour_a]
    if n2._neighbour_a==ix1
        n2._neighbour_a = ix2
    elseif n2._neighbour_b==ix1
        n2._neighbour_b = ix2
    else
        n2._neighbour_c = ix2
    end
end
function _flipedge!(tess::DelaunayTessellation2D{T}, point1::Int64,point2::Int64,point3::Int64
                  ix1::Int64, ix2::Int64,
                  ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where T<:AbstractPoint2D


    old_ot1_geom_point3 = getpoint(ot1,point3)
    setpoint!(ot1,point3, getpoint(ot2,point2))

    old_ot1_neighbor_point2 = getneighbor(ot1,point2)
    setneighbor!(ot1,point2,ix2)
    setneighbor!(ot1,point1,getneighbor(ot2,point3))

    set2points!(ot2,point1,point3,getpoint(ot1,point1), old_ot1_geom_point3)

    setneighbor!(ot2, point2, old_ot1_neighbor_point2)
    setneighbor!(ot2, point3, ix1)

    _endflippoint!(tess,ix1,ix2,ot1,ot2)
end

 


#Is flipedge!(tess,2,1,3,_)
function _flipba!(tess::DelaunayTessellation2D{T},
                  ix1::Int64, ix2::Int64,
                  ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where T<:AbstractPoint2D
    old_ot1_geom_c = getc(ot1)
    setc(ot1, geta(ot2))
    old_ot1_neighbour_a = ot1._neighbour_a
    ot1._neighbour_a = ix2
    ot1._neighbour_b = ot2._neighbour_c

    setbc(ot2, getb(ot1), old_ot1_geom_c)
    ot2._neighbour_a = old_ot1_neighbour_a
    ot2._neighbour_c = ix1

    _endflipb!(tess,ix1,ix2,ot1,ot2)
end

function _flipbb!(tess::DelaunayTessellation2D{T},
                  ix1::Int64, ix2::Int64,
                  ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where T<:AbstractPoint2D
    old_ot1_geom_c = getc(ot1)
    setc(ot1, getb(ot2))
    old_ot1_neighbour_a = ot1._neighbour_a
    ot1._neighbour_a = ix2
    ot1._neighbour_b = ot2._neighbour_a

    newa = getb(ot2)
    setabc(ot2, newa, getb(ot1), old_ot1_geom_c)
    ot2._neighbour_a = old_ot1_neighbour_a
    ot2._neighbour_b = ot2._neighbour_c
    ot2._neighbour_c = ix1

    _endflipb!(tess,ix1,ix2,ot1,ot2)
end

function _flipbc!(tess::DelaunayTessellation2D{T},
                  ix1::Int64, ix2::Int64,
                  ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where T<:AbstractPoint2D
    old_ot1_geom_c = getc(ot1)
    setc(ot1, getc(ot2))
    old_ot1_neighbour_a = ot1._neighbour_a
    ot1._neighbour_a = ix2
    ot1._neighbour_b = ot2._neighbour_b

    newa = getc(ot2)
    setabc(ot2, newa, getb(ot1), old_ot1_geom_c)
    ot2._neighbour_b = ot2._neighbour_a
    ot2._neighbour_a = old_ot1_neighbour_a
    ot2._neighbour_c = ix1

    _endflipb!(tess,ix1,ix2,ot1,ot2)
end

########################

function _flipc!(tess::DelaunayTessellation2D{T},
                 ix1::Int64, ix2::Int64) where T<:AbstractPoint2D
    @inbounds ot1 = tess._trigs[ix1]
    @inbounds ot2 = tess._trigs[ix2]
    if ot2._neighbour_a == ix1
        _flipca!(tess, ix1, ix2, ot1, ot2)
    elseif ot2._neighbour_b == ix1
        _flipcb!(tess, ix1, ix2, ot1, ot2)
    else
        _flipcc!(tess, ix1, ix2, ot1, ot2)
    end
end

function _endflipc!(tess::DelaunayTessellation2D{T},
                    ix1::Int64, ix2::Int64,
                    ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where T<:AbstractPoint2D
    @inbounds n1 = tess._trigs[ot1._neighbour_c]
    if n1._neighbour_a==ix2
        n1._neighbour_a = ix1
    elseif n1._neighbour_b==ix2
        n1._neighbour_b = ix1
    else
        n1._neighbour_c = ix1
    end
    @inbounds n2 = tess._trigs[ot2._neighbour_b]
    if n2._neighbour_a==ix1
        n2._neighbour_a = ix2
    elseif n2._neighbour_b==ix1
        n2._neighbour_b = ix2
    else
        n2._neighbour_c = ix2
    end
end

function _flipca!(tess::DelaunayTessellation2D{T},
                  ix1::Int64, ix2::Int64,
                  ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where T<:AbstractPoint2D
    old_ot1_geom_a = geta(ot1)
    seta(ot1, geta(ot2))
    old_ot1_neighbour_b = ot1._neighbour_b
    ot1._neighbour_b = ix2
    ot1._neighbour_c = ot2._neighbour_c

    newb = geta(ot2)
    setabc(ot2, old_ot1_geom_a, newb, getc(ot1))
    ot2._neighbour_a = ix1
    ot2._neighbour_c = ot2._neighbour_b
    ot2._neighbour_b = old_ot1_neighbour_b

    _endflipc!(tess,ix1,ix2,ot1,ot2)
end

function _flipcb!(tess::DelaunayTessellation2D{T},
                  ix1::Int64, ix2::Int64,
                  ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where T<:AbstractPoint2D
    old_ot1_geom_a = geta(ot1)
    seta(ot1, getb(ot2))
    old_ot1_neighbour_b = ot1._neighbour_b
    ot1._neighbour_b = ix2
    ot1._neighbour_c = ot2._neighbour_a

    setac(ot2, old_ot1_geom_a, getc(ot1))
    ot2._neighbour_a = ix1
    ot2._neighbour_b = old_ot1_neighbour_b

    _endflipc!(tess,ix1,ix2,ot1,ot2)
end

function _flipcc!(tess::DelaunayTessellation2D{T},
                  ix1::Int64, ix2::Int64,
                  ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where T<:AbstractPoint2D
    old_ot1_geom_a = geta(ot1)
    seta(ot1, getc(ot2))
    old_ot1_neighbour_b = ot1._neighbour_b
    ot1._neighbour_b = ix2
    ot1._neighbour_c = ot2._neighbour_b

    newb = getc(ot2)
    setabc(ot2, old_ot1_geom_a, newb, getc(ot1))
    ot2._neighbour_c = ot2._neighbour_a
    ot2._neighbour_a = ix1
    ot2._neighbour_b = old_ot1_neighbour_b

    _endflipc!(tess,ix1,ix2,ot1,ot2)
end

"""
    _restoredelaunayhood!(tess,ix_trig)

Restores the Delaunay property by flipping edges. 
Here `ix_trig` is the index of the triangle whose `a` vertex
is the new point that was added (see `_pushunfixed`) and whose `b`
and `c` neighbors are the other "new" triangles.
"""
function _restoredelaunayhood!(tess::DelaunayTessellation2D{T},
                               ix_trig::Int64) where T<:AbstractPoint2D
    #This is the point that was newly added
    @inbounds center_pt = geta(tess._trigs[ix_trig])

    # `A` - edge
    push!(tess._edges_to_check, ix_trig)
    while length(tess._edges_to_check) > 0
        #Find the corresponding `a` vertex
        @inbounds trix = tess._edges_to_check[end]
        @inbounds tr_i = tess._trigs[trix]
        @inbounds nb_a = tr_i._neighbour_a
        #Check that we are not at the boundary
        if nb_a > 1
            @inbounds tr_f = tess._trigs[nb_a]
            if incircle(tr_f, center_pt) > 0
                _flipa!(tess, trix, nb_a)
                push!(tess._edges_to_check, nb_a)
                continue
            end
        end
        pop!(tess._edges_to_check)
    end

    # `B` - edge
    push!(tess._edges_to_check, tess._last_trig_index-1)
    while length(tess._edges_to_check) > 0
        @inbounds trix = tess._edges_to_check[end]
        @inbounds tr_i = tess._trigs[trix]
        @inbounds nb_b = tr_i._neighbour_b
        if nb_b > 1
            @inbounds tr_f = tess._trigs[nb_b]
            if incircle(tr_f, center_pt) > 0
                _flipb!(tess, trix, nb_b)
                push!(tess._edges_to_check, nb_b)
                continue
            end
        end
        pop!(tess._edges_to_check)
    end

    # `C` - edge
    push!(tess._edges_to_check, tess._last_trig_index)
    while length(tess._edges_to_check) > 0
        @inbounds trix = tess._edges_to_check[end]
        @inbounds tr_i = tess._trigs[trix]
        @inbounds nb_c = tr_i._neighbour_c
        if nb_c > 1
            @inbounds tr_f = tess._trigs[nb_c]
            if incircle(tr_f, center_pt) > 0
                _flipc!(tess, trix, nb_c)
                push!(tess._edges_to_check, nb_c)
                continue
            end
        end
        pop!(tess._edges_to_check)
    end
end

# push a single point. Grows tessellation as needed
function push!(tess::DelaunayTessellation2D{T}, p::T) where T<:AbstractPoint2D
    tess._total_points_added += 1
    sizefit_at_least(tess, tess._total_points_added)
    i = _pushunfixed!(tess, p)
    _restoredelaunayhood!(tess, i)
end

# push an array in given order
function _pushunsorted!(tess::DelaunayTessellation2D{T}, a::Array{T, 1}) where T<:AbstractPoint2D
    sizehint!(tess, length(a))
    for p in a
        push!(tess, p)
    end
end

# push an array but sort it first for better performance
function push!(tess::DelaunayTessellation2D{T}, a::Array{T, 1}) where T<:AbstractPoint2D
    shuffle!(a)
    mssort!(a)
    _pushunsorted!(tess, a)
end

