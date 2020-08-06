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

"""
    _flip!(tess,point_added,ix1,ix2,ot1,ot2)

Flips a single edge in the triangulation. 
Here `point_added` is the index (within `ix1`, so one of 1,2,3) of the 
new point that was added. 
Here `ix1` is the triangle that has this point as its vertex 
and the edge to be flipped opposite to it.
Here `ix2` is the triangle on the other side of this edge.

"""
function _flip!(tess::DelaunayTessellation2D{T},point_added::Int64,
                  ix1::Int64, ix2::Int64,
                  ) where T<:AbstractPoint2D


    @inbounds ot1 = tess._trigs[ix1]
    @inbounds ot2 = tess._trigs[ix2]

    #We want ot2.getneighbor(point_other) == ix1,
    #i.e. to get the local index of ix1
    point_other::Int64 = 0
    if getneighbor(ot2,1) == ix1
        point_other = 1
    elseif getneighbor(ot2,2) ==  ix1
        point_other = 2
    else
        point_other = 3
    end

    #These are the local indices of the other two points in ix1
    point_added_p = mod1(point_added+1,3)
    point_added_pp = mod1(point_added+2,3)

    #These are the local indices of the other two points in ix2
    point_other_p = mod1(point_other+1,3)
    point_other_pp = mod1(point_other+2,3)

    #Flip the edge
    
    #First we fix ot1
    #We can keep `point_added` in the triangle.
    #We can also safely keep `point_added_pp`
    #This means we must replace `point_added_p` (with the point that is accross
    #the edge i.e. `point_other`
    old_ot1_geom_point_added_p = getpoint(ot1,point_added_p)
    setpoint!(ot1,point_added_p, getpoint(ot2,point_other))
    setneighbor!(ot1,point_added, getneighbor(ot2,point_other_pp))
    old_ot1_neighbor_point_added_pp = getneighbor!(ot1,point_added_pp)
    setneighbor!(ot1,point_added_pp,ix2)


    #Fix ot2
    new_pp = getpoint(ot2,point_other)
    setpoint!(ot1,point_other, getpoint(ot1,point_added) )
    setpoint!(ot1,point_other_p, old_ot1_geom_point_added_p)
    setpoint!(ot1,point_other_pp, new_pp )

    setneighbor!(ot2,point_other, getneighbor!(ot2,point_other_p))
    setneighbor!(ot2,point_other_p,ix1)
    setneighbor!(ot2,point_other_pp,old_ot1_neighbor_point_added_pp)


    #Fix neighbors of triangles that were not fixed

    @inbounds n1 = tess._trigs[getneighbor(ot1,point_added)]
    if getneighbor(n1,1) ==ix2
        setneighbor!(n1,1,ix1)
    elseif getneighbor(n1,2) == ix2
        setneighbor!(n1,2,ix1)
    else
        setneighbor!(n1,3,ix1)
    end
    @inbounds n2 = tess._trigs[getneighbor(ot2,point_other_pp)]
    if getneighbor(n2,1) ==ix1
        setneighbor!(n2,1,ix2)
    elseif getneighbor(n2,2) == ix1
        setneighbor!(n2,2,ix2)
    else
        setneighbor!(n2,3,ix2)
    end

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
    for j in [1,2,3]
        push!(tess._edges_to_check, ix_trig)
        while length(tess._edges_to_check) > 0
            #Find the corresponding `a` vertex
            @inbounds trix = tess._edges_to_check[end]
            @inbounds tr_i = tess._trigs[trix]
            @inbounds nb_j = getneighbor(tr_i,j)
            #Check that we are not at the boundary
            if nb_a > 1
                @inbounds tr_f = tess._trigs[nb_j]
                if incircle(tr_f, center_pt) > 0
                    _flip!(tess,j,trix,nb_j)
                    push!(tess._edges_to_check, nb_j)
                    continue
                end
            end
            pop!(tess._edges_to_check)
        end
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

