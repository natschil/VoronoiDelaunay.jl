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


function _visual_restore_delaunayhood!(tess::DelaunayTessellation2D{T},
                               ix_trig::Int64) where T<:AbstractPoint2D
    #This is the point that was newly added
    @inbounds center_pt = geta(tess._trigs[ix_trig])

    to_add = [ix_trig,tess._last_trig_index - 1, tess._last_trig_index]

    res = [getplotxy(delaunayedges(tess,true))]
    for j in [1,2,3]
        push!(tess._edges_to_check, to_add[j])
        while length(tess._edges_to_check) > 0
            #Find the corresponding `j` vertex

            #This is the triangle just added
            @inbounds trix = tess._edges_to_check[end]
            @inbounds tr_i = tess._trigs[trix]

            #This is it's j-th neighbor
            @inbounds nb_j = getneighbor(tr_i,j)

            #Check that we are not at the boundary
            if nb_j > 1
                @inbounds tr_f = tess._trigs[nb_j]
                if _should_flip_edge(trix,tr_i,tr_f,center_pt,j,tess._convex)
                    _flip!(tess,j,trix,nb_j)
                    getplotxy
                    push!(res,getplotxy(delaunayedges(tess,true)))
                    push!(tess._edges_to_check, nb_j)
                    continue
                end
            end
            pop!(tess._edges_to_check)
        end
    end

    _clean_nonexternal_indexes!(tess)
    return res
end


function visual_push!(tess::DelaunayTessellation2D{T}, p::T,clean=true) where T<:AbstractPoint2D
    tess._total_points_added += 1
    sizefit_at_least(tess, tess._total_points_added)

    if tess._convex
        new_pt = rescale(p,tess._scaling)
        if ! _in_ref_tri(new_pt)
            invscale = tess._invscaling
            new_scaling = get_new_scaling(invscale,p,p)
            rescale!(tess,new_scaling)
            new_pt = rescale(p,tess._scaling)

        end
    else
        new_pt = p
        if !in_base_square_domain(new_pt)
            error("Trying to add point that is out of bounds")
        end
    end

    i = _pushunfixed!(tess, new_pt,false)
    println(new_pt)
    res = _visual_restore_delaunayhood!(tess, i)
    if clean
        _clean_nonexternal_indexes!(tess)
    end
    return res
end


function push!(tess::DelaunayTessellation2D{T}, a::Array{T, 1},clean=true) where T<:AbstractPoint2D
    shuffle!(a)
    mssort!(a)
    res = []
    if tess._convex
        maxx = maximum([getx(p) for p in a])
        minx = minimum([getx(p) for p in a])
        maxy = maximum([gety(p) for p in a])
        miny = minimum([gety(p) for p in a])
        p1 = Point2D(minx,miny)
        p2 = Point2D(maxx,maxy)
        if ! _in_ref_tri(rescale(p1,tess._scaling)) || ! _in_ref_tri(rescale(p2,tess._scaling))
            invscale = tess._invscaling
            new_scaling = get_new_scaling(invscale,p1,p2)
            rescale!(tess,new_scaling)
        end
    end
    _pushunsorted!(tess, a,clean)
end
