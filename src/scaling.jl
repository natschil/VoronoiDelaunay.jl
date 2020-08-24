"""
end
   rescale(point,scalefactor)

Creates a new point with x coordinates `(getx(point) - scalefactor[1])/scalefactor[2] ` as x-coordinate and y coordinate
given by `(gety(point) - scalefactor[3])/scalefactor[4]`
"""
function rescale( point::Point2D, scalefactor::NTuple{4,Float64})
    return Point2D((getx(point) - scalefactor[1])/scalefactor[2], (gety(point) - scalefactor[3])/scalefactor[4])
end

"""
    rescale(points,scalefactor)
Applies `rescale` pointwise to the array of points
"""

function rescale( points::AbstractArray{T}, scalefactor::NTuple{4,Float64}) where T <: AbstractPoint2D
    return map(rescale,points)
end

"""
   invertscaling(scalefactor)
Returns a 4-tuple `res` so that `rescale(p,res)` is the inverse to `rescale(p,scalefactor)`
"""
function invertscaling(scalefactor::NTuple{4,Float64})
    return (-scalefactor[1]/scalefactor[2], 1/scalefactor[2], -scalefactor[3]/scalefactor[4],1/scalefactor[4])
end

"""
    composescaling(s, t)
Takes 4-tuples `s` and `t` and returns a 4-tuple `res` that satisfies
rescale(p,`res`) = rescale(s,rescale(t,p))`. 
"""
function composescaling(s2::NTuple{4,Float64},s1::NTuple{4,Float64})
    return (s1[1] + s2[1]*s1[2], s1[2]*s2[2], s1[3] + s2[3]*s1[4], s1[4]*s2[4])
end

"""
    get_new_scaling(invscale_old,p1,p2)

Calculates a scaling so that invscale_old(S), p1 and p2 are all mapped to S, where S is the
square with lower left `min_coord + 1/4, min_coord`  and height/width 1/2.
Includes a small safety factor

"""
function get_new_scaling(invscale_old::NTuple{4,Float64}, p3::T, p4::T ) where T <: AbstractPoint2D
            p0 = rescale(Point2D(min_coord + 1/4 + eps(),min_coord),invscale_old)
            p1 = rescale(Point2D(min_coord + 3/4, min_coord + 1/2 ),invscale_old)
            maxx = max(max(getx(p1),getx(p3),getx(p4)))
            maxy = max(max(gety(p1),gety(p3)),gety(p4))

            minx = min(min(getx(p0),getx(p3)),getx(p4))
            miny = min(min(gety(p0),gety(p3)),gety(p4))

            widthx = maxx - minx
            widthy = maxy - miny
            
            newwidth = max(widthx,widthy)

            #Add a bit of a safety factor so we don't do this too often
            default_x = (min_coord + 1/4 )
            default_y = min_coord
            default_width = 1/2
            tosquare = invertscaling((default_x,default_width,default_y,default_width))

            new_scaling = composescaling(tosquare ,(minx - 0.1*newwidth,1.5*newwidth,miny - 0.1*newwidth,1.5*newwidth))
            return new_scaling
end



#=
function expand( points::Array{Point2D,1}, ranges::NTuple{4,Float64} )
  xmin = ranges[1]
  ymin = ranges[3]
  scale = max( ranges[4] - ranges[3], ranges[2] - ranges[1] ) / 0.98
  offset = 1.01
  scaledPoints = [ Point2D( ( p._x - offset ) * scale + xmin, ( p._y - offset ) * scale + ymin ) for p in points ]
end

function shrink( points::Array{Point2D,1}, ranges::NTuple{4,Float64} )
  h = ranges[4] - ranges[3]
  b = ranges[2] - ranges[1]
  offset = 1.01
  scale = 0.98 / max( h, b )
  scaledPoints = [ Point2D( ( p._x - ranges[1] ) * scale + offset, ( p._y - ranges[3] ) * scale + offset ) for p in points ]
end
=#



#=== auxiliary functions and structs ===#

# convex hull of point set
#=
function quickHull(points::Array{T,1}) where T<:AbstractPoint2D
    return quickHull(Point2D.(points))
end

function quickHull( points::Array{Point2D,1} )
  ConvexHull = Array{Point2D,1}(undef,0)
  A = points[ argmin( getx.( points ) ) ]   # leftmost point
  B = points[ argmax( getx.( points ) ) ]   # rightmost point
  push!( ConvexHull, A )
  push!( ConvexHull, B )
  pr = Array{Point2D,1}(undef,0)            # points to the right of AB
  pl = Array{Point2D,1}(undef,0)            # points to the left of AB
  ( pl, pr ) = dividePointSet( Line2D( A, B ), setdiff( points, ConvexHull ) )
  findHull!( ConvexHull, pr, A, B )         # divide-and-conquer approach
  findHull!( ConvexHull, pl, B, A )
  return ConvexHull
end

function findHull!( ConvexHull::Array{T,1}, points::Array{T,1}, A::T, B::T ) where T<:AbstractPoint2D
  if isempty( points )
    return
  end
  C = findFarthestPoint( A, B, points )
  pos = findfirst( x -> x == A, ConvexHull ) + 1
  insert!( ConvexHull, pos, C )
  pAC = dividePointSet( Line2D(A,C), setdiff(points,[C]) )[2]
  pCB = dividePointSet( Line2D(C,B), setdiff(points,[C]) )[2]
  findHull!( ConvexHull, pAC, A, C )
  findHull!( ConvexHull, pCB, C, B )
end

function dividePointSet( l::Line2D, points::Array{T,1} ) where T <: AbstractPoint2D
  pr = Array{T,1}(undef,0)
  pl = Array{T,1}(undef,0)
  for i in 1:length(points)
    if orientation(l,points[i]) == -1
      push!(pr,points[i])
    else
      push!(pl,points[i])
    end
  end
  return (pl, pr)
end

function findFarthestPoint(a::T,b::T,points::Array{T,1}) where T<:AbstractPoint2D
  distances = Array{Float64,1}(undef,size(points,1))
  for i in 1:length(distances)
    distances[i] = pointLineDistance(a,b,points[i])
  end
  return points[argmax(distances)]
end

function pointLineDistance(a::T,b::T,c::T) where T<:AbstractPoint2D
  return abs( (gety(b)-gety(a))*getx(c) - (getx(b)-getx(a))*gety(c) + getx(b)*gety(a) - gety(b)*getx(a) ) / sqrt( (gety(b)-gety(a))^2 + (getx(b)-getx(a))^2 )
end



mutable struct Circle{T<:AbstractPoint2D}
  c::T
  r::Float64
end

function getcen(circ::Circle{T}) where T<:AbstractPoint2D
    return circ.c
end

function getcenx(circ::Circle{T}) where T<:AbstractPoint2D
    return getx( getcen( circ ) )
end
function getceny(circ::Circle{T}) where T<:AbstractPoint2D
    return gety( getcen( circ ) )
end
function getrad(circ::Circle{T}) where T<: AbstractPoint2D
    return circ.r
end

function circumcircleUnion( hull::Array{T,1}, points::Array{S,1} ) where {T<:AbstractPoint2D,S<:AbstractPoint2D}
  ccU = Array{Circle{Point2D},1}(undef,length(hull))
  hullcirc = copy(hull)
  push!( hullcirc, hull[1] )
  for i in 1:length( hull )
    cc2 = circumcircle( hullcirc[i], hullcirc[i+1] )
    ccU[i] = cc2
    for j in 1:length(points)
      if pointsDistance( points[j], cc2.c ) < cc2.r
        cc3 = circumcircle( hullcirc[i], hullcirc[i+1], points[j] )
        if getrad(cc3) > getrad(ccU[i])
          ccU[i] = cc3
        end
      end
    end
  end
  return ccU
end

function pointsDistance( p1::T, p2::S ) where {T<:AbstractPoint2D, S<:AbstractPoint2D}
  sqrt( ( getx(p1) - getx(p2) ) ^ 2 + ( gety(p1) - gety(p2) ) ^ 2 )
end

function circumcircle( a::T, b::S ) where {T<:AbstractPoint2D, S<:AbstractPoint2D}
  px = (getx(a) + getx(b)) / 2.0
  py = (gety(a) + gety(b)) / 2.0
  dx = (getx(b) - getx(a)) / 2.0
  dy = (gety(b) - gety(a)) / 2.0
  r = sqrt( dx^2 + dy^2 )
  return Circle{Point2D}( Point2D(px,py), r )
end

function circumcircle( a::T, b::S, c::U ) where {T<:AbstractPoint2D,S<:AbstractPoint2D, U<:AbstractPoint2D}
  la = getx(a)^2 + gety(a)^2
  lb = getx(b)^2 + gety(b)^2
  lc = getx(c)^2 + gety(c)^2
  xyy = getx(a) * (gety(b) - gety(c))
  yxx = gety(a) * (getx(b) - getx(c))
  xy = getx(b) * gety(c)
  yx = gety(b) * getx(c)
  z = 2 * ( xyy - yxx + xy - yx )
  px = ( la*(gety(b)-gety(c)) + lb*(gety(c)-gety(a)) + lc*(gety(a)-gety(b)) ) / z
  py = ( la*(getx(c)-getx(b)) + lb*(getx(a)-getx(c)) + lc*(getx(b)-getx(a)) ) / z
  r = sqrt( (px-getx(a))^2 + (py-gety(a))^2 )
  return Circle{Point2D}( Point2D(px,py), r )
end




function frameRanges( ccU::Array{Circle{T},1} ) where T<:AbstractPoint2D
  xmin = minimum( getcenx.(ccU) .- getrad.(ccU) )
  xmax = maximum( getcenx.(ccU) .+ getrad.(ccU) )
  ymin = minimum( getceny.(ccU) .- getrad.(ccU) )
  ymax = maximum( getceny.(ccU) .+ getrad.(ccU) )
  return xmin, xmax, ymin, ymax
end
=#
