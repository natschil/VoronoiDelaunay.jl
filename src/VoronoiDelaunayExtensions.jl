# scale the point set down such that the union of all circumcircles lies within the square frame.
# this does not only allow for arbitrary point distributions, but more importantly
# ensures convexity (and correctness) of the resulting Delaunay triangulation. 

function scaleShiftPoints( points::Array{Point2D,1} )
  # convex hull of point set to find the outer edges
  hull = quickHull( points )
  # union of the outer circumcircles
  ccu = circumcircleUnion( hull, points )
  # frame around the union of circumcircles
  ranges = frameRanges( ccu )
  # scale the point set down to the initial square
  scaledPoints = shrink( points, ranges )
  # do the tessellation on the scaled point set,
  # use the ranges to convert back to original scale afterwards
  return scaledPoints, ranges
end

# convert the edge points back to original scale after the tessellation

function expand( points::Array{Point2D,1}, ranges::NTuple{4,Float64} )
  scale = max( ranges[4] - ranges[3], ranges[2] - ranges[1] ) / 0.98
  offset = 1.01
  points = [ Point2D( ( p._x - ranges[1] ) * scale + offset, ( p._y - ranges[3] ) * scale + offset ) for p in points ]
end



#=== auxiliary functions and structs ===#

# convex hull of point set

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

function findHull!( ConvexHull::Array{Point2D,1}, points::Array{Point2D,1}, A::Point2D, B::Point2D )
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

function dividePointSet( l::Line2D, points::Array{Point2D,1} )
  pr = Array{Point2D,1}(undef,0)
  pl = Array{Point2D,1}(undef,0)
  for i in 1:length(points)
    if orientation(l,points[i]) == -1
      push!(pr,points[i])
    else
      push!(pl,points[i])
    end
  end
  return (pl, pr)
end

function findFarthestPoint(a::Point2D,b::Point2D,points::Array{Point2D,1})
  distances = Array{Float64,1}(undef,size(points,1))
  for i in 1:length(distances)
    distances[i] = pointLineDistance(a,b,points[i])
  end
  return points[argmax(distances)]
end

function pointLineDistance(a::Point2D,b::Point2D,c::Point2D)
  abs( (b._y-a._y)*c._x - (b._x-a._x)*c._y + b._x*a._y - b._y*a._x ) / sqrt( (b._y-a._y)^2 + (b._x-a._x)^2 )
end



mutable struct Circle
  c::Point2D
  r::Float64
end

function circumcircleUnion( hull::Array{Point2D,1}, points::Array{Point2D,1} )
  ccU = Array{Circle,1}(undef,length(hull))
  hullcirc = copy(hull)
  push!( hullcirc, hull[1] )
  for i in 1:length( hull )
    cc2 = circumcircle( hullcirc[i], hullcirc[i+1] )
    ccU[i] = cc2
    for j in 1:length(points)
      if pointsDistance( points[j], cc2.c ) < cc2.r
        cc3 = circumcircle( hullcirc[i], hullcirc[i+1], points[j] )
        if cc3.r > ccU[i].r
          ccU[i] = cc3
        end
      end
    end
  end
  return ccU
end

function pointsDistance( p1::Point2D, p2::Point2D )
  sqrt( ( p1._x - p2._x ) ^ 2 + ( p1._y - p2._y ) ^ 2 )
end

function circumcircle( a::Point2D, b::Point2D )
  px = (a._x + b._x) / 2.0
  py = (a._y + b._y) / 2.0
  dx = (b._x - a._x) / 2.0
  dy = (b._y - a._y) / 2.0
  r = sqrt( dx^2 + dy^2 )
  return Circle( Point2D(px,py), r )
end

function circumcircle( a::Point2D, b::Point2D, c::Point2D )
  la = a._x^2 + a._y^2
  lb = b._x^2 + b._y^2
  lc = c._x^2 + c._y^2
  xyy = a._x * (b._y - c._y)
  yxx = a._y * (b._x - c._x)
  xy = b._x * c._y
  yx = b._y * c._x
  z = 2 * ( xyy - yxx + xy - yx )
  px = ( la*(b._y-c._y) + lb*(c._y-a._y) + lc*(a._y-b._y) ) / z
  py = ( la*(c._x-b._x) + lb*(a._x-c._x) + lc*(b._x-a._x) ) / z
  r = sqrt( (px-a._x)^2 + (py-a._y)^2 )
  return Circle( Point2D(px,py), r )
end

getcen(circ::Circle)  = circ.c
getcenx(circ::Circle) = getx( getcen( circ ) )
getceny(circ::Circle) = gety( getcen( circ ) )
getrad(circ::Circle)  = circ.r



function frameRanges( ccU::Array{Circle,1} )
  xmin = minimum( getcenx.(ccU) .- getrad.(ccU) )
  xmax = maximum( getcenx.(ccU) .+ getrad.(ccU) )
  ymin = minimum( getceny.(ccU) .- getrad.(ccU) )
  ymax = maximum( getceny.(ccU) .+ getrad.(ccU) )  
  return xmin, xmax, ymin, ymax
end
    


function shrink( points::Array{Point2D,1}, ranges::NTuple{4,Float64} )
  h = ranges[4] - ranges[3]
  b = ranges[2] - ranges[1]
  offset = 1.01
  scale = 0.98 / max( h, b )
  points = [ Point2D( ( p._x - ranges[1] ) * scale + offset, ( p._y - ranges[3] ) * scale + offset ) for p in points ]
end
