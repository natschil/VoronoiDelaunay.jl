#Helper type for keeping track of point numbers
struct NumberedPoint2D <: AbstractPoint2D
    x::Float64
    y::Float64
    id::Int64
    NumberedPoint2D(x::Float64,y::Float64,k::Int64) = new(x,y,k)
    NumberedPoint2D(x::Float64,y::Float64) = new(x, y, 0)
    NumberedPoint2D(p::Point2D) = new(p.x, p.y, 0)
    NumberedPoint2D(p) = new(p[1], p[2], 0)
end
GP.Point(x::Real, y::Real, k::Int64) = NumberedPoint2D(x, y, k)
GP.Point2D(p::NumberedPoint2D) = GP.Point2D(p.x,p.y)
GP.gety(p::NumberedPoint2D) = p.y
GP.getx(p::NumberedPoint2D) = p.x

function rescale( point::NumberedPoint2D, scalefactor::NTuple{4,Float64})
    return NumberedPoint2D((getx(point) - scalefactor[1])/scalefactor[2], (gety(point) - scalefactor[3])/scalefactor[4],point.id)
end


function _delaunay2(x)
    n = length(x)
    a = [NumberedPoint2D(x[i][1],x[i][2], i) for i in 1:n]
    tess = DelaunayTessellation2D{NumberedPoint2D}(n,convex=true)
    push!(tess,a)
    trigs = triangles(tess)
    res = zeros(Int64,length(trigs), 3)
    for (i,t) in enumerate(trigs)
        for j in 1:3
            res[i,j] = t[j].id
        end
    end
    return tess,res
end

struct Delaunay2
    tess::DelaunayTessellation2D{NumberedPoint2D}
    trigs::Array{Int64,2}

    function Delaunay2(x)
        tess, res = _delaunay2(x)
        new(tess,res)
    end

end


function findindex(d::Delaunay2, p)
    ind = _findindex(d.tess,rescale(p,d.tess._scaling))
    if isexternal(d.tess._trigs[ind])
        error("Point is outside of triangulation")
    end
    return d.tess._nonexternal_indexes_inv[ind]
end
