
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

