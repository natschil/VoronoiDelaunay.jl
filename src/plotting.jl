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

