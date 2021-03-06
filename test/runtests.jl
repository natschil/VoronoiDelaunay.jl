using VoronoiDelaunay
import VoronoiDelaunay: _pushunfixed!, _flip!, getneighbor,_findindex,_locate
import GeometricalPredicates
import GeometricalPredicates: incircle, intriangle
using Test


@testset "VoronoiDelaunay tests" begin
    @testset begin
        tess = DelaunayTessellation2D(100)
        @test _findindex(tess, Point2D(1.1, 1.1)) == 2
        @test _findindex(tess, Point2D(1.9, 1.9)) == 3
        @test _findindex(tess, Point2D(1.9, 1.9)) == 3
    end

    pa = Point2D(GeometricalPredicates.min_coord, GeometricalPredicates.min_coord)
    pb = Point2D(GeometricalPredicates.min_coord, GeometricalPredicates.max_coord)
    pc = Point2D(GeometricalPredicates.max_coord, GeometricalPredicates.min_coord)
    pd = Point2D(GeometricalPredicates.max_coord, GeometricalPredicates.max_coord)
    pp = Point2D(1.1,1.1)

    @testset begin
        tess = DelaunayTessellation2D(100)
        @test _findindex(tess, Point2D(1.9, 1.9)) == 3
        @test _findindex(tess, Point2D(1.1, 1.1)) == 2
        @test _findindex(tess, Point2D(1.1, 1.1)) == 2

        @test getneighbor(tess._trigs[2],1) == 3

        _pushunfixed!(tess, Point2D(1.1, 1.1))

        @test getneighbor(tess._trigs[2],1) == 3
        @test getneighbor(tess._trigs[2],2) == 4
        @test getneighbor(tess._trigs[2],3) == 5

        @test getneighbor(tess._trigs[3],1) == 2
        @test getneighbor(tess._trigs[3],2) == 1
        @test getneighbor(tess._trigs[3],3) == 1

        @test getneighbor(tess._trigs[4],1) == 2
        @test getneighbor(tess._trigs[4],2) == 1
        @test getneighbor(tess._trigs[4],3) == 5

        @test getneighbor(tess._trigs[5],1) == 2
        @test getneighbor(tess._trigs[5],2) == 4
        @test getneighbor(tess._trigs[5],3) == 1

        @test getc(tess._trigs[5]) == pp
        @test geta(tess._trigs[2]) == pp
        @test getb(tess._trigs[4]) == pp

        @test getb(tess._trigs[5]) == pb
        @test getb(tess._trigs[2]) == pb
        @test getc(tess._trigs[3]) == pb

        @test getc(tess._trigs[4]) == pc
        @test getc(tess._trigs[2]) == pc
        @test getb(tess._trigs[3]) == pc

        @test geta(tess._trigs[5]) == pa
        @test geta(tess._trigs[4]) == pa

        @test geta(tess._trigs[3]) == pd

        @test _findindex(tess, Point2D(1.01, 1.1)) == 5
        @test _findindex(tess, Point2D(1.1, 1.01)) == 4

        @test _findindex(tess, Point2D(1.11, 1.11)) == 2
        @test _findindex(tess, Point2D(1.6, 1.6)) == 3
        @test _findindex(tess, Point2D(1.11, 1.1101)) == 2
        @test _findindex(tess, Point2D(1.6, 1.601)) == 3
        @test _findindex(tess, Point2D(1.11, 1.11)) == 2
        @test _findindex(tess, Point2D(1.6, 1.6)) == 3

        p2 = Point2D(1.9,1.9)
        _pushunfixed!(tess, p2)

        @test geta(tess._trigs[7]) == pd
        @test getb(tess._trigs[7]) == pc
        @test getc(tess._trigs[7]) == p2

        @test getneighbor(tess._trigs[7],1) == 3
        @test getneighbor(tess._trigs[7],2) == 6
        @test getneighbor(tess._trigs[7],3) == 1

        @test geta(tess._trigs[6]) == pd
        @test getb(tess._trigs[6]) == p2
        @test getc(tess._trigs[6]) == pb

        @test getneighbor(tess._trigs[6],1) == 3
        @test getneighbor(tess._trigs[6],2) == 1
        @test getneighbor(tess._trigs[6],3) == 7

        @test geta(tess._trigs[3]) == p2
        @test getb(tess._trigs[3]) == pc
        @test getc(tess._trigs[3]) == pb

        @test getneighbor(tess._trigs[3],1) == 2
        @test getneighbor(tess._trigs[3],2) == 6
        @test getneighbor(tess._trigs[3],3) == 7

        @test geta(tess._trigs[5]) == pa
        @test getb(tess._trigs[5]) == pb
        @test getc(tess._trigs[5]) == pp

        @test getneighbor(tess._trigs[5],1) == 2
        @test getneighbor(tess._trigs[5],2) == 4
        @test getneighbor(tess._trigs[5],3) == 1

        @test geta(tess._trigs[4]) == pa
        @test getb(tess._trigs[4]) == pp
        @test getc(tess._trigs[4]) == pc

        @test getneighbor(tess._trigs[4],1) == 2
        @test getneighbor(tess._trigs[4],2) == 1
        @test getneighbor(tess._trigs[4],3) == 5

        @test geta(tess._trigs[2]) == pp
        @test getb(tess._trigs[2]) == pb
        @test getc(tess._trigs[2]) == pc

        @test getneighbor(tess._trigs[2],1) == 3
        @test getneighbor(tess._trigs[2],2) == 4
        @test getneighbor(tess._trigs[2],3) == 5

        @test tess._last_trig_index == 7


        _flip!(tess,1, 2, 3)

        @test geta(tess._trigs[2]) == pp
        @test getb(tess._trigs[2]) == p2
        @test getc(tess._trigs[2]) == pc

        @test getneighbor(tess._trigs[2],1) == 7
        @test getneighbor(tess._trigs[2],2) == 4
        @test getneighbor(tess._trigs[2],3) == 3

        @test geta(tess._trigs[3]) == pp
        @test getb(tess._trigs[3]) == pb
        @test getc(tess._trigs[3]) == p2

        @test getneighbor(tess._trigs[3],1) == 6
        @test getneighbor(tess._trigs[3],2) == 2
        @test getneighbor(tess._trigs[3],3) == 5

        @test geta(tess._trigs[5]) == pa
        @test getb(tess._trigs[5]) == pb
        @test getc(tess._trigs[5]) == pp

        @test geta(tess._trigs[6]) == pd
        @test getb(tess._trigs[6]) == p2
        @test getc(tess._trigs[6]) == pb

        @test geta(tess._trigs[4]) == pa
        @test getb(tess._trigs[4]) == pp
        @test getc(tess._trigs[4]) == pc

        @test geta(tess._trigs[7]) == pd
        @test getb(tess._trigs[7]) == pc
        @test getc(tess._trigs[7]) == p2

        @test getneighbor(tess._trigs[4],1) == 2
        @test getneighbor(tess._trigs[4],2) == 1
        @test getneighbor(tess._trigs[4],3) == 5

        @test getneighbor(tess._trigs[5],1) == 3
        @test getneighbor(tess._trigs[5],2) == 4
        @test getneighbor(tess._trigs[5],3) == 1

        @test getneighbor(tess._trigs[6],1) == 3
        @test getneighbor(tess._trigs[6],2) == 1
        @test getneighbor(tess._trigs[6],3) == 7

        @test getneighbor(tess._trigs[7],1) == 2
        @test getneighbor(tess._trigs[7],2) == 6
        @test getneighbor(tess._trigs[7],3) == 1

        _flip!(tess,2, 3, 2)

        @test geta(tess._trigs[3]) == pp
        @test getb(tess._trigs[3]) == pb
        @test getc(tess._trigs[3]) == pc

        @test geta(tess._trigs[2]) == pc
        @test getb(tess._trigs[2]) == pb
        @test getc(tess._trigs[2]) == p2

        @test getneighbor(tess._trigs[3],1) == 2
        @test getneighbor(tess._trigs[3],2) == 4
        @test getneighbor(tess._trigs[3],3) == 5

        @test getneighbor(tess._trigs[2],1) == 6
        @test getneighbor(tess._trigs[2],2) == 7
        @test getneighbor(tess._trigs[2],3) == 3

        @test getneighbor(tess._trigs[6],1) == 2
        @test getneighbor(tess._trigs[7],1) == 2
        @test getneighbor(tess._trigs[5],1) == 3
        @test getneighbor(tess._trigs[4],1) == 3

        _flip!(tess,3,2,3)

        @test geta(tess._trigs[2]) == pp
        @test getb(tess._trigs[2]) == pb
        @test getc(tess._trigs[2]) == p2

        @test geta(tess._trigs[3]) == pc
        @test getb(tess._trigs[3]) == pp
        @test getc(tess._trigs[3]) == p2

        @test getneighbor(tess._trigs[2],1) == 6
        @test getneighbor(tess._trigs[2],2) == 3
        @test getneighbor(tess._trigs[2],3) == 5

        @test getneighbor(tess._trigs[3],1) == 2
        @test getneighbor(tess._trigs[3],2) == 7
        @test getneighbor(tess._trigs[3],3) == 4

        @test getneighbor(tess._trigs[6],1) == 2
        @test getneighbor(tess._trigs[7],1) == 3
        @test getneighbor(tess._trigs[5],1) == 2
        @test getneighbor(tess._trigs[4],1) == 3
    end

    @testset begin
        tess = DelaunayTessellation2D(100)
        pp = Point2D(1.45, 1.49)
        _pushunfixed!(tess, pp)
        _flip!(tess,1,2,3)

        @test geta(tess._trigs[2]) == pp
        @test getb(tess._trigs[2]) == pd
        @test getc(tess._trigs[2]) == pc

        @test geta(tess._trigs[3]) == pp
        @test getb(tess._trigs[3]) == pb
        @test getc(tess._trigs[3]) == pd
    end

    @testset begin
        point_arr = Point2D[]
        n=1000
        tess = DelaunayTessellation2D(n*10)
        for i in 1:n
            push!(point_arr, Point2D(rand()+1.0, rand()+1.0))
        end
        push!(tess, point_arr)
        @test tess._last_trig_index == n*2+3

        for p in point_arr
            for t in tess._trigs[2:tess._last_trig_index]
                i = incircle(t, p)
                if i > 0 && ((p == geta(t)) || (p == getb(t)) || (p == getc(t)))
                    i = 0.
                end
                @test i <= 0
            end
        end
    end

    @testset begin
        point_arr = Point2D[]
        n=10
        tess = DelaunayTessellation2D(n*n*10)
        for x in range(1.001,stop=1.999,length=n)
            for y in range(1.001,stop=1.999,length=n)
                push!(point_arr, Point2D(x,y))
            end
        end
        push!(tess, point_arr)
        @test tess._last_trig_index == n*n*2+3

        for p in point_arr
            for t in tess._trigs[2:tess._last_trig_index]
                i = incircle(t, p)
                if i > 0 && ((p == geta(t)) || (p == getb(t)) || (p == getc(t)))
                    i = 0.
                end
                @test i <= 0
            end
        end
    end
    
    # Iterator test
    @testset begin
        point_arr = Point2D[]
        n=1000
        tess = DelaunayTessellation2D(n*10)
        for i in 1:n
            push!(point_arr, Point2D(rand()+1.0, rand()+1.0))
        end
        push!(tess, point_arr)
        p = Point2D(rand()+1.0, rand()+1.0)
        counter = 0
        for t in tess
            if intriangle(t, p) == 1
                counter += 1
            end
        end
        @test counter == 1 # p can be contained only in one triangle
    end

    @testset begin
        points = [  Point2D(-1.0563841812533212, -1.4606363138997696) 
                    Point2D(0.06312982975128989, -0.48031801152366027)
                    Point2D(0.1624918689993189, -0.19919450833195906) 
                    Point2D(-1.5293344878962758, -0.7657808444340142) 
                    Point2D(0.5319064220493406, 0.6107808132476504)   
                    Point2D(-0.3670342825169435, 0.8207427582546951)  
                    Point2D(-1.9797019290444364, -0.5066353099040788) 
                    Point2D(-1.5811584920674324, 1.0346409888830976)  
                    Point2D(1.2185165319349451, 1.4177209167374605)   
                    Point2D(-1.5991536318191626, -1.3063986775765466) ];
        tess = DelaunayTessellation( points )
        t = locate( tess, Point2D(0.6, 0.6) )
        @test Point2D(0.16249186899931745, -0.19919450833195929) ∈ t
    end

end
# that's it for today!
