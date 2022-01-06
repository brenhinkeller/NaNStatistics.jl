## --- 1D Histograms
    n = 10000
    a = rand(n)
    h = histcounts(a, 0:0.1:1)
    @test h == histcounts(a, 0, 1, 10)
    @test isa(h, Array{Int64,1})
    @test length(h) == 10
    @test sum(h) == n

    @test histcounts(1:100,0:10:100) == fill(10,10)
    @test histcounts(1:100,0,100,10) == fill(10,10)
    @test histcounts(1:100.,0.,100.,10) == fill(10,10)

    N, bin = histcounts(1:100,0:10:100)


    # Test results and warnings when N is too small to hold results
    w = "length(N) < nbins; any bins beyond length(N) will not be filled"
    @test (@test_logs (:warn, w) histcounts!(zeros(5), 1:100 ,0:10:100)) == fill(10,5)
    @test (@test_logs (:warn, w) histcounts!(zeros(5), 1:100. ,0:10:100)) == fill(10,5)

## --- 2D histograms

    x = y = 0.5:9.5
    xedges = yedges = 0:10
    N = histcounts(x,y,xedges,yedges)
    @test N == I(10)

    # Test results and warnings when N is too small to hold results
    w = "size(N,1) < nybins; any y bins beyond size(N,1) will not be filled"
    @test (@test_logs (:warn, w) histcounts!(zeros(5,10), x, y, xedges, yedges)) == I(10)[1:5,:]
    w = "size(N,2) < nxbins; any x bins beyond size(N,2) will not be filled"
    @test (@test_logs (:warn, w) histcounts!(zeros(10,5), x, y, xedges, yedges)) == I(10)[:,1:5]


## --- End of File
