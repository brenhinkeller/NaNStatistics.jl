## --- Binning.jl

    # Prepare test arrays
    M3 = zeros(3)
    N3 = fill(0,3)
    M33 = zeros(3,3)
    N33 = fill(0,3,3)
    r = [17.0 117.0 217.0; 50.0 150.0 250.0; 83.0 183.0 283.0]

    # Means
    @test nanbinmean([1:100..., 1], [1:100..., NaN], 0, 100, 3) == [17, 50, 83]
    @test nanbinmean([1:100..., 1], [1:100..., NaN], range(0,100,length=4)) == [17, 50, 83]
    nanbinmean!(M3, N3, [1:100..., 1],[1:100..., NaN],0,100,3)
    @test M3 == [17, 50, 83]
    @test N3 == fill(33,3)
    @test nanbinmean(1:100, reshape(1:300,100,3), 0, 100, 3) == r
    @test nanbinmean(1:100, reshape(1:300,100,3), range(0,100,length=4)) == r

    # Weighted means
    @test NaNStatistics.nanbinwmean([1:100..., 1], [1:100..., NaN], ones(100), 0,100,3) == [17, 50, 83]
    @test NaNStatistics.nanbinwmean([1:100..., 1], [1:100..., NaN], ones(100), range(0,100,length=4)) == [17, 50, 83]
    @test NaNStatistics.nanbinwmean(1:100, reshape(1:300,100,3), ones(100), 0, 100, 3) == r
    @test NaNStatistics.nanbinwmean(1:100, reshape(1:300,100,3), ones(100), range(0,100,length=4)) == r
    NaNStatistics.nanbinwmean!(M33, N33, 1:100, reshape(1:300,100,3), ones(100), 0, 100, 3)
    @test M33 == r
    @test N33 == fill(33,3,3)

    # Medians
    @test nanbinmedian([1:100..., 1],[1:100..., NaN],0,100,3) == [17, 50, 83]
    @test nanbinmedian([1:100..., 1],[1:100..., NaN],range(0,100,length=4)) == [17, 50, 83]
    nanbinmedian!(M3, N3, [1:100..., 1],[1:100..., NaN],0,100,3)
    @test M3 == [17, 50, 83]
    @test N3 == fill(33,3)
    @test nanbinmedian(1:100, reshape(1:300,100,3), 0, 100, 3) == r
    @test nanbinmedian(1:100, reshape(1:300,100,3), range(0,100,length=4)) == r
    nanbinmedian!(M33, N33, 1:100, reshape(1:300,100,3), 0, 100, 3)
    @test M33 == r
    @test N33 == fill(33,3,3)

    # 2D Means
    x = y = z = 0.5:9.5
    xedges = yedges = 0:10
    mu = nanbinmean(x,y,z,xedges,yedges)
    @test isequal(mu, [((i==j) ? i-0.5 : NaN) for i in 1:10, j in 1:10])

## -- End of File
