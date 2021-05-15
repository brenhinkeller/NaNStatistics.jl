## --- Binning.jl

    # Prepare test arrays
    M3 = zeros(3)
    N3 = fill(0,3)
    M33 = zeros(3,3)
    N33 = fill(0,3,3)
    r = [17.0 117.0 217.0; 50.0 150.0 250.0; 83.0 183.0 283.0]

    # Means
    @test nanbinmean([1:100..., 1], [1:100..., NaN], 0, 100, 3) == [17, 50, 83]
    nanbinmean!(M3, N3, [1:100..., 1],[1:100..., NaN],0,100,3)
    @test M3 == [17, 50, 83]
    @test N3 == fill(33,3)
    @test nanbinmean(1:100, reshape(1:300,100,3), 0, 100, 3) == r

    # Weighted means
    @test NaNStatistics.nanbinwmean([1:100..., 1], [1:100..., NaN], ones(100), 0,100,3) == [17, 50, 83]
    @test NaNStatistics.nanbinwmean(1:100, reshape(1:300,100,3), ones(100), 0, 100, 3) == r
    NaNStatistics.nanbinwmean!(M33, N33, 1:100, reshape(1:300,100,3), ones(100), 0, 100, 3)
    @test M33 == r
    @test N33 == fill(33,3,3)

    # Medians
    @test nanbinmedian([1:100..., 1],[1:100..., NaN],0,100,3) == [17, 50, 83]
    nanbinmedian!(M3, N3, [1:100..., 1],[1:100..., NaN],0,100,3)
    @test M3 == [17, 50, 83]
    @test N3 == fill(33,3)
    @test nanbinmedian(1:100, reshape(1:300,100,3), 0, 100, 3) == r
    nanbinmedian!(M33, N33, 1:100, reshape(1:300,100,3), 0, 100, 3)
    @test M33 == r
    @test N33 == fill(33,3,3)

## -- End of File
