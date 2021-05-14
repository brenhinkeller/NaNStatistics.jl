## --- Binning.jl

    # Prepare test arrays
    M3 = zeros(3)
    N3 = fill(0,3)
    M33 = zeros(3,3)
    N33 = fill(0,3,3)
    r = [17.0 117.0 217.0; 50.0 150.0 250.0; 83.0 183.0 283.0]

    # Means
    @test nanmean([1:100..., 1], [1:100..., NaN], 0, 100, 3) == [17, 50, 83]
    nanmean!(M3, N3, [1:100..., 1],[1:100..., NaN],0,100,3)
    @test M3 == [17, 50, 83]
    @test N3 == fill(33,3)
    @test nanmean(1:100, reshape(1:300,100,3), 0, 100, 3) == r

    # Weighted means
    @test nanmean([1:100..., 1], [1:100..., NaN], ones(100), 0,100,3) == [17, 50, 83]
    @test nanmean(1:100, reshape(1:300,100,3), ones(100), 0, 100, 3) == r
    nanmean!(M33, N33, 1:100, reshape(1:300,100,3), 0, 100, 3)
    @test M33 == r
    @test N33 == fill(33,3,3)

    # Medians
    @test nanmedian([1:100..., 1],[1:100..., NaN],0,100,3) == [17, 50, 83]
    nanmedian!(M3, N3, [1:100..., 1],[1:100..., NaN],0,100,3)
    @test M3 == [17, 50, 83]
    @test N3 == fill(33,3)
    @test nanmedian(1:100, reshape(1:300,100,3), 0, 100, 3) == r
    nanmedian!(M33, N33, 1:100, reshape(1:300,100,3), 0, 100, 3)
    @test M33 == r
    @test N33 == fill(33,3,3)

## -- End of File
