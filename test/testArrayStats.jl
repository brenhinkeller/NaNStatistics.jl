## --- ArrayStats.jl

    # Filtering
    A = 1:100
    @test A[inpctile(A,80)] == 11:90

    # NaN handling functions, simple cases
    A = [1:10; fill(NaN,10)]
    B = [fill(NaN,10); 11:20]
    @test A[nanmask(A)] == 1:10
    @test nanadd(A,B) == 1:20
    @test nanadd!(A,B) == 1:20
    @test zeronan!(B) == [fill(0,10); 11:20]

    # Summary statistics: simple cases, Float64
    A = [1:10.0..., NaN]
    @test nansum(A) == 55.0
    @test nanmean(A) == 5.5
    @test nanmean(A, ones(11)) == 5.5 # weighted
    @test nanrange(A) == 9.0
    @test nanminimum(A) == 1.0
    @test nanmaximum(A) == 10.0
    @test nanextrema(A) == (1.0, 10.0)
    @test nanstd([1,2,3,NaN]) == 1.0
    @test nanstd([1,2,3,NaN], ones(4)) == 1.0 # weighted
    @test nanmad([1,2,3,NaN]) == 1.0
    @test nanaad([1,2,3,NaN]) ≈ 2/3
    @test nanmedian([1,2,3,NaN]) == 2.0
    @test pctile([0:100...,NaN],99) == 99.0
    @test StatGeochem.nanmin(1.,2.) == 1.
    @test StatGeochem.nanmax(1.,2.) == 2.

    # Arrays containing only NaNs should yield NaN
    A = fill(NaN,10)
    @test nansum(A) == 0
    @test isnan(nanmean(A))
    @test isnan(nanmean(A, ones(10))) # weighted
    @test isnan(nanrange(A))
    @test isnan(nanminimum(A))
    @test isnan(nanmaximum(A))
    @test all(isnan.(nanextrema(A)))
    @test isnan(nanstd(A))
    @test isnan(nanstd(A, ones(10))) # weighted
    @test isnan(nanaad(A))
    @test isnan(StatGeochem.nanmin(NaN,NaN))
    @test isnan(StatGeochem.nanmax(NaN,NaN))

    # Summary statistics: simple cases, Int64
    A = collect(1:10)
    @test nansum(A) == 55.0
    @test nanmean(A) == 5.5
    @test nanmean(A, ones(10)) == 5.5 # weighted
    @test nanrange(A) == 9.0
    @test nanminimum(A) == 1.0
    @test nanmaximum(A) == 10.0
    @test nanextrema(A) == (1.0, 10.0)
    @test nanstd([1,2,3]) == 1.0
    @test nanstd([1,2,3], ones(3)) == 1.0 # weighted
    @test nanmad([1,2,3]) == 1.0
    @test nanaad([1,2,3]) ≈ 2/3
    @test nanmedian([1,2,3]) == 2.0
    @test pctile([0:100...],99) == 99.0
    @test StatGeochem.nanmin(1,2) == 1
    @test StatGeochem.nanmax(1,2) == 2

    # Summary statistics: simple cases, ranges
    A = 1:10
    @test nansum(A) == 55.0
    @test nanmean(A) == 5.5
    @test nanmean(A, ones(10)) == 5.5 # weighted
    @test nanrange(A) == 9.0
    @test nanminimum(A) == 1.0
    @test nanmaximum(A) == 10.0
    @test nanextrema(A) == (1.0, 10.0)
    @test nanstd(1:3) == 1.0
    @test nanstd(1:3, ones(3)) == 1.0 # weighted
    @test nanmad(1:3) == 1.0
    @test nanaad(1:3) ≈ 2/3
    @test nanmedian(1:3) == 2.0
    @test pctile(0:100,99) == 99.0


    A = 1:10.
    @test nansum(A) == 55.0
    @test nanmean(A) == 5.5
    @test nanmean(A, ones(10)) == 5.5 # weighted
    @test nanrange(A) == 9.0
    @test nanminimum(A) == 1.0
    @test nanmaximum(A) == 10.0
    @test nanextrema(A) == (1.0, 10.0)
    @test nanstd(1:3.) == 1.0
    @test nanstd(1:3., ones(3)) == 1.0 # weighted
    @test nanmad(1:3.) == 1.0
    @test nanaad(1:3.) ≈ 2/3
    @test nanmedian(1:3.) == 2.0
    @test pctile(0:100.,99) == 99.0

    # Summary statistics: dimensional tests, Int64
    A = reshape(1:300,100,3)
    @test nanminimum(A, dims=1) == minimum(A, dims=1)
    @test nanminimum(A, dims=2) == minimum(A, dims=2)
    @test nanmaximum(A, dims=1) == maximum(A, dims=1)
    @test nanmaximum(A, dims=2) == maximum(A, dims=2)
    @test nanmean(A, dims=1) == mean(A, dims=1)
    @test nanmean(A, dims=2) == mean(A, dims=2)
    @test nanmean(A, ones(size(A)), dims=1) == mean(A, dims=1) # weighted
    @test nanmean(A, ones(size(A)), dims=2) == mean(A, dims=2) # weighted
    @test nanstd(A, dims=1) ≈ std(A, dims=1)
    @test nanstd(A, dims=2) ≈ std(A, dims=2)
    @test nanstd(A, ones(size(A)), dims=1) ≈ std(A, dims=1) # weighted
    @test nanstd(A, ones(size(A)), dims=2) ≈ std(A, dims=2) # weighted
    @test nanmad(A, dims=1) == [25.0 25.0 25.0]
    @test nanmad(A, dims=2) == fill(100.0,100,1)
    @test nanaad(A, dims=1) == [25.0 25.0 25.0]
    @test nanaad(A, dims=2) ≈ fill(200/3,100,1)
    @test nanmedian(A, dims=1) == median(A, dims=1)
    @test nanmedian(A, dims=2) == median(A, dims=2)
    @test pctile(A, 10, dims=1) == [10.9  110.9  210.9]
    @test pctile(A, 10, dims=2) ≈ 21:120


    # Summary statistics: dimensional tests, Float64
    A = reshape(1:300.,100,3)
    @test nanminimum(A, dims=1) == minimum(A, dims=1)
    @test nanminimum(A, dims=2) == minimum(A, dims=2)
    @test nanmaximum(A, dims=1) == maximum(A, dims=1)
    @test nanmaximum(A, dims=2) == maximum(A, dims=2)
    @test nanmean(A, dims=1) == mean(A, dims=1)
    @test nanmean(A, dims=2) == mean(A, dims=2)
    @test nanmean(A, ones(size(A)), dims=1) == mean(A, dims=1) # weighted
    @test nanmean(A, ones(size(A)), dims=2) == mean(A, dims=2) # weighted
    @test nanstd(A, dims=1) ≈ std(A, dims=1)
    @test nanstd(A, dims=2) ≈ std(A, dims=2)
    @test nanstd(A, ones(size(A)), dims=1) ≈ std(A, dims=1) # weighted
    @test nanstd(A, ones(size(A)), dims=2) ≈ std(A, dims=2) # weighted
    @test nanmad(A, dims=1) == [25.0 25.0 25.0]
    @test nanmad(A, dims=2) == fill(100.0,100,1)
    @test nanaad(A, dims=1) == [25.0 25.0 25.0]
    @test nanaad(A, dims=2) ≈ fill(200/3,100,1)
    @test nanmedian(A, dims=1) == median(A, dims=1)
    @test nanmedian(A, dims=2) == median(A, dims=2)
    @test pctile(A, 10, dims=1) == [10.9 110.9 210.9]
    @test pctile(A, 10, dims=2) ≈ 21:120

    # Same but with other keyword, for 1-d output
    @test nanminimum(A, dim=1) == [1,101,201]
    @test nanminimum(A, dim=2) == 1:100
    @test nanmaximum(A, dim=1) == [100,200,300]
    @test nanmaximum(A, dim=2) == 201:300
    @test nanmean(A, dim=1) == [50.5, 150.5, 250.5]
    @test nanmean(A, dim=2) == 101:200
    @test nanmean(A, ones(size(A)), dim=1) == [50.5, 150.5, 250.5] # weighted
    @test nanmean(A, ones(size(A)), dim=2) == 101:200 # weighted
    @test nanstd(A, dim=1) ≈ fill(29.011491975882016,3)
    @test nanstd(A, dim=2) ≈ fill(100,100)
    @test nanstd(A, ones(size(A)), dim=1) ≈ fill(29.011491975882016,3) # weighted
    @test nanstd(A, ones(size(A)), dim=2) ≈ fill(100,100) # weighted
    @test nanmedian(A, dim=1) == [50.5, 150.5, 250.5]
    @test nanmedian(A, dim=2) == 101:200
    @test pctile(A, 10, dim=1) == [10.9, 110.9, 210.9]
    @test pctile(A, 10, dim=2) ≈ 21:120

    # Summary statistics: binning
    # Means
    @test nanmean([1:100..., 1], [1:100..., NaN], 0, 100, 3) == [17, 50, 83]
    @test nanmean(1:100, reshape(1:300,100,3), 0, 100, 3) ==
                [17.0 117.0 217.0; 50.0 150.0 250.0; 83.0 183.0 283.0]
    # Weighted means
    @test nanmean([1:100..., 1], [1:100..., NaN], ones(100), 0,100,3) == [17, 50, 83]
    @test nanmean(1:100, reshape(1:300,100,3), ones(100), 0, 100, 3) ==
                [17.0 117.0 217.0; 50.0 150.0 250.0; 83.0 183.0 283.0]
    # Medians
    @test nanmedian([1:100..., 1],[1:100..., NaN],0,100,3) == [17, 50, 83.5]
    @test nanmedian(1:100, reshape(1:300,100,3), 0, 100, 3) ==
                [17.0 117.0 217.0; 50.0 150.0 250.0; 83.5 183.5 283.5]

    # Standardization
    @test standardize!(collect(1:10.)) ≈ ((1:10) .- mean(1:10)) / std(1:10)

    # Sorting, counting, matching
    A = rand(1:100.,100); B = sort(A)
    @test A[1:count_unique!(A)] == unique(B)
    @test findclosest(3.3:5.3,1:10) == 3:5
    @test findclosestbelow(3.3:5.3,1:10) == 3:5
    @test findclosestabove(3.3:5.3,1:10) == 4:6
    @test findmatches(40:60,1:100) == 40:60
    @test findnth(fill(true,50), 25) == 25

    # Interpolation
    @test linterp1(1:10,21:30,5:0.5:6) == [25.0, 25.5, 26.0]
    @test linterp1s(10:-1:1,21:30,5:0.5:6) == [26.0, 25.5, 25.0]
    @test linterp_at_index(1:100,10) == 10
    @test movmean(collect(1:10.),5) == movmean(1:10,5)
    @test movmean(1:10,4) == [2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.5, 9.0]
    @test movmean(repeat(1:10,1,10),4) == repeat([2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.5, 9.0],1,10)
    @test cntr(0:2:100) == 1:2:99

    # Integration
    @test trapz(1:10,fill(1,10)) == 9
    @test trapz(collect(1:10.),ones(10)) == 9
    @test midpointintegrate(1:10,ones(10)) == 10

    # Distributions
    A = draw_from_distribution(ones(100), 10000)::AbstractArray
    @test length(A) == 10000
    @test isapprox(mean(A), 0.5, atol=0.08)
    @test isapprox(std(A), sqrt(1/12), atol=0.08)

## ---
