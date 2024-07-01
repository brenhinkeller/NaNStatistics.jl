## --- ArrayStats.jl

## --- Filtering

    A = 1:100
    @test A[inpctile(A,80)] == 11:90
    @test nanmask(A) == trues(100)
    @test nanmask([1//1, 1//2, 1//3]) == trues(3)
    @test nanmask([NaN, NaN, NaN]) == falses(3)

## --- NaN handling functions, simple cases

    A = [1:10; fill(NaN,10)]
    B = [fill(NaN,10); 11:20]
    @test countnans(A) === countnotnans(A) == 10
    @test A[nanmask(A)] == 1:10
    @test nanadd(A,B) == 1:20
    @test nanadd!(A,B) == 1:20
    @test zeronan!(B) == [fill(0,10); 11:20]

    # Test arrays that are mostly NaN
    x = [1.:100;]
    x[34:end] .= NaN
    @test nanminimum(x) === 1.0
    @test nanmaximum(x) === 33.0


## --- Summary statistics: simple cases, Float64

    A = [1:10.0..., NaN]
    @test nansum(A) === 55.0
    @test nancumsum(A) == [1.0, 3.0, 6.0, 10.0, 15.0, 21.0, 28.0, 36.0, 45.0, 55.0, 55.0]
    @test nanlogsumexp(A) ≈ 10.45862974442671
    @test nanlogcumsumexp(A) ≈ [1.0, 2.313261687518223, 3.4076059644443806, 4.440189698561196, 5.451914395937593, 6.456193316018123, 7.457762847404243, 8.458339626479004, 9.458551727967379, 10.45862974442671, 10.45862974442671]
    @test isequal(nanlogcumsumexp(A, reverse=true), [10.45862974442671, 10.458551727967379, 10.458339626479004, 10.457762847404243, 10.456193316018123, 10.451914395937594, 10.440189698561195, 10.40760596444438, 10.313261687518223, 10.0, NaN])
    @test nanmean(A) === 5.5
    @test nanmean(A, ones(11)) === 5.5 # weighted
    @test nanrange(A) === 9.0
    @test nanminimum(A) === 1.0
    @test nanmaximum(A) === 10.0
    @test nanargmin(A) === 1
    @test nanargmax(A) === 10
    @test nanextrema(A) === (1.0, 10.0)
    @test nanvar([1,2,3,NaN]) === 1.0
    @test nanstd([1,2,3,NaN]) === 1.0
    @test nansem([1,2,3,NaN]) ≈ 1/sqrt(3)
    @test nanstd([1,2,3,NaN], ones(4)) === 1.0 # weighted
    @test nanmad([1,2,3,NaN]) === 1.0
    @test nanaad([1,2,3,NaN]) ≈ 2/3
    @test nanmedian([1,2,3,NaN]) === 2.0
    @test nanpctile([0:100...,NaN],99) === 99.0
    @test nanmin(1.,2.) === 1.0
    @test nanmax(1.,2.) === 2.0
    @test nanskewness([1,2,3,NaN]) ≈ 0.0
    @test nankurtosis([1,2,3,NaN]) ≈ -1.5

## --- Non-monotonic arrays

    A = [1:5.0..., NaN, 5:-1:1...]
    @test nancumsum(A) ≈ [1.0, 3.0, 6.0, 10.0, 15.0, 15.0, 20.0, 24.0, 27.0, 29.0, 30.0]
    @test nanlogcumsumexp(A) ≈ [1.0, 2.313261687518223, 3.4076059644443806, 4.440189698561196, 5.451914395937593, 5.451914395937593, 5.944418386887494, 6.078136371527456, 6.123152745316754, 6.139216411276987, 6.145061576497539]
    @test nanlogcumsumexp(A, reverse=true) ≈ [6.145061576497539, 6.139216411276987, 6.123152745316754, 6.078136371527456, 5.944418386887494, 5.451914395937593, 5.451914395937593, 4.440189698561196, 3.4076059644443806, 2.313261687518223, 1.0]

    A = [1:10.0..., NaN, 10:-1:1...]
    @test nansum(A) === 110.0
    @test nanlogsumexp(A) ≈ 11.151776924986656
    @test nanrange(A) === 9.0
    @test nanminimum(A) === 1.0
    @test nanmaximum(A) === 10.0
    @test nanextrema(A) === (1.0, 10.0)
    @test nanvar(A) ≈ 8.68421052631579
    @test nanstd(A) ≈ 2.946898458772509
    @test nansem(A) ≈ 2.946898458772509/sqrt(20)
    @test nanmad(A) ≈ 2.5
    @test nanaad(A) ≈ 2.5
    @test nanmedian(A) ≈ 5.5
    @test nanpctile(A,50) ≈ 5.5
    @test nanskewness(A) ≈ 0.0
    @test nankurtosis(A) ≈ -1.2242424242424244

## --- Arrays containing only NaNs should yield NaN (or 0 for sums)

    A = fill(NaN,10)
    @test nansum(A) == 0
    @test nancumsum(A) == zeros(10)
    @test isnan(nanlogsumexp(A))
    @test all(isnan, nanlogcumsumexp(A))
    @test all(isnan, nanlogcumsumexp(A, reverse=true))
    @test isnan(nanmean(A))
    @test isnan(nanmean(A, ones(10))) # weighted
    @test isnan(nanrange(A))
    @test isnan(nanminimum(A))
    @test isnan(nanmaximum(A))
    @test all(isnan, nanextrema(A))
    @test isnan(nanvar(A))
    @test isnan(nanvar(A, mean=NaN))
    @test isnan(nanstd(A))
    @test isnan(nansem(A))
    @test isnan(nanstd(A, ones(10))) # weighted
    @test isnan(nanaad(A))
    @test isnan(nanmad(A))
    @test isnan(nanmedian(A))
    @test isnan(nanpctile(A, 90))
    @test isnan(nanquantile(A, 0.9))
    @test isnan(nanskewness(A))
    @test isnan(nankurtosis(A))
    @test isnan(nanmin(NaN,NaN))
    @test isnan(nanmax(NaN,NaN))

## --- empty arrays should yield NaN (or 0 for sums)

    A = Float64[]
    @test nansum(A) == 0
    @test nancumsum(A) == Float64[]
    @test isnan(nanlogsumexp(A))
    @test nanlogcumsumexp(A) == Float64[]
    @test nanlogcumsumexp(A, reverse=true) == Float64[]
    @test isnan(nanmean(A))
    @test isnan(nanmean(A, copy(A))) # weighted
    @test isnan(nanvar(A))
    @test isnan(nanvar(A, mean=0))
    @test isnan(nanstd(A))
    @test isnan(nansem(A))
    @test isnan(nanstd(A, copy(A))) # weighted
    @test isnan(nanaad(A))
    @test isnan(nanmad(A))
    @test isnan(nanmedian(A))
    @test isnan(nanpctile(A, 90))
    @test isnan(nanquantile(A, 0.9))
    @test isnan(nanskewness(A))
    @test isnan(nankurtosis(A))
    @test nanargmin(A) === 1
    @test nanargmax(A) === 1


## --- Summary statistics: simple cases, Int64

    A = collect(1:10)
    @test nansum(A) === 55
    @test nancumsum(A) == [1, 3, 6, 10, 15, 21, 28, 36, 45, 55]
    @test nanlogsumexp(A) ≈ 10.45862974442671
    @test nanlogcumsumexp(A) ≈ [1.0, 2.313261687518223, 3.4076059644443806, 4.440189698561196, 5.451914395937593, 6.456193316018123, 7.457762847404243, 8.458339626479004, 9.458551727967379, 10.45862974442671]
    @test nanlogcumsumexp(A, reverse=true) ≈ [10.45862974442671, 10.458551727967379, 10.458339626479004, 10.457762847404243, 10.456193316018123, 10.451914395937594, 10.440189698561195, 10.40760596444438, 10.313261687518223, 10.0]
    @test nanmean(A) === 5.5
    @test nanmean(A, ones(10)) === 5.5 # weighted
    @test nanrange(A) === 9
    @test nanminimum(A) === 1
    @test nanmaximum(A) === 10
    @test nanargmin(A) === 1
    @test nanargmax(A) === 10
    @test nanextrema(A) === (1, 10)
    @test nanvar([1,2,3]) === 1.0
    @test nanvar([1,2,3], mean=2) === 1.0
    @test nanstd([1,2,3]) === 1.0
    @test nansem([1,2,3]) ≈ 1/sqrt(3)
    @test nanstd([1,2,3], ones(3)) === 1.0 # weighted
    @test nanmad([1,2,3]) === 1.0
    @test nanaad([1,2,3]) ≈ 2/3
    @test nanmedian([1,2,3]) === 2.0
    @test nanpctile([0:100...],99) === 99.0
    @test nanskewness([1,2,3]) ≈ 0.0
    @test nankurtosis([1,2,3]) ≈ -1.5
    @test nanmin(1,2) === 1
    @test nanmax(1,2) === 2

## --- Summary statistics: simple cases, ranges

    A = 1:10
    @test nansum(A) === 55
    @test nancumsum(A) == [1, 3, 6, 10, 15, 21, 28, 36, 45, 55]
    @test nanlogsumexp(A) ≈ 10.45862974442671
    @test nanlogcumsumexp(A) ≈ [1.0, 2.313261687518223, 3.4076059644443806, 4.440189698561196, 5.451914395937593, 6.456193316018123, 7.457762847404243, 8.458339626479004, 9.458551727967379, 10.45862974442671]
    @test nanlogcumsumexp(A, reverse=true) ≈ [10.45862974442671, 10.458551727967379, 10.458339626479004, 10.457762847404243, 10.456193316018123, 10.451914395937594, 10.440189698561195, 10.40760596444438, 10.313261687518223, 10.0]
    @test nanmean(A) === 5.5
    @test nanmean(A, ones(10)) === 5.5 # weighted
    @test nanrange(A) === 9
    @test nanminimum(A) === 1
    @test nanmaximum(A) === 10
    @test nanargmin(A) === 1
    @test nanargmax(A) === 10
    @test nanextrema(A) === (1, 10)
    @test nanvar(1:3) === 1.0
    @test nanvar(1:3, mean=2) === 1.0
    @test nanstd(1:3) === 1.0
    @test nansem(1:3) ≈ 1/sqrt(3)
    @test nanstd(1:3, ones(3)) === 1.0 # weighted
    @test nanmad(1:3) === 1.0
    @test nanaad(1:3) ≈ 2/3
    @test nanmedian(1:3) === 2.0
    @test nanpctile(0:100,99) === 99.0
    @test nanskewness(1:3) ≈ 0.0
    @test nankurtosis(1:3) ≈ -1.5

    A = 1:10.
    @test nansum(A) === 55.0
    @test nancumsum(A) == [1.0, 3.0, 6.0, 10.0, 15.0, 21.0, 28.0, 36.0, 45.0, 55.0]
    @test nanlogsumexp(A) ≈ 10.45862974442671
    @test nanlogcumsumexp(A) ≈ [1.0, 2.313261687518223, 3.4076059644443806, 4.440189698561196, 5.451914395937593, 6.456193316018123, 7.457762847404243, 8.458339626479004, 9.458551727967379, 10.45862974442671]
    @test nanlogcumsumexp(A, reverse=true) ≈ [10.45862974442671, 10.458551727967379, 10.458339626479004, 10.457762847404243, 10.456193316018123, 10.451914395937594, 10.440189698561195, 10.40760596444438, 10.313261687518223, 10.0]
    @test nanmean(A) === 5.5
    @test nanmean(A, ones(10)) === 5.5 # weighted
    @test nanrange(A) === 9.0
    @test nanminimum(A) === 1.0
    @test nanmaximum(A) === 10.0
    @test nanargmin(A) === 1
    @test nanargmax(A) === 10
    @test nanextrema(A) === (1.0, 10.0)
    @test nanvar(1:3.) === 1.0
    @test nanvar(1:3., mean=2.0) === 1.0
    @test nanstd(1:3.) === 1.0
    @test nansem(1:3.) ≈ 1/sqrt(3)
    @test nanstd(1:3., ones(3)) === 1.0 # weighted
    @test nanmad(1:3.) === 1.0
    @test nanaad(1:3.) ≈ 2/3
    @test nanmedian(1:3.) === 2.0
    @test nanpctile(0:100.,99) === 99.0
    @test nanskewness(1:3.) ≈ 0.0
    @test nankurtosis(1:3.) ≈ -1.5

## --- Summary statistics: dimensional tests, Int64

    A = collect(reshape(1:300,100,3))
    @test nansum(A, dims=1) == sum(A, dims=1)
    @test nansum(A, dims=2) == sum(A, dims=2)
    @test nancumsum(A, dims=1) == cumsum(A, dims=1)
    @test nancumsum(A, dims=2) == cumsum(A, dims=2)
    @test nanminimum(A, dims=1) == minimum(A, dims=1)
    @test nanminimum(A, dims=2) == minimum(A, dims=2)
    @test nanmaximum(A, dims=1) == maximum(A, dims=1)
    @test nanmaximum(A, dims=2) == maximum(A, dims=2)
    @test nanextrema(A, dims=1) == [(1.0, 100.0) (101.0, 200.0) (201.0, 300.0)]
    @test nanmean(A, dims=1) == mean(A, dims=1)
    @test nanmean(A, dims=2) == mean(A, dims=2)
    @test nanmean(A, ones(size(A)), dims=1) == mean(A, dims=1) # weighted
    @test nanmean(A, ones(size(A)), dims=2) == mean(A, dims=2) # weighted
    @test nanvar(A, dims=1) ≈ var(A, dims=1)
    @test nanvar(A, dims=2) ≈ var(A, dims=2)
    @test nanstd(A, dims=1) ≈ std(A, dims=1)
    @test nanstd(A, dims=2) ≈ std(A, dims=2)
    @test nanstd(A, dims=1, mean=nanmean(A,dims=1)) ≈ std(A, dims=1)
    @test nanstd(A, dims=2, mean=nanmean(A,dims=2)) ≈ std(A, dims=2)
    @test nansem(A, dims=1) ≈ std(A, dims=1)./sqrt(size(A,1))
    @test nansem(A, dims=2) ≈ std(A, dims=2)./sqrt(size(A,2))
    @test nansem(A, dims=1, mean=nanmean(A,dims=1)) ≈ std(A, dims=1)./sqrt(size(A,1))
    @test nansem(A, dims=2, mean=nanmean(A,dims=2)) ≈ std(A, dims=2)./sqrt(size(A,2))
    @test nanstd(A, ones(size(A)), dims=1) ≈ std(A, dims=1) # weighted
    @test nanstd(A, ones(size(A)), dims=2) ≈ std(A, dims=2) # weighted
    @test nanmad(A, dims=1) == [25.0 25.0 25.0]
    @test nanmad(A, dims=2) == fill(100.0,100,1)
    @test nanaad(A, dims=1) == [25.0 25.0 25.0] #
    @test nanaad(A, dims=2) ≈ fill(200/3,100,1)
    @test nanmedian(A, dims=1) == median(A, dims=1)
    @test nanmedian(A, dims=2) == median(A, dims=2)
    @test nanpctile(A, 10, dims=1) ≈ [10.9  110.9  210.9]
    @test nanpctile(A, 10, dims=2) ≈ 21:120
    @test vec(nanskewness(A, dims=1)) ≈ skewness.(eachcol(A))
    @test vec(nanskewness(A, dims=2)) ≈ skewness.(eachrow(A))
    @test vec(nankurtosis(A, dims=1)) ≈ kurtosis.(eachcol(A))
    @test vec(nankurtosis(A, dims=2)) ≈ kurtosis.(eachrow(A))

## --- Summary statistics: dimensional tests, Float64

    A = collect(reshape(1:300.,100,3))
    @test nansum(A, dims=1) == sum(A, dims=1)
    @test nansum(A, dims=2) == sum(A, dims=2)
    @test nancumsum(A, dims=1) == cumsum(A, dims=1)
    @test nancumsum(A, dims=2) == cumsum(A, dims=2)
    @test nanminimum(A, dims=1) == minimum(A, dims=1)
    @test nanminimum(A, dims=2) == minimum(A, dims=2)
    @test nanmaximum(A, dims=1) == maximum(A, dims=1)
    @test nanmaximum(A, dims=2) == maximum(A, dims=2)
    @test nanextrema(A, dims=1) == [(1.0, 100.0) (101.0, 200.0) (201.0, 300.0)]
    @test nanmean(A, dims=1) == mean(A, dims=1)
    @test nanmean(A, dims=2) == mean(A, dims=2)
    @test nanmean(A, ones(size(A)), dims=1) == mean(A, dims=1) # weighted
    @test nanmean(A, ones(size(A)), dims=2) == mean(A, dims=2) # weighted
    @test nanvar(A, dims=1) ≈ var(A, dims=1)
    @test nanvar(A, dims=2) ≈ var(A, dims=2)
    @test nanstd(A, dims=1) ≈ std(A, dims=1)
    @test nanstd(A, dims=2) ≈ std(A, dims=2)
    @test nanstd(A, dims=1, mean=nanmean(A,dims=1)) ≈ std(A, dims=1)
    @test nanstd(A, dims=2, mean=nanmean(A,dims=2)) ≈ std(A, dims=2)
    @test nansem(A, dims=1) ≈ std(A, dims=1)./sqrt(size(A,1))
    @test nansem(A, dims=2) ≈ std(A, dims=2)./sqrt(size(A,2))
    @test nansem(A, dims=1, mean=nanmean(A,dims=1)) ≈ std(A, dims=1)./sqrt(size(A,1))
    @test nansem(A, dims=2, mean=nanmean(A,dims=2)) ≈ std(A, dims=2)./sqrt(size(A,2))
    @test nanstd(A, ones(size(A)), dims=1) ≈ std(A, dims=1) # weighted
    @test nanstd(A, ones(size(A)), dims=2) ≈ std(A, dims=2) # weighted
    @test nanmedian(A, dims=1) == median(A, dims=1)
    @test nanmedian(A, dims=2) == median(A, dims=2)
    @test nanminimum(A, dims=1) == [1 101 201]
    @test dropdims(nanminimum(A, dims=2), dims=2) == 1:100
    @test nanmaximum(A, dims=1) == [100 200 300]
    @test dropdims(nanmaximum(A, dims=2), dims=2) == 201:300
    @test nanmean(A, dims=1) == [50.5 150.5 250.5]
    @test dropdims(nanmean(A, dims=2), dims=2) == 101:200
    @test nanmean(A, ones(size(A)), dims=1) == [50.5 150.5 250.5] # weighted
    @test dropdims(nanmean(A, ones(size(A)), dims=2), dims=2) == 101:200 # weighted
    @test nanstd(A, dims=1) ≈ fill(29.011491975882016, 1, 3)
    @test nanstd(A, dims=2) ≈ fill(100, 100, 1)
    @test nanstd(A, ones(size(A)), dims=1) ≈ fill(29.011491975882016, 1, 3) # weighted
    @test nanstd(A, ones(size(A)), dims=2) ≈ fill(100, 100, 1) # weighted
    @test nanmedian(A, dims=1) == [50.5 150.5 250.5]
    @test dropdims(nanmedian(A, dims=2), dims=2) == 101:200
    @test nanpctile(A, 10, dims=1) ≈ [10.9 110.9 210.9]
    @test nanpctile(A, 10, dims=2) ≈ 21:120
    @test nanmad(A, dims=1) == [25.0 25.0 25.0]
    @test nanmad(A, dims=2) == fill(100.0, 100, 1)
    @test nanaad(A, dims=1) == [25.0 25.0 25.0]
    @test nanaad(A, dims=2) ≈ fill(200/3, 100, 1)
    @test vec(nanskewness(A, dims=1)) ≈ skewness.(eachcol(A))
    @test vec(nanskewness(A, dims=2)) ≈ skewness.(eachrow(A))
    @test vec(nankurtosis(A, dims=1)) ≈ kurtosis.(eachcol(A))
    @test vec(nankurtosis(A, dims=2)) ≈ kurtosis.(eachrow(A))

## --- Summary statistics: dimensional tests, Float64, dim

    A = collect(reshape(1:300.,100,3))
    @test nansum(A, dim=1) == vec(sum(A, dims=1))
    @test nansum(A, dim=2) == vec(sum(A, dims=2))
    @test nanminimum(A, dim=1) == vec(minimum(A, dims=1))
    @test nanminimum(A, dim=2) == vec(minimum(A, dims=2))
    @test nanmaximum(A, dim=1) == vec(maximum(A, dims=1))
    @test nanmaximum(A, dim=2) == vec(maximum(A, dims=2))
    @test nanextrema(A, dim=1) == [(1.0, 100.0), (101.0, 200.0), (201.0, 300.0)]
    @test nanmean(A, dim=1) == vec(mean(A, dims=1))
    @test nanmean(A, dim=2) == vec(mean(A, dims=2))
    @test nanmean(A, ones(size(A)), dim=1) == vec(mean(A, dims=1)) # weighted
    @test nanmean(A, ones(size(A)), dim=2) == vec(mean(A, dims=2)) # weighted
    @test nanvar(A, dim=1) ≈ vec(var(A, dims=1))
    @test nanvar(A, dim=2) ≈ vec(var(A, dims=2))
    @test nanvar(A, dim=1, mean=nanmean(A,dims=1)) ≈ vec(var(A, dims=1))
    @test nanvar(A, dim=2, mean=nanmean(A,dims=2)) ≈ vec(var(A, dims=2))
    @test nanstd(A, dim=1) ≈ vec(std(A, dims=1))
    @test nanstd(A, dim=2) ≈ vec(std(A, dims=2))
    @test nansem(A, dim=1) ≈ vec(std(A, dims=1)./sqrt(size(A, 1)))
    @test nansem(A, dim=2) ≈ vec(std(A, dims=2)./sqrt(size(A, 2)))
    @test nanstd(A, ones(size(A)), dim=1) ≈ vec(std(A, dims=1)) # weighted
    @test nanstd(A, ones(size(A)), dim=2) ≈ vec(std(A, dims=2)) # weighted
    @test nanmedian(A, dim=1) == vec(median(A, dims=1))
    @test nanmedian(A, dim=2) == vec(median(A, dims=2))
    @test nanminimum(A, dim=1) == [1, 101, 201]
    @test nanminimum(A, dim=2) == 1:100
    @test nanmaximum(A, dim=1) == [100, 200, 300]
    @test nanmaximum(A, dim=2) == 201:300
    @test nanmean(A, dim=1) == [50.5, 150.5, 250.5]
    @test nanmean(A, dim=2) == 101:200
    @test nanmean(A, ones(size(A)), dim=1) == [50.5, 150.5, 250.5] # weighted
    @test nanmean(A, ones(size(A)), dim=2) == 101:200 # weighted
    @test nanstd(A, dim=1) ≈ fill(29.011491975882016, 3)
    @test nanstd(A, dim=2) ≈ fill(100, 100)
    @test nanstd(A, ones(size(A)), dim=1) ≈ fill(29.011491975882016, 3) # weighted
    @test nanstd(A, ones(size(A)), dim=2) ≈ fill(100, 100) # weighted
    @test nanmedian(A, dim=1) == [50.5, 150.5, 250.5]
    @test nanmedian(A, dim=2) == 101:200
    @test nanpctile(A, 10, dim=1) ≈ [10.9, 110.9, 210.9]
    @test nanpctile(A, 10, dim=2) ≈ 21:120
    @test nanmad(A, dim=1) == [25.0, 25.0, 25.0]
    @test nanmad(A, dim=2) == fill(100.0, 100)
    @test nanaad(A, dim=1) == [25.0, 25.0, 25.0]
    @test nanaad(A, dim=2) ≈ fill(200/3, 100)
    @test nanskewness(A, dim=1) ≈ skewness.(eachcol(A))
    @test nanskewness(A, dim=2) ≈ skewness.(eachrow(A))
    @test nankurtosis(A, dim=1) ≈ kurtosis.(eachcol(A))
    @test nankurtosis(A, dims=2) ≈ kurtosis.(eachrow(A))

## --- Summary statistics: dimensional tests, Float64, nonstandard array type

    A = reshape(1:300.,100,3)
    @test nansum(A, dims=1) == sum(A, dims=1)
    @test nansum(A, dims=2) == sum(A, dims=2)
    @test nancumsum(A, dims=1) == cumsum(A, dims=1)
    @test nancumsum(A, dims=2) == cumsum(A, dims=2)
    @test nanminimum(A, dims=1) == minimum(A, dims=1)
    @test nanminimum(A, dims=2) == minimum(A, dims=2)
    @test nanmaximum(A, dims=1) == maximum(A, dims=1)
    @test nanmaximum(A, dims=2) == maximum(A, dims=2)
    @test nanextrema(A, dims=1) == [(1.0, 100.0) (101.0, 200.0) (201.0, 300.0)]
    @test nanmean(A, dims=1) == mean(A, dims=1)
    @test nanmean(A, dims=2) == mean(A, dims=2)
    @test nanmean(A, ones(size(A)), dims=1) == mean(A, dims=1) # weighted
    @test nanmean(A, ones(size(A)), dims=2) == mean(A, dims=2) # weighted
    @test nanvar(A, dims=1) ≈ var(A, dims=1)
    @test nanvar(A, dims=2) ≈ var(A, dims=2)
    @test nanstd(A, dims=1) ≈ std(A, dims=1)
    @test nanstd(A, dims=2) ≈ std(A, dims=2)
    @test nanstd(A, dims=1, mean=nanmean(A,dims=1)) ≈ std(A, dims=1)
    @test nanstd(A, dims=2, mean=nanmean(A,dims=2)) ≈ std(A, dims=2)
    @test nansem(A, dims=1) ≈ std(A, dims=1)./sqrt(size(A,1))
    @test nansem(A, dims=2) ≈ std(A, dims=2)./sqrt(size(A,2))
    @test nansem(A, dims=1, mean=nanmean(A,dims=1)) ≈ std(A, dims=1)./sqrt(size(A,1))
    @test nansem(A, dims=2, mean=nanmean(A,dims=2)) ≈ std(A, dims=2)./sqrt(size(A,2))
    @test nanstd(A, ones(size(A)), dims=1) ≈ std(A, dims=1) # weighted
    @test nanstd(A, ones(size(A)), dims=2) ≈ std(A, dims=2) # weighted
    @test nanmedian(A, dims=1) == median(A, dims=1)
    @test nanmedian(A, dims=2) == median(A, dims=2)
    @test nanminimum(A, dims=1) == [1 101 201]
    @test dropdims(nanminimum(A, dims=2), dims=2) == 1:100
    @test nanmaximum(A, dims=1) == [100 200 300]
    @test dropdims(nanmaximum(A, dims=2), dims=2) == 201:300
    @test nanmean(A, dims=1) == [50.5 150.5 250.5]
    @test dropdims(nanmean(A, dims=2), dims=2) == 101:200
    @test nanmean(A, ones(size(A)), dims=1) == [50.5 150.5 250.5] # weighted
    @test dropdims(nanmean(A, ones(size(A)), dims=2), dims=2) == 101:200 # weighted
    @test nanstd(A, dims=1) ≈ fill(29.011491975882016, 1, 3)
    @test nanstd(A, dims=2) ≈ fill(100, 100, 1)
    @test nanstd(A, ones(size(A)), dims=1) ≈ fill(29.011491975882016, 1, 3) # weighted
    @test nanstd(A, ones(size(A)), dims=2) ≈ fill(100, 100, 1) # weighted
    @test nanmedian(A, dims=1) == [50.5 150.5 250.5]
    @test dropdims(nanmedian(A, dims=2), dims=2) == 101:200
    @test nanpctile(A, 10, dims=1) ≈ [10.9 110.9 210.9]
    @test nanpctile(A, 10, dims=2) ≈ 21:120
    @test nanmad(A, dims=1) == [25.0 25.0 25.0]
    @test nanmad(A, dims=2) == fill(100.0, 100, 1)
    @test nanaad(A, dims=1) == [25.0 25.0 25.0]
    @test nanaad(A, dims=2) ≈ fill(200/3, 100, 1)
    @test vec(nanskewness(A, dims=1)) ≈ skewness.(eachcol(A))
    @test vec(nanskewness(A, dims=2)) ≈ skewness.(eachrow(A))
    @test vec(nankurtosis(A, dims=1)) ≈ kurtosis.(eachcol(A))
    @test vec(nankurtosis(A, dims=2)) ≈ kurtosis.(eachrow(A))

## --- Test fallbacks for complex reductions

    A = randn((2 .+ (1:6))...);
    @test nansum(A, dims=(4,5,6)) ≈ sum(A, dims=(4,5,6))
    @test nanmean(A, dims=(4,5,6)) ≈ mean(A, dims=(4,5,6))
    @test nanstd(A, dims=(4,5,6)) ≈ std(A, dims=(4,5,6))
    @test nanstd(A, dims=(4,5,6)) ≈ nanstd(A, dims=(4,5,6), mean=nanmean(A, dims=(4,5,6)))
    @test nansem(A, dims=(4,5,6)) ≈ std(A, dims=(4,5,6))./sqrt(size(A,4)*size(A,5)*size(A,6))
    @test nansem(A, dims=(4,5,6)) ≈ nansem(A, dims=(4,5,6), mean=nanmean(A, dims=(4,5,6)))

## --- Test in-place vs. out-of-place versions

    A = rand(100)
    @test nanpctile(A, 25) == nanpctile!(A, 25) == nanquantile!(A, 0.25)
    @test nanmedian(A) == nanmedian!(A)
    @test nanmad(A) == nanmad!(A)


## --- A few tests with other types
    for T in (Bool, Float16, Float32)
        let A = rand(T, 100)
            @test nansum(A) ≈ sum(A)
            @test nanmean(A) ≈ mean(A)
            @test nanstd(A) ≈ std(A)
        end
    end

    for T in (UInt8, UInt16, UInt32, Int8, Int16, Int32)
        let A = collect(T, 1:100)
            @test nansum(A) == 5050
            @test nanmean(A) == 50.5
            @test nanstd(A) ≈ 29.011491975882016
        end
    end

## --- Test stats on nontraditional types

    @test nanminimum((1,2,3,4,5)) === 1
    @test nanmaximum((1,2,3,4,5)) === 5
    @test nanargmin((1,2,3,4,5)) === 1
    @test nanargmax((1,2,3,4,5)) === 5
    @test nanrange((1,2,3,4,5)) === 4
    @test nansum((1,2,3,4,5)) === 15
    @test nanmean((1,2,3,4,5)) === 3.0
    @test nanstd((1,2,3)) === 1.0
    @test nanstd((1,2,3), mean=2.0) === 1.0
    @test nansem((1,2,3)) ≈ 1/sqrt(3)
    @test nansem((1,2,3), mean=2.0) ≈ 1/sqrt(3)

## --- Standardization

    @test nanstandardize!(collect(1:10.)) ≈ ((1:10) .- mean(1:10)) / std(1:10)
    @test nanstandardize(1:10.) ≈ ((1:10) .- mean(1:10)) / std(1:10)

## --- Moving sums and averages

    # Moving sum: 1D
    @test movsum(collect(1:10.),5) == movsum(1:10,5)
    @test movsum(1:10,5) == [6, 10, 15, 20, 25, 30, 35, 40, 34, 27]
    # Moving sum: 2D
    @test movsum(repeat(1:5,1,5),3) == [6 9 9 9 6; 12 18 18 18 12; 18 27 27 27 18; 24 36 36 36 24; 18 27 27 27 18]
    @test movsum(repeat(1:5,1,5)',3) == [6 12 18 24 18; 9 18 27 36 27; 9 18 27 36 27; 9 18 27 36 27; 6 12 18 24 18]

    # Moving average: 1D
    @test movmean(collect(1:10.),5) == movmean(1:10,5)
    @test movmean(1:10,5) == [2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.5, 9.0]
    # Moving average: 2D
    @test movmean(repeat(1:10,1,10),5) == repeat([2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.5, 9.0],1,10)
    @test movmean(repeat(1:10,1,10)',5) == repeat([2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.5, 9.0],1,10)'
    @test movmean([3 5 10 3; 4 2 5 8; 5 6 8 8; 2 6 10 6], 3) ≈ [7//2 29//6 11//2 13//2; 25//6 16//3 55//9 7//1; 25//6 16//3 59//9 15//2; 19//4 37//6 22//3 8//1]

## --- Bugfixes

    a = rand(10,5)
    a[3,1] = NaN
    av = selectdim(a,2,[true,true,true,false,false])
    as = a[:,[true,true,true,false,false]]
    @test !any(isnan, nanmean(av,dims=2))
    @test nanmean(av,dims=2) ≈ nanmean(as, dims=2)

    # ForwardDiff / custom type integration
    using ForwardDiff
    A = fill(ForwardDiff.Dual(-Inf, 0.0), 2, 3)
    @test !any(isnan, nanmean(A; dims = 1))

    # Ensure non-contiguous-dims don't segfault
    x = ones(100, 100, 100)
    @test nanmedian(x, dim=(1, 2)) == ones(100)
    @test nanmedian(x, dim=(2, 3)) == ones(100)
    @test nanmedian(x, dim=(1, 3)) == ones(100)
    @test nanmedian(x, dim=(1, 2, 3)) == fill(1.0)


## --- End of File
