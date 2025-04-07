module BaseTest
    using Test, NaNStatistics, StatsBase, Statistics, LinearAlgebra, Distributions

    @testset "ArrayStats" begin include("testArrayStats.jl") end
    @testset "Covariance and Correlation" begin include("testCovCor.jl") end
    @testset "Sorting" begin include("testSorting.jl") end

    @testset "Histograms" begin include("testHistograms.jl") end
    @testset "Binning" begin include("testBinning.jl") end
end

module DimensionalDataTest
    using Test, NaNStatistics, DimensionalData

    @testset "DimensionalData extension" begin include("testDimensionalDataExt.jl") end
end

module UnitfulTest
    using Test, NaNStatistics, Unitful

    @testset "Unitful extension" begin
        x = [1,2,3] * u"nT"

        # Scalar cases
        @test nansum(x) ≈ 6.0 * u"nT"
        @test nancumsum(x) ≈ [1,3,6] * u"nT"
        @test nanmean(x) ≈ 2.0 * u"nT"
        @test nanmedian(x) ≈ 2.0 * u"nT"
        @test nancor(x,x) ≈ 1.0
        @test nanvar(x) ≈ 1.0 * u"nT"^2
        @test nanstd(x) ≈ 1.0 * u"nT"
        @test nansem(x) ≈ 1/sqrt(3) * u"nT"
        @test nanskewness(x) ≈ 0.0
        @test nankurtosis(x) ≈ -1.5
        @test nancov(x,x) ≈ 1.0 * u"nT"^2
        @test nancovem(x,x) ≈ 1/3 * u"nT"^2

        X = [1 2 3; 2 4 6; 3 6 9] * u"nT"
        @test nancor(X) ≈ ones(3,3)
        @test nancov(X) ≈ X * u"nT"
        @test nancovem(X) ≈ X/3 * u"nT"

        # Dimensional tests
        @test nanvar(x, dim=1) ≈ [1.,] .* u"nT"^2
        @test nanvar(X, dim=1) ≈ [1., 4., 9.] .* u"nT"^2
        @test nanstd(x, dim=1) ≈ [1.,] .* u"nT"
        @test nanstd(X, dim=1) ≈ [1., 2., 3.] .* u"nT"
        @test nansem(x, dim=1) ≈ [1.,]/sqrt(3) .* u"nT"
        @test nansem(X, dim=1) ≈ [1., 2., 3.]/sqrt(3) .* u"nT"
        @test nanskewness(x, dim=1) ≈ [0.,]
        @test nanskewness(X, dim=1) ≈ [0., 0., 0.]
        @test nankurtosis(x, dim=1) ≈ [-1.5,]
        @test nankurtosis(X, dim=1) ≈ [-1.5, -1.5, -1.5]
    end
end