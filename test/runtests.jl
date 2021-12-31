using NaNStatistics
using Test, Statistics, LinearAlgebra

@testset "ArrayStats" begin include("testArrayStats.jl") end
@testset "Covariance and Correlation" begin include("testCovCor.jl") end
@testset "Sorting" begin include("testSorting.jl") end

@testset "Histograms" begin include("testHistograms.jl") end
@testset "Binning" begin include("testBinning.jl") end
