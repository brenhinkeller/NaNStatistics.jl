module NaNStatistics

    # AVX vectorziation tools
    using IfElse: ifelse
    using LoopVectorization
    using VectorizationBase: Vec

    using Statistics
    using StatsBase: percentile, mean, std

    include("ArrayStats.jl")
    include("Histograms.jl")
    include("Binning.jl")

end
