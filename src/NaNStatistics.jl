module NaNStatistics

    using Statistics
    using StatsBase: percentile, mean, std, ProbabilityWeights

    # AVX vectorziation tools
    using LoopVectorization
    using VectorizationBase

    include("ArrayStats.jl")

end
