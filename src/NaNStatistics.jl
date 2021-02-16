module NaNStatistics

    using Statistics
    using StatsBase: percentile, mean, std, ProbabilityWeights

    # AVX vectorziation tools
    using LoopVectorization
    using VectorizationBase: Vec



    include("ArrayStats.jl")

end
