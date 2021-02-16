module NaNStatistics

    # AVX vectorziation tools
    using IfElse: ifelse
    using LoopVectorization
    using VectorizationBase: Vec

    using Statistics
    using StatsBase: percentile, mean, std, ProbabilityWeights

    include("ArrayStats.jl")

end
