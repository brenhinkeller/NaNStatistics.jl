module NaNStatistics

    using Statistics
    using StatsBase: percentile, mean, std, ProbabilityWeights

    # AVX vectorziation tools
    using LoopVectorization
    using SIMDPirates: vifelse, verf
    using VectorizationBase: SVec

    import SpecialFunctions.erf
    erf(x::SVec) = verf(x)

    include("ArrayStats.jl")

end
