module NaNStatistics

    # AVX vectorziation tools
    using IfElse: ifelse
    using LoopVectorization
    using VectorizationBase: Vec

    using Static
    _dim(::Type{StaticInt{N}}) where {N} = N::Int

    using Statistics
    using StatsBase: percentile, mean, std

    include("ArrayStats/ArrayStats.jl")
    include("ArrayStats/nanmean.jl")
    include("ArrayStats/nansum.jl")
    include("ArrayStats/nanvar.jl")
    include("ArrayStats/nanstd.jl")
    include("ArrayStats/nancov.jl")
    include("Histograms.jl")
    include("Binning.jl")

end
