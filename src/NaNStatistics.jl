module NaNStatistics

    struct _StaticInt{N} end
    _StaticInt(N::Int) = _StaticInt{N}()
    const _IntOrStaticInt = Union{Integer, _StaticInt}
    _dim(::Type{_StaticInt{N}}) where {N} = N::Int

    struct _True end
    struct _False end
    _static(b::Bool) = b ? _True() : _False()

    include("ArrayStats/ArrayStats.jl")
    include("ArrayStats/nanmean.jl")
    include("ArrayStats/nansum.jl")
    include("ArrayStats/nancumsum.jl")
    include("ArrayStats/nanlogsumexp.jl")
    include("ArrayStats/nanvar.jl")
    include("ArrayStats/nanstd.jl")
    include("ArrayStats/nanskewness.jl")
    include("ArrayStats/nankurtosis.jl")
    include("ArrayStats/nansem.jl")
    include("ArrayStats/nancov.jl")
    include("Sorting/quicksort.jl")
    include("Sorting/nanmedian.jl")
    include("Sorting/nanpctile.jl")
    include("Histograms.jl")
    include("Binning.jl")

    using PrecompileTools
    include("precompile.jl")

end
