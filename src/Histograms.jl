
"""
```julia
histcounts(x, xedges::AbstractRange; T=Int64)::Vector{T}
```
A 1D histogram, ignoring `NaN`s: calculate the number of `x` values that fall into
each of `length(xedges)-1` equally spaced bins along the `x` axis with bin edges
specified by `xedges`.

By default, the counts are returned as `Int64`s, though this can be changed by
specifying an output type with the optional keyword argument `T`.

## Examples
```julia
julia> b = 10 * rand(100000);

julia> histcounts(b, 0:1:10)
10-element Vector{Int64}:
 10054
  9987
  9851
  9971
  9832
 10033
 10250
 10039
  9950
 10033
```
"""
function histcounts(x, xedges::AbstractRange; T=Int64)
    N = fill(zero(T), length(xedges)-1)
    histcounts!(N, x, xedges)
    return N
end
histcounts(x, xmin::Number, xmax::Number, nbins::Integer; T=Int64) = histcounts(x, range(xmin, xmax, length=nbins+1); T=T)


"""
```julia
histmean(counts, bincenters)
```
Estimate the mean of the data represented by a histogram, 
specified as `counts` in equally spaced bins centered at `bincenters`.

## Examples
```julia
julia> binedges = -10:0.01:10;

julia> counts = histcounts(randn(10000), binedges);

julia> bincenters = (binedges[1:end-1] + binedges[2:end])/2
-9.995:0.01:9.995

julia> histmean(counts, bincenters)
0.0039890000000003135
```
"""
function histmean(counts, bincenters)
    N = ∅ₙ = zero(eltype(counts))
    Σ = ∅ = zero(Base.promote_op(*, eltype(counts), eltype(bincenters)))
    @inbounds @simd ivdep for i in eachindex(counts, bincenters)
        pᵢ = counts[i] * bincenters[i]
        Σ += ifelse(isnan(pᵢ), ∅, pᵢ)
        N += ifelse(isnan(pᵢ), ∅ₙ, counts[i])
    end
    return Σ / N
end
export histmean


"""
```julia
histvar(counts, bincenters; corrected::Bool=true)
```
Estimate the variance of the data represented by a histogram, 
specified as `counts` in equally spaced bins centered at `bincenters`.

If `counts` have been normalized, or represent an analytical estimate of a PDF 
rather than a histogram representing counts of a dataset, Bessel's correction 
to the variance should likely not be performed - i.e., set the 
`corrected` keyword argument to `false`. 

## Examples
```julia
julia> binedges = -10:0.01:10;

julia> counts = histcounts(randn(10000), binedges);

julia> bincenters = (binedges[1:end-1] + binedges[2:end])/2
-9.995:0.01:9.995
t
julia> histvar(counts, bincenters)
0.9991854064196424
```
"""
function histvar(counts, bincenters; corrected::Bool=true)
    N = ∅ₙ = zero(eltype(counts))
    Σ = ∅ = zero(Base.promote_op(*, eltype(counts), eltype(bincenters)))
    @inbounds @simd ivdep for i in eachindex(counts, bincenters)
        pᵢ = counts[i] * bincenters[i]
        Σ += ifelse(isnan(pᵢ), ∅, pᵢ)
        N += ifelse(isnan(pᵢ), ∅ₙ, counts[i])
    end
    μ = Σ/N
    Σ = ∅
    @inbounds @simd ivdep for i in eachindex(counts, bincenters)
        dᵢ = bincenters[i] - μ
        pᵢ = counts[i] * dᵢ * dᵢ
        Σ += ifelse(isnan(pᵢ), ∅, pᵢ)
    end
    return Σ / max(N-corrected, ∅ₙ)
end
export histvar


"""
```julia
histskewness(counts, bincenters)
```
Estimate the skewness of the data represented by a histogram, 
specified as `counts` in equally spaced bins centered at `bincenters`.

## Examples
```julia
julia> binedges = -10:0.01:10;

julia> counts = histcounts(randn(10000), binedges);

julia> bincenters = (binedges[1:end-1] + binedges[2:end])/2
-9.995:0.01:9.995

julia> histskewness(counts, bincenters)
0.011075369240851738
```
"""
function histskewness(counts, bincenters; corrected::Bool=false)
    N = ∅ₙ = zero(eltype(counts))
    Σ = ∅ = zero(Base.promote_op(*, eltype(counts), eltype(bincenters)))
    @inbounds @simd ivdep for i in eachindex(counts, bincenters)
        pᵢ = counts[i] * bincenters[i]
        Σ += ifelse(isnan(pᵢ), ∅, pᵢ)
        N += ifelse(isnan(pᵢ), ∅ₙ, counts[i])
    end
    μ = Σ/N
    μ₂ = μ₃ = ∅
    @inbounds @simd ivdep for i in eachindex(counts, bincenters)
        dᵢ = bincenters[i] - μ
        δ² = dᵢ * dᵢ
        pᵢ = counts[i] * δ²
        μ₂ += ifelse(isnan(pᵢ), ∅, pᵢ)
        μ₃ += ifelse(isnan(pᵢ), ∅, pᵢ*dᵢ)
    end
    σ = sqrt(μ₂ / max(N-corrected, ∅ₙ))
    return (μ₃ / N)/σ^3
end
export histskewness

"""
```julia
histkurtosis(counts, bincenters)
```
Estimate the excess kurtosis [1] of the data represented by a histogram, 
specified as `counts` in equally spaced bins centered at `bincenters`.

[1] We follow Distributions.jl in returning excess kurtosis rather than raw kurtosis.
Excess kurtosis is defined as as kurtosis - 3, such that a Normal distribution
has zero excess kurtosis. 

## Examples
```julia
julia> binedges = -10:0.01:10;

julia> counts = histcounts(randn(10000), binedges);

julia> bincenters = (binedges[1:end-1] + binedges[2:end])/2
-9.995:0.01:9.995
t
julia> histkurtosis(counts, bincenters)
0.028863400305099596
```
"""
function histkurtosis(counts, bincenters; corrected::Bool=false)
    N = ∅ₙ = zero(eltype(counts))
    Σ = ∅ = zero(Base.promote_op(*, eltype(counts), eltype(bincenters)))
    @inbounds @simd ivdep for i in eachindex(counts, bincenters)
        pᵢ = counts[i] * bincenters[i]
        Σ += ifelse(isnan(pᵢ), ∅, pᵢ)
        N += ifelse(isnan(pᵢ), ∅ₙ, counts[i])
    end
    μ = Σ/N
    μ₂ = μ₄ = ∅
    @inbounds @simd ivdep for i in eachindex(counts, bincenters)
        dᵢ = bincenters[i] - μ
        δ² = dᵢ * dᵢ
        pᵢ = counts[i] * δ²
        μ₂ += ifelse(isnan(pᵢ), ∅, pᵢ)
        μ₄ += ifelse(isnan(pᵢ), ∅, pᵢ*δ²)
    end
    σ = sqrt(μ₂ / max(N-corrected, ∅ₙ))
    return (μ₄ / N)/σ^4 - 3
end
export histkurtosis


"""
```julia
histstd(counts, bincenters; corrected::Bool=true)
```
Estimate the standard deviation of the data represented by a histogram, 
specified as `counts` in equally spaced bins centered at `bincenters`.

If `counts` have been normalized, or represent an analytical estimate of a PDF 
rather than a histogram representing counts of a dataset, Bessel's correction 
to the standard deviation should likely not be performed - i.e., set the 
`corrected` keyword argument to `false`. 

## Examples
```julia
julia> binedges = -10:0.01:10;

julia> counts = histcounts(randn(10000), binedges);

julia> bincenters = (binedges[1:end-1] + binedges[2:end])/2
-9.995:0.01:9.995
t
julia> histstd(counts, bincenters)
0.999592620230683
```
"""
histstd(counts, bincenters; corrected::Bool=true) = sqrt(histvar(counts, bincenters; corrected))
export histstd



"""
```julia
histcountindices(x, xedges::AbstractRange; T=Int64)::Vector{T}
```
A 1D histogram, ignoring `NaN`s; as `histcounts` but also returning a vector of
the bin index of each `x` value.

## Examples
```julia
julia> b = 10 * rand(100000);

julia> N, bin = histcountindices(b, 0:2:10)
([20082, 19971, 20049, 19908, 19990], [2, 3, 2, 2, 4, 2, 3, 3, 2, 4  …  1, 3, 3, 3, 3, 5, 2, 3, 3, 1])
```
"""
function histcountindices(x, xedges::AbstractRange; T=Int64)
    N = fill(zero(T), length(xedges)-1)
    bin = fill(0, size(x))
    histcountindices!(N, bin, x, xedges)
    return N, bin
end
histcountindices(x, xmin::Number, xmax::Number, nbins::Integer; T=Int64) = histcountindices(x, range(xmin, xmax, length=nbins+1); T=T)
export histcountindices


"""
```julia
histcounts(x, y, xedges::AbstractRange, yedges::AbstractRange; T=Int64)::Matrix{T}
```
A 2D histogram, ignoring `NaN`s: calculate the number of `x, y` pairs that fall into
each square of a 2D grid of equally-spaced square bins with edges specified by
`xedges` and `yedges`.

The resulting matrix `N` of counts is oriented with the lowest x and y bins in
`N[1,1]`, where the first (vertical / row) dimension of `N` corresponds to the y
axis (with `size(N,1) == length(yedges)-1`) and the second (horizontal / column)
dimension of `N` corresponds to the x axis (with `size(N,2) == length(xedges)-1`).

By default, the counts are returned as `Int64`s, though this can be changed by
specifying an output type with the optional keyword argument `T`.

## Examples
```julia
julia> x = y = 0.5:9.5;

julia> xedges = yedges = 0:10;

julia> N = histcounts(x,y,xedges,yedges)
10×10 Matrix{Int64}:
 1  0  0  0  0  0  0  0  0  0
 0  1  0  0  0  0  0  0  0  0
 0  0  1  0  0  0  0  0  0  0
 0  0  0  1  0  0  0  0  0  0
 0  0  0  0  1  0  0  0  0  0
 0  0  0  0  0  1  0  0  0  0
 0  0  0  0  0  0  1  0  0  0
 0  0  0  0  0  0  0  1  0  0
 0  0  0  0  0  0  0  0  1  0
 0  0  0  0  0  0  0  0  0  1
```
"""
function histcounts(x, y, xedges::AbstractRange, yedges::AbstractRange; T=Int64)
    N = fill(zero(T), length(yedges)-1, length(xedges)-1)
    histcounts!(N, x, y, xedges, yedges)
    return N
end
export histcounts

"""
```julia
histcounts!(N, x, xedges::AbstractRange)
```
Simple 1D histogram; as `histcounts`, but in-place, adding counts to the first
`length(xedges)-1` elements of Array `N`.

Note that counts will be added to `N`, not overwrite `N`, allowing you to produce
cumulative histograms. However, this means you will have to initialize `N` with
zeros before first use.
"""
function histcounts!(N::AbstractArray, x::AbstractArray, xedges::AbstractRange)
    @assert firstindex(N) === 1

    # What is the size of each bin?
    nbins = length(xedges) - 1
    xmin, xmax = extrema(xedges)
    δ𝑖δx = nbins/(xmax-xmin)

    # Make sure we don't have a segfault by filling beyond the length of N
    # in the @inbounds loop below
    if length(N) < nbins
        nbins = length(N)
        @warn "length(N) < nbins; any bins beyond length(N) will not be filled"
    end

    # Loop through each element of x
    @inbounds for n ∈ eachindex(x)
        xᵢ = x[n]
        𝑖 = (xᵢ - xmin) * δ𝑖δx
        if 0 < 𝑖 <= nbins
            i = ceil(Int, 𝑖)
            N[i] += 1
        end
    end
    return N
end


"""
```julia
histcountindices!(N, bin, x, xedges::AbstractRange)
```
Simple 1D histogram; as `histcounts!`, but also recording the bin index of each
`x` value.
"""
function histcountindices!(N::AbstractArray, bin::AbstractArray, x::AbstractArray, xedges::AbstractRange)
    @assert firstindex(N) === 1
    @assert firstindex(bin) === 1

    # What is the size of each bin?
    nbins = length(xedges) - 1
    xmin, xmax = extrema(xedges)
    δ𝑖δx = nbins/(xmax-xmin)

    # Make sure we don't have a segfault by filling beyond the length of N
    # in the @inbounds loop below
    if length(N) < nbins
        nbins = length(N)
        @warn "length(N) < nbins; any bins beyond length(N) will not be filled"
    end

    # Loop through each element of x
    @inbounds for n ∈ eachindex(x)
        xᵢ = x[n]
        𝑖 = (xᵢ - xmin) * δ𝑖δx
        if 0 < 𝑖 <= nbins
            i = ceil(Int, 𝑖)
            N[i] += 1
            bin[n] = i
        end
    end
    return N, bin
end

"""
```julia
histcounts!(N, x, y, xedges::AbstractRange, yedges::AbstractRange)
```
Simple 2D histogram; as `histcounts`, but in-place, adding counts to the first
`length(xedges)-1` columns and the first `length(yedges)-1` rows of `N`
elements of Array `N`.

Note that counts will be added to `N`, not overwrite `N`, allowing you to produce
cumulative histograms. However, this means you will have to initialize `N` with
zeros before first use.
"""
function histcounts!(N::AbstractMatrix, x::AbstractVector, y::AbstractVector, xedges::AbstractRange, yedges::AbstractRange)
    @assert firstindex(N) === 1

    # Calculate bin index from x value
    nxbins = length(xedges)-1
    xmin, xmax = extrema(xedges)
    δ𝑗δx = nxbins / (xmax - xmin)

    nybins = length(yedges)-1
    ymin, ymax = extrema(yedges)
    δ𝑖δy = nybins / (ymax - ymin)

    # Make sure we don't have a segfault by filling beyond the length of N
    # in the @inbounds loop below
    if size(N,1) < nybins
        nybins = size(N,1)
        @warn "size(N,1) < nybins; any y bins beyond size(N,1) will not be filled"
    end
    if size(N,2) < nxbins
        nxbins = size(N,2)
        @warn "size(N,2) < nxbins; any x bins beyond size(N,2) will not be filled"
    end

    # Calculate the means for each bin, ignoring NaNs
    @inbounds for n ∈ eachindex(x)
        𝑖 = (y[n] - ymin) * δ𝑖δy
        𝑗 = (x[n] - xmin) * δ𝑗δx
        if (0 < 𝑖 <= nybins) && (0 < 𝑗 <= nxbins)
            i = ceil(Int, 𝑖)
            j = ceil(Int, 𝑗)
            N[i,j] += 1
        end
    end

    return N
end
export histcounts!
