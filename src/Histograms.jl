
"""
```julia
histcounts(x, xedges::AbstractRange; T=Int64)::Vector{T}
```
A 1D histogram, ignoring NaNs: calculate the number of `x` values that fall into
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
histcounts(x, y, xedges::AbstractRange, yedges::AbstractRange; T=Int64)::Matrix{T}
```
A 2D histogram, ignoring NaNs: calculate the number of `x, y` pairs that fall into
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
function histcounts!(N::Array, x::AbstractArray, xedges::AbstractRange)
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
