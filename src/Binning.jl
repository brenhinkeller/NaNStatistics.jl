## -- Binning functions, in-place

    """
    ```julia
    nanbinmean!(MU, [N], x, y, xedges::AbstractRange)
    ```
    Ignoring `NaN`s, fill the array `MU` with the means (and optionally `N` with
    the counts) of non-NAN `y` values that fall into each of `length(xedges)-1`
    equally spaced bins along the `x` axis with bin edges specified by `xedges`.

    Bins are open on the left and closed on the right, i.e., 
    (`xedges[i]`, `exedges[i+1]`] for the `i`th x bin.

    The array of `x` data should given as a one-dimensional array (any subtype
    of AbstractVector) and `y` as either a 1-d or 2-d array (any subtype of
    AbstractVecOrMat).

    The output arrays `MU` and `N` must be the same size, and must have the same
    number of columns as `y`; if `y` is a 2-d array (matrix), then each column of
    `y` will be treated as a separate variable.
    """
    function nanbinmean!(MU::AbstractVecOrMat, x::AbstractVector, y::AbstractVecOrMat, xedges::AbstractRange)
        N = Array{Int}(undef, size(MU))
        return nanbinmean!(MU, N, x, y, xedges)
    end
    # As above, but also return an array of counts, N; y, N, and MU as 1D vectors
    function nanbinmean!(MU::AbstractVector, N::AbstractVector, x::AbstractVector, y::AbstractVector, xedges::AbstractRange)
        # Check inputs
        @assert eachindex(x) == axes(y,1) "Indices of `x` and `y` must match, got $(eachindex(x)), $(axes(y,1))"
        @assert eachindex(MU) == eachindex(N) "Indices axis of `MU` and `N` must match, got $(eachindex(MU)), $(eachindex(N))"
        imin, imax = firstindex(MU)::Int-1, lastindex(MU)

        # What is the size of each bin?
        nbins = length(xedges) - 1
        @assert length(MU) == nbins "`MU`` has $(length(MU)) bins but `xedges` has $nbins bins"
        xmin, xmax = extrema(xedges)
        δ𝑖δx = nbins / (xmax - xmin)

        # Calculate the means for each bin, ignoring NaNs
        fill!(N, 0)
        fill!(MU, zero(eltype(MU))) # Fill the output array with zeros to start
        @inbounds for n ∈ eachindex(x)
            𝑖 = (x[n] - xmin) * δ𝑖δx + imin
            yₙ = y[n]
            if (yₙ==yₙ) && (imin < 𝑖 <= imax)
                i = ceil(Int, 𝑖)
                N[i] += 1
                MU[i] += yₙ
            end
        end
        MU ./= N # Divide by N to calculate means. Empty bin = 0/0 = NaN

        return MU
    end
    # In-place binning; as above but with y, N, MU as 2D matrices instead of 1D vectors
    function nanbinmean!(MU::AbstractMatrix, N::AbstractMatrix, x::AbstractVector, y::AbstractMatrix, xedges::AbstractRange)
        # Check inputs
        @assert eachindex(x) == axes(y,1) "Indices of `x` and `y` must match, got $(eachindex(x)), $(axes(y,1))"
        @assert axes(MU,1) == axes(N,1) "First axis of `MU` and `N` must match, got $(axes(MU,1)), $(axes(N,1))"
        @assert axes(MU,2) == axes(N,2) == axes(y,2) "Second axes of `MU`, `N`, and `y` must match, got $(axes(MU,2)), $(axes(N,2)), $(axes(y,2))"
        imin, imax = firstindex(MU,1)::Int-1, lastindex(MU,1)

        # What is the size of each bin?
        nbins = length(xedges) - 1
        @assert size(MU,1) == nbins "`MU`` has $(size(MU,1)) bins but `xedges` has $nbins bins"
        xmin, xmax = extrema(xedges)
        δ𝑖δx = nbins / (xmax - xmin)

        # Calculate the means for each bin, ignoring NaNs
        fill!(N, 0)
        fill!(MU, zero(eltype(MU))) # Fill the output array with zeros to start
        @inbounds for n ∈ eachindex(x)
            𝑖 = (x[n] - xmin) * δ𝑖δx + imin
            if (imin < 𝑖 <= imax)
                i = ceil(Int, 𝑖)
                for j ∈ axes(y,2)
                    if y[n,j]==y[n,j]
                        N[i,j] += 1
                        MU[i,j] += y[n,j]
                    end
                end
            end
        end
        MU ./= N # Divide by N to calculate means. Empty bin = 0/0 = NaN

        return MU
    end
    nanbinmean!(MU, x, y, xmin::Number, xmax::Number, nbins::Integer) = nanbinmean!(MU, x, y, range(xmin,xmax,length=nbins+1))
    nanbinmean!(MU, N, x, y, xmin::Number, xmax::Number, nbins::Integer) = nanbinmean!(MU, N, x, y, range(xmin,xmax,length=nbins+1))
    export nanbinmean!

    # 2-d binning
    """
    ```julia
    nanbinmean!(MU, N, x, y, z, xedges::AbstractRange, yedges::AbstractRange)
    ```
    Ignoring `NaN`s, fill the matrix `MU` with the means and `N` with
    the counts of non-NAN `z` values that fall into a 2D grid of x and y bins
    defined by `xedges` and `yedges`. The independent variables `x` and `y`,
    as well as the dependent variable `z`, are all expected as 1D vectors (any
    subtype of AbstractVector).

    Bins are open on the left and closed on the right, i.e., 
    (`xedges[i]`, `exedges[i+1]`] for the `i`th x bin.

    The output matrices `MU` and `N` must be the same size, and must each have
    `length(yedges)-1` rows and `length(xedges)-1` columns.
    """
    function nanbinmean!(MU::AbstractMatrix, N::AbstractMatrix, x::AbstractVector, y::AbstractVector, z::AbstractVector, xedges::AbstractRange, yedges::AbstractRange)
        # Check inputs
        @assert eachindex(x) == eachindex(y) == eachindex(z) "Indices of `x`, `y` and `z` must match, got $(eachindex(x)), $(eachindex(y)), $(eachindex(z))"
        @assert axes(MU) == axes(N) "Indices of `MU` and `N` must match, got $(axes(MU)), $(axes(N))"
        imin, imax = firstindex(MU,1)::Int-1, lastindex(MU,1)
        jmin, jmax = firstindex(MU,2)::Int-1, lastindex(MU,2)

        # What is the size of each bin?
        nxbins = length(xedges)-1
        @assert size(MU,2) == nxbins "`MU`` has $(size(MU,2)) x bins but `xedges` has $nxbins"
        xmin, xmax = extrema(xedges)
        δjδx = nxbins / (xmax - xmin)
        nybins = length(yedges)-1
        @assert size(MU,1) == nybins "`MU`` has $(size(MU,1)) y bins but `yedges` has $nybins"
        ymin, ymax = extrema(yedges)
        δiδy = nybins / (ymax - ymin)

        # Calculate the means for each bin, ignoring NaNs
        fill!(N, 0)
        fill!(MU, zero(eltype(MU))) # Fill the output array with zeros to start
        @inbounds for n ∈ eachindex(x)
            𝑖 = (y[n] - ymin) * δiδy + imin
            𝑗 = (x[n] - xmin) * δjδx + jmin
            zₙ = z[n]
            if (zₙ==zₙ) && (imin < 𝑖 <= imax) && (jmin < 𝑗 <= jmax)
                i = ceil(Int, 𝑖)
                j = ceil(Int, 𝑗)
                N[i,j] += 1
                MU[i,j] += zₙ
            end
        end
        MU ./= N # Divide by N to calculate means. Empty bin = 0/0 = NaN

        return MU
    end

    """
    ```julia
    nanbinmean!(MU, W, x, y, w, xedges::AbstractRange)
    ```
    Ignoring `NaN`s, fill the array `MU` with the weighted means (and `W` with
    the sum of weights) of non-NAN `y` values that fall into each of
    `length(xedges)-1` equally spaced bins along the `x` axis with bin edges
    specified by `xedges`.

    Bins are open on the left and closed on the right, i.e., 
    (`xedges[i]`, `exedges[i+1]`] for the `i`th x bin.

    The array of `x` data should given as a one-dimensional array (any subtype
    of AbstractVector) and `y` as either a 1-d or 2-d array (any subtype of
    AbstractVecOrMat).

    The output arrays `MU` and `N` must be the same size, and must have the same
    number of columns as `y`; if `y` is a 2-d array (matrix), then each column of
    `y` will be treated as a separate variable.
    """
    # In-place binning with weights; y, W, MU as 1D vectors
    function nanbinmean!(MU::AbstractVector, W::AbstractVector, x::AbstractVector, y::AbstractVector, w::AbstractVector, xedges::AbstractRange)
        # Check inputs
        @assert eachindex(x) == eachindex(y) == eachindex(w) "Indices of `x`, `y` and `w` must match, got $(eachindex(x)), $(eachindex(y)), $(eachindex(w))"
        @assert eachindex(MU) == eachindex(W) "Indices of `MU` and `W` must match, got $(axes(MU)), $(axes(W))"
        imin, imax = firstindex(MU)::Int-1, lastindex(MU)

        # What is the size of each bin?
        nbins = length(xedges) - 1
        @assert length(MU) == nbins "`MU`` has $(length(MU)) bins but `xedges` has $nbins bins"
        xmin, xmax = extrema(xedges)
        δ𝑖δx = nbins / (xmax - xmin)

        # Calculate the means for each bin, ignoring NaNs
        fill!(W, zero(eltype(W)))
        fill!(MU, zero(eltype(MU))) # Fill the output array with zeros to start
        @inbounds for n ∈ eachindex(x)
            𝑖 = (x[n] - xmin) * δ𝑖δx + imin
            yₙ = y[n]
            if (yₙ==yₙ) && (imin < 𝑖 <= imax)
                i = ceil(Int, 𝑖)
                W[i] += w[n]
                MU[i] += yₙ * w[n]
            end
        end
        MU ./= W # Divide by sum of weights to calculate means. Empty bin = 0/0 = NaN

        return MU
    end
    # In-place binning with weights; y, W, MU as 2D matrices
    function nanbinmean!(MU::AbstractMatrix, W::AbstractMatrix, x::AbstractVector, y::AbstractMatrix, w::AbstractVector, xedges::AbstractRange)
        # Check inputs
        @assert eachindex(x) == eachindex(w) == axes(y,1) "Indices of `x` and `w` must match with first axis of `y`, got $(eachindex(x)), $(eachindex(w)), $(axes(y,1))"
        @assert axes(MU,1) == axes(W,1) "First axis of `MU` and `W` must match, got $(axes(MU,1)), $(axes(W,1))"
        @assert axes(MU,2) == axes(W,2) == axes(y,2) "Second axes of `MU`, `N`, and `y` must match, got $(axes(MU,2)), $(axes(W,2)), $(axes(y,2))"
        imin, imax = firstindex(MU,1)::Int-1, lastindex(MU,1)
        
        # What is the size of each bin?
        nbins = length(xedges) - 1
        @assert size(MU,1) == nbins "`MU`` has $(size(MU,1)) bins but `xedges` has $nbins bins"
        xmin, xmax = extrema(xedges)
        δ𝑖δx = nbins / (xmax - xmin)

        # Calculate the means for each bin, ignoring NaNs
        fill!(W, zero(eltype(W)))
        fill!(MU, zero(eltype(MU))) # Fill the output array with zeros to start
        @inbounds for n ∈ eachindex(x)
            𝑖 = (x[n] - xmin) * δ𝑖δx + imin
            if (imin < 𝑖 <= imax)
                i = ceil(Int, 𝑖)
                for j ∈ axes(y,2)
                    if y[n,j]==y[n,j]
                        W[i,j] += w[n]
                        MU[i,j] += y[n,j]*w[n]
                    end
                end
            end
        end
        MU ./= W # Divide by sum of weights to calculate means. Empty bin = 0/0 = NaN

        return MU
    end
    nanbinmean!(MU, W, x, y, w, xmin::Number, xmax::Number, nbins::Integer) = nanbinmean!(MU, W, x, y, w, range(xmin,xmax,length=nbins+1))


    """
    ```julia
    nanbinmedian!(M, [N], x, y, xedges::AbstractRange)
    ```
    Fill the array `M` with the medians (and optionally `N` with the counts) of
    non-NaN `y` values that fall into each of `length(xedges)-1` equally spaced
    bins along the `x` axis with bin edges specified by `xedges`.

    If `y` is a 2-d array (matrix), each column will be treated as a separate variable
    """
    function nanbinmedian!(M::AbstractVector, x::AbstractVector, y::AbstractVector, xedges::AbstractRange)
        # Check inputs
        @assert eachindex(x) == eachindex(y) "Indices of `x` and `y` must match, got $(eachindex(x)) and $(eachindex(y))"
        @assert eachindex(M) == eachindex(xedges)[1:end-1] "Indices of bins and bin edges must be compatible"
        # Median for each bin
        t = BitVector(undef, length(x))
        @inbounds for i = eachindex(M)
            t .= (xedges[i] .< x .<= xedges[i+1])
            M[i] = _nanmedian!(y[t],:)
        end
        return M
    end
    function nanbinmedian!(M::AbstractMatrix, x::AbstractVector, y::AbstractMatrix, xedges::AbstractRange)
        # Check inputs
        @assert eachindex(x) == axes(y,1) "First axis of `x` and `y` must match, got $(eachindex(x)) and $(axes(y,1))"
        @assert axes(M,1) == eachindex(xedges)[1:end-1] "Indices of bins and bin edges must be compatible"
        @assert axes(M,2) == axes(y,2) "Second axis of `M` and `y` must match, got $(axes(M,2)) and $(axes(y,2))"
        # Median for each bin
        nbins = length(xedges) - 1
        t = BitVector(undef, length(x))
        @inbounds for i ∈ axes(M,1)
            t .= xedges[i] .< x .<= xedges[i+1]
            if any(t)
                for j ∈ axes(y,2)
                    M[i,j] = _nanmedian!(y[t,j],:)
                end
            else
                M[i,:] .= float(eltype(y))(NaN)
            end
        end
        return M
    end
    function nanbinmedian!(M::AbstractVector, N::AbstractVector, x::AbstractVector, y::AbstractVector, xedges::AbstractRange)
        # Check inputs
        @assert eachindex(x) == eachindex(y) "Indices of `x` and `y` must match, got $(eachindex(x)) and $(eachindex(y))"
        @assert eachindex(M) == eachindex(N) == eachindex(xedges)[1:end-1] "Indices of bins and bin edges must be compatible"
        # Median for each bin
        t = BitVector(undef, length(x))
        @inbounds for i ∈ eachindex(M)
            t .= (xedges[i] .< x .<= xedges[i+1])
            yₜ = y[t]
            M[i] = _nanmedian!(yₜ,:)
            N[i] = countnotnans(yₜ)
        end
        return M
    end
    function nanbinmedian!(M::AbstractMatrix, N::AbstractMatrix, x::AbstractVector, y::AbstractMatrix, xedges::AbstractRange)
        # Check inputs
        @assert eachindex(x) == axes(y,1) "First axis of `x` and `y` must match, got $(eachindex(x)) and $(axes(y,1))"
        @assert axes(M,1) == axes(N,1) == eachindex(xedges)[1:end-1] "Indices of bins and bin edges must be compatible"
        @assert axes(M,2) == axes(N,2) == axes(y,2) "Second axis of `M`, `N`, and `y` must match, got $(axes(M,2)), $(axes(N,2)), and $(axes(y,2))"
        # Median for each bin
        t = BitVector(undef, length(x))
        @inbounds for i ∈ axes(M,1)
            t .= xedges[i] .< x .<= xedges[i+1]
            if any(t)
                for j ∈ axes(y,2)
                    yₜ = y[t,j]
                    M[i,j] = _nanmedian!(yₜ,:)
                    N[i,j] = countnotnans(yₜ)
                end
            else
                M[i,:] .= NaN
                N[i,:] .= 0
            end
        end
        return M
    end
    nanbinmedian!(MU, x, y, xmin::Number, xmax::Number, nbins::Integer) = nanbinmedian!(MU, x, y, range(xmin,xmax,length=nbins+1))
    nanbinmedian!(MU, N, x, y, xmin::Number, xmax::Number, nbins::Integer) = nanbinmedian!(MU, N, x, y, range(xmin,xmax,length=nbins+1))
    export nanbinmedian!


## -- Binning functions, convenience versions (not in place)

    """
    ```julia
    nanbinmean(x, y, xedges::AbstractRange)
    ```
    Ignoring `NaN`s, calculate the mean of `y` values that fall into each of
    `length(xedges)-1` equally spaced bins along the `x` axis with bin edges
    specified by `xedges`.

    Bins are open on the left and closed on the right, i.e., 
    (`xedges[i]`, `exedges[i+1]`] for the `i`th x bin.

    The array of `x` data should be given as a one-dimensional array (any subtype
    of AbstractVector) and `y` as either a 1-d or 2-d array (any subtype of
    AbstractVecOrMat). If `y` is a 2-d array, then each column of `y` will be
    treated as a separate variable.

    ## Examples
    ```julia
    julia> nanbinmean([1:100..., 1], [1:100..., NaN], 0:25:100)
    4-element Vector{Float64}:
     13.0
     38.0
     63.0
     88.0

    julia> nanbinmean(1:100, reshape(1:300,100,3), 0:25:100)
    4×3 Matrix{Float64}:
     13.0  113.0  213.0
     38.0  138.0  238.0
     63.0  163.0  263.0
     88.0  188.0  288.0
    ```
    """
    function nanbinmean(x::AbstractVector, y::AbstractVecOrMat, xedges::AbstractRange)
        N = Array{Int}(undef, length(xedges)-1, size(y)[2:end]...)
        MU = Array{float(eltype(y))}(undef, length(xedges)-1, size(y)[2:end]...)
        return nanbinmean!(MU, N, x, y, xedges)
    end
    nanbinmean(x, y, xmin::Number, xmax::Number, nbins::Integer) = nanbinmean(x, y, range(xmin,xmax,length=nbins+1))
    export nanbinmean

    """
    ```julia
    nanbinmean(x, y, z, xedges, yedges)
    ```
    Ignoring `NaN`s, calculate the mean of `z` values that fall into a 2D grid of
    x and y bins with bin edges defined by `xedges` and `yedges`. The independent
    variables `x` and `y`, as well as the dependent variable `z`, are all expected
    as 1D vectors (any subtype of AbstractVector).

    Bins are open on the left and closed on the right, i.e., 
    (`xedges[i]`, `exedges[i+1]`] for the `i`th x bin.

    ## Examples
    ```julia
    julia> x = y = z = 0.5:9.5;

    julia> xedges = yedges = 0:10;

    julia> nanbinmean(x,y,z,xedges,yedges)
    10×10 Matrix{Float64}:
       0.5  NaN    NaN    NaN    NaN    NaN    NaN    NaN    NaN    NaN
     NaN      1.5  NaN    NaN    NaN    NaN    NaN    NaN    NaN    NaN
     NaN    NaN      2.5  NaN    NaN    NaN    NaN    NaN    NaN    NaN
     NaN    NaN    NaN      3.5  NaN    NaN    NaN    NaN    NaN    NaN
     NaN    NaN    NaN    NaN      4.5  NaN    NaN    NaN    NaN    NaN
     NaN    NaN    NaN    NaN    NaN      5.5  NaN    NaN    NaN    NaN
     NaN    NaN    NaN    NaN    NaN    NaN      6.5  NaN    NaN    NaN
     NaN    NaN    NaN    NaN    NaN    NaN    NaN      7.5  NaN    NaN
     NaN    NaN    NaN    NaN    NaN    NaN    NaN    NaN      8.5  NaN
     NaN    NaN    NaN    NaN    NaN    NaN    NaN    NaN    NaN      9.5
    ```
    """
    function nanbinmean(x::AbstractVector, y::AbstractVector, z::AbstractVector, xedges::AbstractRange, yedges::AbstractRange)
        N = Array{Int}(undef, length(yedges)-1, length(xedges)-1)
        MU = Array{float(eltype(z))}(undef, length(yedges)-1, length(xedges)-1)
        return nanbinmean!(MU, N, x, y, z, xedges, yedges)
    end

    """
    ```julia
    nanbinmean(x, y, xedges::AbstractRange)
    ```
    Ignoring `NaN`s, calculate the weighted mean of `y` values that
    fall into each of `length(xedges)-1` equally spaced bins along the `x`
    axis with bin edges specified by `xedges`.

    Bins are open on the left and closed on the right, i.e., 
    (`xedges[i]`, `exedges[i+1]`] for the `i`th x bin.

    The array of `x` data should given as a one-dimensional array (any subtype
    of AbstractVector) and `y` as either a 1-d or 2-d array (any subtype of
    AbstractVecOrMat). If `y` is a 2-d array, then each column of `y` will be
    treated as a separate variable.
    """
    function nanbinmean(x::AbstractVector, y::AbstractVecOrMat, w::AbstractVector, xedges::AbstractRange)
        W = Array{eltype(w)}(undef, length(xedges)-1, size(y)[2:end]...)
        MU = Array{float(eltype(y))}(undef, length(xedges)-1, size(y)[2:end]...)
        return nanbinmean!(MU, W, x, y, w, xedges)
    end
    nanbinmean(x, y, w, xmin::Number, xmax::Number, nbins::Integer) = nanbinmean(x, y, w, range(xmin,xmax,length=nbins+1))


    """
    ```julia
    nanbinmedian(x, y, xedges::AbstractRange)
    ```
    Calculate the median, ignoring `NaN`s, of y values that fall into each of
    `length(xedges)-1` equally spaced bins along the `x` axis with bin edges
    specified by `xedges`.

    If `y` is a 2-d array (matrix), each column will be treated as a separate variable

    ## Examples
    ```julia
    4-element Vector{Float64}:
     13.0
     38.0
     63.0
     88.0

    julia> nanbinmedian(1:100, reshape(1:300,100,3), 0:25:100)
    4×3 Matrix{Float64}:
     13.0  113.0  213.0
     38.0  138.0  238.0
     63.0  163.0  263.0
     88.0  188.0  288.0
    ```
    """
    function nanbinmedian(x::AbstractVector, y::AbstractVecOrMat, xedges::AbstractRange)
        M = Array{float(eltype(y))}(undef, length(xedges)-1, size(y)[2:end]...)
        return nanbinmedian!(M, x, y, xedges)
    end
    nanbinmedian(x, y, xmin::Number, xmax::Number, nbins::Integer) = nanbinmedian(x, y, range(xmin,xmax,length=nbins+1))
    export nanbinmedian


## -- End of File
