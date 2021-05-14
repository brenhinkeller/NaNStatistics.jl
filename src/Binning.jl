## -- Binning functions, in-place

    """
    ```julia
    nanmean!(MU, [N], x, y, [w], xmin::Number, xmax::Number, nbins::Integer)
    ```
    Ignoring NaNs, fill the array `MU` with the means (and optionally `N` with
    the counts) of non-NAN `y` values that fall into each of `nbins` equally
    spaced `x` bins between `xmin` and `xmax`, aligned with bin edges as
    `xmin`:(`xmax`-`xmin`)/`nbins`:`xmax`

    If an optional array of weights [`w`] is specified, then `N` is required, and
    will be filled with the sum of weights for each bin.

    The array of `x` data should given as a one-dimensional array (any subtype
    of AbstractVector) and `y` as either a 1-d or 2-d array (any subtype of
    AbstractVecOrMat).

    The output arrays `MU` and `N` must be the same size, and must have the same
    number of columns as `y`; if `y` is a 2-d array (matrix), then each column of
    `y` will be treated as a separate variable.
    """
    function nanmean!(MU::AbstractVecOrMat, x::AbstractVector, y::AbstractVecOrMat, xmin::Number, xmax::Number, nbins::Integer)
        N = Array{Int}(undef, size(MU))
        return nanmean!(MU, N, x, y, xmin, xmax, nbins)
    end
    # As above, but also return an array of counts, N; y, N, and MU as 1D vectors
    function nanmean!(MU::AbstractVector, N::AbstractVector, x::AbstractVector, y::AbstractVector, xmin::Number, xmax::Number, nbins::Integer)
        # Calculate bin index from x value
        scalefactor = nbins / (xmax - xmin)

        # Calculate the means for each bin, ignoring NaNs
        fill!(N, 0)
        fill!(MU, 0) # Fill the output array with zeros to start
        @inbounds for i = 1:length(x)
            bin_index_float = (x[i] - xmin) * scalefactor
            if (0 < bin_index_float < nbins) && !isnan(y[i])
                bin_index = ceil(Int, bin_index_float)
                N[bin_index] += 1
                MU[bin_index] += y[i]
            end
        end
        MU ./= N # Divide by N to calculate means. Empty bin = 0/0 = NaN

        return MU
    end
    # In-place binning; as above but with y, N, MU as 2D matrices instead of 1D vectors
    function nanmean!(MU::AbstractMatrix, N::AbstractMatrix, x::AbstractVector, y::AbstractMatrix, xmin::Number, xmax::Number, nbins::Integer)
        # Calculate bin index from x value
        scalefactor = nbins / (xmax - xmin)

        # Calculate the means for each bin, ignoring NaNs
        fill!(N, 0)
        fill!(MU, 0) # Fill the output array with zeros to start
        @inbounds for i = 1:length(x)
            bin_index_float = (x[i] - xmin) * scalefactor
            if (0 < bin_index_float < nbins)
                bin_index = ceil(Int, bin_index_float)
                for j = 1:size(y,2)
                    if !isnan(y[i,j])
                        N[bin_index,j] += 1
                        MU[bin_index,j] += y[i,j]
                    end
                end
            end
        end
        MU ./= N # Divide by N to calculate means. Empty bin = 0/0 = NaN

        return MU
    end
    # In-place binning with weights; y, W, MU as 1D vectors
    function nanmean!(MU::AbstractVector, W::AbstractVector, x::AbstractVector, y::AbstractVector, w::AbstractVector, xmin::Number, xmax::Number, nbins::Integer)
        # Calculate bin index from x value
        scalefactor = nbins / (xmax - xmin)

        # Calculate the means for each bin, ignoring NaNs
        fill!(W, 0)
        fill!(MU, 0) # Fill the output array with zeros to start
        @inbounds for i = 1:length(x)
            bin_index_float = (x[i] - xmin) * scalefactor
            if (0 < bin_index_float < nbins) && !isnan(y[i])
                bin_index = ceil(Int, bin_index_float)
                W[bin_index] += w[i]
                MU[bin_index] += y[i]*w[i]
            end
        end
        MU ./= W # Divide by sum of weights to calculate means. Empty bin = 0/0 = NaN

        return MU
    end
    # In-place binning with weights; y, W, MU as 2D matrices
    function nanmean!(MU::AbstractMatrix, W::AbstractMatrix, x::AbstractVector, y::AbstractMatrix, w::AbstractVector, xmin::Number, xmax::Number, nbins::Integer)
        # Calculate bin index from x value
        scalefactor = nbins / (xmax - xmin)

        # Calculate the means for each bin, ignoring NaNs
        fill!(W, 0)
        fill!(MU, 0) # Fill the output array with zeros to start
        @inbounds for i = 1:length(x)
            bin_index_float = (x[i] - xmin) * scalefactor
            if (0 < bin_index_float < nbins)
                bin_index = ceil(Int, bin_index_float)
                for j = 1:size(y,2)
                    if !isnan(y[i,j])
                        W[bin_index,j] += w[i]
                        MU[bin_index,j] += y[i,j]*w[i]
                    end
                end
            end
        end
        MU ./= W # Divide by sum of weights to calculate means. Empty bin = 0/0 = NaN

        return MU
    end
    export nanmean!


    """
    ```julia
    nanmedian!(M::AbstractVecOrMat, x::AbstractVector, y::AbstractVecOrMat, xmin::Number, xmax::Number, nbins::Integer)
    ```
    Fill the array `M` with the medians of non-NaN `y` values that fall into
    each of `nbins` equally spaced `x` bins between `xmin` and `xmax`, aligned
    with bin edges as `xmin`:(`xmax`-`xmin`)/`nbins`:`xmax`

    If `y` is a 2-d array (matrix), each column will be treated as a separate variable
    """
    function nanmedian!(M::AbstractVector, x::AbstractVector, y::AbstractVector, xmin::Number, xmax::Number, nbins::Integer)
        binedges = range(xmin, xmax, length=nbins+1)
        t = Array{Bool}(undef, length(x))
        for i = 1:nbins
            t .= (binedges[i] .<= x .< binedges[i+1]) .& (y.==y)
            M[i] = any(t) ? median(y[t]) : float(eltype(A))(NaN)
        end
        return M
    end
    function nanmedian!(M::AbstractMatrix, x::AbstractVector, y::AbstractMatrix, xmin::Number, xmax::Number, nbins::Integer)
        binedges = range(xmin, xmax, length=nbins+1)
        t = Array{Bool}(undef, length(x))
        tj = Array{Bool}(undef, length(x))
        for i = 1:nbins
            t .= binedges[i] .<= x .< binedges[i+1]
            for j = 1:size(y,2)
                tj .= t .& .!isnan.(y[:,j])
                M[i,j] = any(tj) ? median(y[tj,j]) : float(eltype(A))(NaN)
            end
        end
        return M
    end
    function nanmedian!(M::AbstractVector, N::AbstractVector, x::AbstractVector, y::AbstractVector, xmin::Number, xmax::Number, nbins::Integer)
        binedges = range(xmin, xmax, length=nbins+1)
        t = Array{Bool}(undef, length(x))
        for i = 1:nbins
            t .= (binedges[i] .<= x .< binedges[i+1]) .& (y.==y)
            M[i] = any(t) ? median(y[t]) : float(eltype(A))(NaN)
            N[i] = count(t)
        end
        return M
    end
    function nanmedian!(M::AbstractMatrix, N::AbstractMatrix, x::AbstractVector, y::AbstractMatrix, xmin::Number, xmax::Number, nbins::Integer)
        binedges = range(xmin, xmax, length=nbins+1)
        t = Array{Bool}(undef, length(x))
        tj = Array{Bool}(undef, length(x))
        for i = 1:nbins
            t .= binedges[i] .<= x .< binedges[i+1]
            for j = 1:size(y,2)
                tj .= t .& .!isnan.(y[:,j])
                M[i,j] = any(tj) ? median(y[tj,j]) : float(eltype(A))(NaN)
                N[i,j] = count(tj)
            end
        end
        return M
    end
    export nanmedian!


## -- Binning functions, convenience versions (not in place)

    """
    ```julia
    nanmean(x, y, [w], xmin::Number, xmax::Number, nbins::Integer)
    ```
    Ignoring NaNs, calculate the mean (optionally weighted) of `y` values that
    fall into each of `nbins` equally spaced `x` bins between `xmin` and `xmax`,
    aligned with bin edges as `xmin`:(`xmax`-`xmin`)/`nbins`:`xmax`

    The array of `x` data should given as a one-dimensional array (any subtype
    of AbstractVector) and `y` as either a 1-d or 2-d array (any subtype of
    AbstractVecOrMat). If `y` is a 2-d array, then each column of `y` will be
    treated as a separate variable.
    """
    function nanmean(x::AbstractVector, y::AbstractVecOrMat, xmin::Number, xmax::Number, nbins::Integer)
        N = Array{Int}(undef, nbins, size(y)[2:end]...)
        MU = Array{float(eltype(y))}(undef, nbins, size(y)[2:end]...)
        return nanmean!(MU, N, x, y, xmin, xmax, nbins)
    end
    function nanmean(x::AbstractVector, y::AbstractVecOrMat, w::AbstractVector, xmin::Number, xmax::Number, nbins::Integer)
        W = Array{float(eltype(y))}(undef, nbins, size(y)[2:end]...)
        MU = Array{float(eltype(y))}(undef, nbins, size(y)[2:end]...)
        return nanmean!(MU, W, x, y, w, xmin, xmax, nbins)
    end
    export nanmean

    """
    ```julia
    nanmedian(x::AbstractVector, y::AbstractVecOrMat, xmin::Number, xmax::Number, nbins::Integer)
    ```
    Calculate the median, ignoring NaNs, of y values that fall into each of `nbins`
    equally spaced `x` bins between `xmin` and `xmax`, aligned with bin edges as
    `xmin`:(`xmax`-`xmin`)/`nbins`:`xmax`

    If `y` is a 2-d array (matrix), each column will be treated as a separate variable
    """
    function nanmedian(x::AbstractVector, y::AbstractVecOrMat, xmin::Number, xmax::Number, nbins::Integer)
        M = Array{float(eltype(y))}(undef, nbins, size(y)[2:end]...)
        return nanmedian!(M, x, y, xmin, xmax, nbins)
    end
    export nanmedian


## -- End of File
