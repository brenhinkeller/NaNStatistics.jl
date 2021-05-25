## -- Binning functions, in-place

    """
    ```julia
    nanbinmean!(MU, [N], x, y, xedges::AbstractRange)
    ```
    Ignoring NaNs, fill the array `MU` with the means (and optionally `N` with
    the counts) of non-NAN `y` values that fall into each of `length(xedges)-1`
    equally spaced bins along the `x` axis with bin edges specified by `xedges`.

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
        # What is the size of each bin?
        nbins = length(xedges) - 1
        xmin, xmax = extrema(xedges)
        δiδx = nbins / (xmax - xmin)

        # Make sure we don't have a segfault by filling beyond the length of N
        # in the @inbounds loop below
        if length(MU) < nbins
            nbins = length(MU)
            @warn "length(MU) < nbins; any bins beyond length(MU) will not be filled"
        end

        # Calculate the means for each bin, ignoring NaNs
        fill!(N, 0)
        fill!(MU, 0) # Fill the output array with zeros to start
        @inbounds for i = 1:length(x)
            bin_index_float = (x[i] - xmin) * δiδx
            if (0 < bin_index_float < nbins) && (y[i]==y[i])
                bin_index = ceil(Int, bin_index_float)
                N[bin_index] += 1
                MU[bin_index] += y[i]
            end
        end
        MU ./= N # Divide by N to calculate means. Empty bin = 0/0 = NaN

        return MU
    end
    # In-place binning; as above but with y, N, MU as 2D matrices instead of 1D vectors
    function nanbinmean!(MU::AbstractMatrix, N::AbstractMatrix, x::AbstractVector, y::AbstractMatrix, xedges::AbstractRange)
        # What is the size of each bin?
        nbins = length(xedges) - 1
        xmin, xmax = extrema(xedges)
        δiδx = nbins / (xmax - xmin)
        ncols = size(y,2)

        # Make sure we don't have a segfault by filling beyond the length of N
        # in the @inbounds loop below
        if size(MU,1) < nbins
            nbins = size(MU,1)
            @warn "size(MU,1) < nbins; any bins beyond size(MU,1) will not be filled"
        end
        if size(MU,2) < ncols
            ncols = size(MU,2)
            @warn "size(MU,2) < size(y,2); any columns beyond size(MU,2) will not be filled"
        end

        # Calculate the means for each bin, ignoring NaNs
        fill!(N, 0)
        fill!(MU, 0) # Fill the output array with zeros to start
        @inbounds for i = 1:length(x)
            bin_index_float = (x[i] - xmin) * δiδx
            if (0 < bin_index_float < nbins)
                bin_index = ceil(Int, bin_index_float)
                for j = 1:ncols
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
    nanbinmean!(MU, x, y, xmin::Number, xmax::Number, nbins::Integer) = nanbinmean!(MU, x, y, range(xmin,xmax,length=nbins+1))
    nanbinmean!(MU, N, x, y, xmin::Number, xmax::Number, nbins::Integer) = nanbinmean!(MU, N, x, y, range(xmin,xmax,length=nbins+1))
    export nanbinmean!

    # 2-d binning
    """
    ```julia
    nanbinmean!(MU, N, x, y, z, xedges::AbstractRange, yedges::AbstractRange)
    ```
    Ignoring NaNs, fill the matrix `MU` with the means and `N` with
    the counts of non-NAN `z` values that fall into a 2D grid of x and y bins
    defined by `xedges` and `yedges`. The independent variables `x` and `y`,
    as well as the dependent variable `z`, are all expected as 1D vectors (any
    subtype of AbstractVector).

    The output matrices `MU` and `N` must be the same size, and must each have
    `length(yedges)-1` rows and `length(xedges)-1` columns.
    """
    function nanbinmean!(MU::AbstractMatrix, N::AbstractMatrix, x::AbstractVector, y::AbstractVector, z::AbstractVector, xedges::AbstractRange, yedges::AbstractRange)
        # Calculate bin index from x value
        nxbins = length(xedges)-1
        xmin, xmax = extrema(xedges)
        δjδx = nxbins / (xmax - xmin)

        nybins = length(yedges)-1
        ymin, ymax = extrema(yedges)
        δiδy = nybins / (ymax - ymin)

        # Make sure we don't have a segfault by filling beyond the length of N
        # in the @inbounds loop below
        if size(MU, 1) < nybins
            nybins = size(MU, 1)
            @warn "size(MU, 1) < nybins; any y bins beyond size(MU, 1) will not be filled"
        end
        if size(MU, 2) < nxbins
            nxbins = size(MU, 2)
            @warn "size(MU, 2) < nxbins; any x bins beyond size(MU, 2) will not be filled"
        end

        # Calculate the means for each bin, ignoring NaNs
        fill!(N, 0)
        fill!(MU, 0) # Fill the output array with zeros to start
        @inbounds for n = 1:length(x)
            i_float = (y[n] - ymin) * δiδy
            j_float = (x[n] - xmin) * δjδx
            zₙ = z[n]
            if (zₙ==zₙ) && (0 < i_float < nybins) && (0 < j_float < nxbins)
                i = ceil(Int, i_float)
                j = ceil(Int, j_float)
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
    Ignoring NaNs, fill the array `MU` with the weighted means (and `W` with
    the sum of weights) of non-NAN `y` values that fall into each of
    `length(xedges)-1` equally spaced bins along the `x` axis with bin edges
    specified by `xedges`.

    The array of `x` data should given as a one-dimensional array (any subtype
    of AbstractVector) and `y` as either a 1-d or 2-d array (any subtype of
    AbstractVecOrMat).

    The output arrays `MU` and `N` must be the same size, and must have the same
    number of columns as `y`; if `y` is a 2-d array (matrix), then each column of
    `y` will be treated as a separate variable.
    """
    # In-place binning with weights; y, W, MU as 1D vectors
    function nanbinwmean!(MU::AbstractVector, W::AbstractVector, x::AbstractVector, y::AbstractVector, w::AbstractVector, xedges::AbstractRange)
        # What is the size of each bin?
        nbins = length(xedges) - 1
        xmin, xmax = extrema(xedges)
        δiδx = nbins / (xmax - xmin)

        # Make sure we don't have a segfault by filling beyond the length of N
        # in the @inbounds loop below
        if length(MU) < nbins
            nbins = length(MU)
            @warn "length(MU) < nbins; any bins beyond length(MU) will not be filled"
        end

        # Calculate the means for each bin, ignoring NaNs
        fill!(W, 0)
        fill!(MU, 0) # Fill the output array with zeros to start
        @inbounds for i = 1:length(x)
            bin_index_float = (x[i] - xmin) * δiδx
            if (0 < bin_index_float < nbins) && (y[i]==y[i])
                bin_index = ceil(Int, bin_index_float)
                W[bin_index] += w[i]
                MU[bin_index] += y[i]*w[i]
            end
        end
        MU ./= W # Divide by sum of weights to calculate means. Empty bin = 0/0 = NaN

        return MU
    end
    # In-place binning with weights; y, W, MU as 2D matrices
    function nanbinwmean!(MU::AbstractMatrix, W::AbstractMatrix, x::AbstractVector, y::AbstractMatrix, w::AbstractVector, xedges::AbstractRange)
        nbins = length(xedges) - 1
        xmin, xmax = extrema(xedges)
        δiδx = nbins / (xmax - xmin)
        ncols = size(y,2)

        # Make sure we don't have a segfault by filling beyond the length of N
        # in the @inbounds loop below
        if size(MU,1) < nbins
            nbins = size(MU,1)
            @warn "size(MU,1) < nbins; any bins beyond size(MU,1) will not be filled"
        end
        if size(MU,2) < ncols
            ncols = size(MU,2)
            @warn "size(MU,2) < size(y,2); any columns beyond size(MU,2) will not be filled"
        end

        # Calculate the means for each bin, ignoring NaNs
        fill!(W, 0)
        fill!(MU, 0) # Fill the output array with zeros to start
        @inbounds for i = 1:length(x)
            bin_index_float = (x[i] - xmin) * δiδx
            if (0 < bin_index_float < nbins)
                bin_index = ceil(Int, bin_index_float)
                for j = 1:ncols
                    if y[i,j]==y[i,j]
                        W[bin_index,j] += w[i]
                        MU[bin_index,j] += y[i,j]*w[i]
                    end
                end
            end
        end
        MU ./= W # Divide by sum of weights to calculate means. Empty bin = 0/0 = NaN

        return MU
    end
    nanbinwmean!(MU, W, x, y, w, xmin::Number, xmax::Number, nbins::Integer) = nanbinwmean!(MU, W, x, y, w, range(xmin,xmax,length=nbins+1))


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
        nbins = length(xedges) - 1
        t = Array{Bool}(undef, length(x))
        @inbounds for i = 1:nbins
            t .= (xedges[i] .<= x .< xedges[i+1]) .& (y.==y)
            M[i] = any(t) ? median(y[t]) : float(eltype(A))(NaN)
        end
        return M
    end
    function nanbinmedian!(M::AbstractMatrix, x::AbstractVector, y::AbstractMatrix, xedges::AbstractRange)
        nbins = length(xedges) - 1
        t = Array{Bool}(undef, length(x))
        tj = Array{Bool}(undef, length(x))
        @inbounds for i = 1:nbins
            t .= xedges[i] .<= x .< xedges[i+1]
            for j = 1:size(y,2)
                tj .= t .& .!isnan.(y[:,j])
                M[i,j] = any(tj) ? median(y[tj,j]) : float(eltype(A))(NaN)
            end
        end
        return M
    end
    function nanbinmedian!(M::AbstractVector, N::AbstractVector, x::AbstractVector, y::AbstractVector, xedges::AbstractRange)
        nbins = length(xedges) - 1
        t = Array{Bool}(undef, length(x))
        @inbounds for i = 1:nbins
            t .= (binedges[i] .<= x .< binedges[i+1]) .& (y.==y)
            M[i] = any(t) ? median(y[t]) : float(eltype(A))(NaN)
            N[i] = count(t)
        end
        return M
    end
    function nanbinmedian!(M::AbstractMatrix, N::AbstractMatrix, x::AbstractVector, y::AbstractMatrix, xedges::AbstractRange)
        nbins = length(xedges) - 1
        t = Array{Bool}(undef, length(x))
        tj = Array{Bool}(undef, length(x))
        @inbounds for i = 1:nbins
            t .= binedges[i] .<= x .< binedges[i+1]
            for j = 1:size(y,2)
                tj .= t .& .!isnan.(y[:,j])
                M[i,j] = any(tj) ? median(y[tj,j]) : float(eltype(A))(NaN)
                N[i,j] = count(tj)
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
    Ignoring NaNs, calculate the mean of `y` values that
    fall into each of `N = length(xedges)-1` equally spaced bins along the `x`
    axis with bin edges specified by `xedges`.

    The array of `x` data should given as a one-dimensional array (any subtype
    of AbstractVector) and `y` as either a 1-d or 2-d array (any subtype of
    AbstractVecOrMat). If `y` is a 2-d array, then each column of `y` will be
    treated as a separate variable.
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
    Ignoring NaNs, calculate the mean of `z` values that fall into a 2D grid of
    x and y bins with bin edges defined by `xedges` and `yedges`. The independent
    variables `x` and `y`, as well as the dependent variable `z`, are all expected
    as 1D vectors (any subtype of AbstractVector).
    """
    function nanbinmean(x::AbstractVector, y::AbstractVector, z::AbstractVector, xedges::AbstractRange, yedges::AbstractRange)
        N = Array{Int}(undef, length(yedges)-1, length(xedges)-1)
        MU = Array{float(eltype(y))}(undef, length(yedges)-1, length(xedges)-1)
        return nanbinmean!(MU, N, x, y, z, xedges, yedges)
    end

    """
    ```julia
    nanbinwmean(x, y, xedges::AbstractRange)
    ```
    Ignoring NaNs, calculate the weighted mean of `y` values that
    fall into each of `length(xedges)-1` equally spaced bins along the `x`
    axis with bin edges specified by `xedges`.

    The array of `x` data should given as a one-dimensional array (any subtype
    of AbstractVector) and `y` as either a 1-d or 2-d array (any subtype of
    AbstractVecOrMat). If `y` is a 2-d array, then each column of `y` will be
    treated as a separate variable.
    """
    function nanbinwmean(x::AbstractVector, y::AbstractVecOrMat, w::AbstractVector, xedges::AbstractRange)
        W = Array{float(eltype(y))}(undef, length(xedges)-1, size(y)[2:end]...)
        MU = Array{float(eltype(y))}(undef, length(xedges)-1, size(y)[2:end]...)
        return nanbinwmean!(MU, W, x, y, w, xedges)
    end
    nanbinwmean(x, y, w, xmin::Number, xmax::Number, nbins::Integer) = nanbinwmean(x, y, w, range(xmin,xmax,length=nbins+1))

    """
    ```julia
    nanbinmedian(x, y, xedges::AbstractRange)
    ```
    Calculate the median, ignoring NaNs, of y values that fall into each of
    `length(xedges)-1` equally spaced bins along the `x` axis with bin edges
    specified by `xedges`.

    If `y` is a 2-d array (matrix), each column will be treated as a separate variable
    """
    function nanbinmedian(x::AbstractVector, y::AbstractVecOrMat, xedges::AbstractRange)
        M = Array{float(eltype(y))}(undef, length(xedges)-1, size(y)[2:end]...)
        return nanbinmedian!(M, x, y, xedges)
    end
    nanbinmedian(x, y, xmin::Number, xmax::Number, nbins::Integer) = nanbinmedian(x, y, range(xmin,xmax,length=nbins+1))
    export nanbinmedian


## -- End of File
