## --- Transformations of arrays with NaNs

    function _reducedims(A, region)
        if (ndims(A) > 1) && (size(A,region)==1)
            return dropdims(A, dims=region)
        else
            return A
        end
    end


    """
    ```julia
    nanmask(A)
    ```
    Create a Boolean mask of dimensions `size(A)` that is false wherever `A` is `NaN`
    """
    nanmask(A) = nanmask!(Array{Bool}(undef,size(A)), A)
    export nanmask

    """
    ```julia
    nanmask!(mask, A)
    ```
    Fill a Boolean mask of dimensions `size(A)` that is false wherever `A` is `NaN`
    """
    function nanmask!(mask, A)
        @avx for i=1:length(A)
            mask[i] = A[i]==A[i]
        end
        return mask
    end
    # Special methods for arrays that cannot contain NaNs
    nanmask!(mask, A::AbstractArray{<:Integer}) = fill!(mask, true)
    nanmask!(mask, A::AbstractArray{<:Rational}) = fill!(mask, true)
    export nanmask!

    """
    ```julia
    zeronan!(A)
    ```
    Replace all `NaN`s in A with zeros of the same type
    """
    function zeronan!(A::Array)
        @avx for i ∈ eachindex(A)
            Aᵢ = A[i]
            A[i] = ifelse(Aᵢ==Aᵢ, Aᵢ, 0)
        end
        return A
    end
    export zeronan!


## --- Min & max ignoring NaNs

    """
    ```julia
    nanmax(a,b)
    ```
    As `max(a,b)`, but if either argument is `NaN`, return the other one
    """
    nanmax(a, b) = ifelse(a > b, a, b)
    nanmax(a, b::AbstractFloat) = ifelse(a==a, ifelse(b > a, b, a), b)
    nanmax(a, b::Vec{N,<:AbstractFloat}) where N = ifelse(a==a, ifelse(b > a, b, a), b)
    export nanmax

    """
    ```julia
    nanmin(a,b)
    ```
    As `min(a,b)`, but if either argument is `NaN`, return the other one
    """
    nanmin(a, b) = ifelse(a < b, a, b)
    nanmin(a, b::AbstractFloat) = ifelse(a==a, ifelse(b < a, b, a), b)
    nanmin(a, b::Vec{N,<:AbstractFloat}) where N = ifelse(a==a, ifelse(b < a, b, a), b)
    export nanmin

## --- Percentile statistics, excluding NaNs

    """
    ```julia
    nanpctile(A, p; dims
    ```
    Find the `p`th percentile of an indexable collection `A`, ignoring NaNs,
    optionally along a dimension specified by `dims`.

    A valid percentile value must satisfy 0 <= `p` <= 100.

    Also supports the `dim` keyword, which behaves identically to `dims`, but
    also drops any singleton dimensions that have been reduced over (as is the
    convention in some other languages).
    """
    nanpctile(A, p; dims=:, dim=:) = __nanpctile(A, p, dims, dim)
    __nanpctile(A, p, ::Colon, ::Colon) = _nanpctile(A, p, :)
    __nanpctile(A, p, region, ::Colon) = _nanpctile(A, p, region)
    __nanpctile(A, p, ::Colon, region) = _reducedims(_nanpctile(A, p, region), region)
    function _nanpctile(A, p, ::Colon)
        t = nanmask(A)
        return any(t) ? percentile(A[t],p) : NaN
    end
    function _nanpctile(A, p, region)
        s = size(A)
        if region == 2
            t = Array{Bool}(undef, s[2])
            result = Array{float(eltype(A))}(undef, s[1], 1)
            for i=1:s[1]
                nanmask!(t, A[i,:])
                result[i] = any(t) ? percentile(A[i,t],p) : NaN
            end
        elseif region == 1
            t = Array{Bool}(undef, s[1])
            result = Array{float(eltype(A))}(undef, 1, s[2])
            for i=1:s[2]
                nanmask!(t, A[:,i])
                result[i] = any(t) ? percentile(A[t,i],p) : NaN
            end
        else
            result = _nanpctile(A, p, :)
        end
        return result
    end
    export nanpctile


    """
    ```julia
    inpctile(A, p::Number; dims)
    ```
    Return a boolean array that identifies which values of the iterable
    collection `A` fall within the central `p`th percentile, optionally along a
    dimension specified by `dims`.

    A valid percentile value must satisfy 0 <= `p` <= 100.
    """
    function inpctile(A, p)
        offset = (100 - p) / 2
        return _nanpctile(A, offset, :) .< A .< _nanpctile(A, 100-offset, :)
    end
    export inpctile

## --- Combine arrays containing NaNs

    """
    ```julia
    nanadd(A, B)
    ```
    Add the non-NaN elements of A and B, treating NaNs as zeros
    """
    nanadd(a,b) = ifelse(a==a, a, zero(typeof(a))) + ifelse(b==b, b, zero(typeof(b)))
    function nanadd(A::AbstractArray, B::AbstractArray)
        result_type = promote_type(eltype(A), eltype(B))
        result = similar(A, result_type)
        @inbounds @simd for i ∈ eachindex(A)
            Aᵢ = A[i]
            Bᵢ = B[i]
            result[i] = (Aᵢ * (Aᵢ==Aᵢ)) + (Bᵢ * (Bᵢ==Bᵢ))
        end
        return result
    end
    export nanadd

    """
    ```julia
    nanadd!(A, B)
    ```
    Add the non-NaN elements of `B` to `A`, treating NaNs as zeros
    """
    function nanadd!(A::Array, B::AbstractArray)
        @inbounds @simd for i ∈ eachindex(A)
            Aᵢ = A[i]
            Bᵢ = B[i]
            A[i] = (Aᵢ * (Aᵢ==Aᵢ)) + (Bᵢ * (Bᵢ==Bᵢ))
        end
        return A
    end
    export nanadd!


## --- Summary statistics of arrays with NaNs

    """
    ```julia
    nanminimum(A; dims)
    ```
    As `minimum` but ignoring `NaN`s: Find the smallest non-`NaN` value of an
    indexable collection `A`, optionally along a dimension specified by `dims`.

    Also supports the `dim` keyword, which behaves identically to `dims`, but
    also drops any singleton dimensions that have been reduced over (as is the
    convention in some other languages).
    """
    nanminimum(A; dims=:, dim=:) = __nanminimum(A, dims, dim)
    __nanminimum(A, ::Colon, ::Colon) = _nanminimum(A, :)
    __nanminimum(A, region, ::Colon) = _nanminimum(A, region)
    __nanminimum(A, ::Colon, region) = _reducedims(_nanminimum(A, region), region)
    _nanminimum(A, region) = reduce(nanmin, A, dims=region, init=float(eltype(A))(NaN))
    _nanminimum(A::Array{<:Number}, ::Colon) = vreduce(nanmin, A)
    export nanminimum


    """
    ```julia
    nanmaximum(A; dims)
    ```
    Find the largest non-NaN value of an indexable collection `A`, optionally
    along a dimension specified by `dims`.

    Also supports the `dim` keyword, which behaves identically to `dims`, but
    also drops any singleton dimensions that have been reduced over (as is the
    convention in some other languages).
    """
    nanmaximum(A; dims=:, dim=:) = __nanmaximum(A, dims, dim)
    __nanmaximum(A, ::Colon, ::Colon) = _nanmaximum(A, :)
    __nanmaximum(A, region, ::Colon) = _nanmaximum(A, region)
    __nanmaximum(A, ::Colon, region) = _reducedims(_nanmaximum(A, region), region)
    _nanmaximum(A, region) = reduce(nanmax, A, dims=region, init=float(eltype(A))(NaN))
    _nanmaximum(A::Array{<:Number}, ::Colon) = vreduce(nanmax, A)
    export nanmaximum


    """
    ```julia
    nanextrema(A; dims)
    ```
    Find the extrema (maximum & minimum) of an indexable collection `A`,
    ignoring NaNs, optionally along a dimension specified by `dims`.

    Also supports the `dim` keyword, which behaves identically to `dims`, but
    also drops any singleton dimensions that have been reduced over (as is the
    convention in some other languages).
    """
    nanextrema(A; dims=:, dim=:) = __nanextrema(A, dims, dim)
    __nanextrema(A, ::Colon, ::Colon) = (_nanminimum(A, :), _nanmaximum(A, :))
    __nanextrema(A, region, ::Colon) = collect(zip(__nanminimum(A, region, :), __nanmaximum(A, region, :)))
    __nanextrema(A, ::Colon, region) = collect(zip(__nanminimum(A, :, region), __nanmaximum(A, :, region)))
    export nanextrema


    """
    ```julia
    nanrange(A; dims)
    ```
    Calculate the range (maximum - minimum) of an indexable collection `A`,
    ignoring NaNs, optionally along a dimension specified by `dims`.

    Also supports the `dim` keyword, which behaves identically to `dims`, but
    also drops any singleton dimensions that have been reduced over (as is the
    convention in some other languages).
    """
    nanrange(A; dims=:, dim=:) = __nanmaximum(A, dims, dim) - __nanminimum(A, dims, dim)
    export nanrange


    """
    ```julia
    nanmean(A, W; dims)
    ```
    Ignoring NaNs, calculate the weighted mean of an indexable
    collection `A`, optionally along dimensions specified by `dims`.

    Also supports the `dim` keyword, which behaves identically to `dims`, but
    also drops any singleton dimensions that have been reduced over (as is the
    convention in some other languages).
    """
    nanmean(A, W; dims=:, dim=:) = __nanmean(A, W, dims, dim)
    __nanmean(A, W, ::Colon, ::Colon) = _nanmean(A, W, :)
    __nanmean(A, W, region, ::Colon) = _nanmean(A, W, region)
    __nanmean(A, W, ::Colon, region) = _reducedims(_nanmean(A, W, region), region)
    function _nanmean(A, W, region)
        mask = nanmask(A)
        return sum(A.*W.*mask, dims=region) ./ sum(W.*mask, dims=region)
    end
    # Fallback method for non-Arrays
    function _nanmean(A, W, ::Colon)
        n = zero(eltype(W))
        m = zero(promote_type(eltype(W), eltype(A)))
        @inbounds @simd for i ∈ eachindex(A)
            Aᵢ = A[i]
            Wᵢ = W[i]
            t = Aᵢ == Aᵢ
            n += Wᵢ * t
            m += Wᵢ * Aᵢ * t
        end
        return m / n
    end
    # Can't have NaNs if array is all Integers
    function _nanmean(A::Array{<:Integer}, W, ::Colon)
        n = zero(eltype(W))
        m = zero(promote_type(eltype(W), eltype(A)))
        @avx for i ∈ eachindex(A)
            Wᵢ = W[i]
            n += Wᵢ
            m += Wᵢ * A[i]
        end
        return m / n
    end
    # Optimized AVX method for floats
    function _nanmean(A::AbstractArray{<:AbstractFloat}, W, ::Colon)
        T1 = eltype(W)
        T2 = promote_type(eltype(W), eltype(A))
        n = zero(T1)
        m = zero(T2)
        @avx for i ∈ eachindex(A)
            Aᵢ = A[i]
            Wᵢ = W[i]
            t = Aᵢ==Aᵢ
            n += ifelse(t, Wᵢ, zero(T1))
            m += ifelse(t, Wᵢ * Aᵢ, zero(T2))
        end
        return m / n
    end
    export nanmean


    """
    ```julia
    nanstd(A, W; dims)
    ```
    Calculate the weighted standard deviation, ignoring NaNs, of an
    indexable collection `A`, optionally along a dimension specified by `dims`.

    Also supports the `dim` keyword, which behaves identically to `dims`, but
    also drops any singleton dimensions that have been reduced over (as is the
    convention in some other languages).
    """
    nanstd(A, W; dims=:, dim=:) = __nanstd(A, W, dims, dim)
    __nanstd(A, W, ::Colon, ::Colon) = _nanstd(A, W, :)
    __nanstd(A, W, region, ::Colon) = _nanstd(A, W, region)
    __nanstd(A, W, ::Colon, region) = _reducedims(_nanstd(A, W, region), region)
    function _nanstd(A, W, region)
        mask = nanmask(A)
        n = sum(mask, dims=region)
        w = sum(W.*mask, dims=region)
        s = sum(A.*W.*mask, dims=region) ./ w
        d = A .- s # Subtract mean, using broadcasting
        @avx for i ∈ eachindex(d)
            dᵢ = d[i]
            d[i] = (dᵢ * dᵢ * W[i]) * mask[i]
        end
        s .= sum(d, dims=region)
        @avx for i ∈ eachindex(s)
            s[i] = sqrt((s[i] * n[i]) / (w[i] * (n[i] - 1)))
        end
        return s
    end
    function _nanstd(A, W, ::Colon)
        n = 0
        w = zero(eltype(W))
        m = zero(promote_type(eltype(W), eltype(A)))
        @inbounds @simd for i ∈ eachindex(A)
            Aᵢ = A[i]
            Wᵢ = W[i]
            t = Aᵢ == Aᵢ
            n += t
            w += Wᵢ * t
            m += Wᵢ * Aᵢ * t
        end
        mu = m / w
        s = zero(typeof(mu))
        @inbounds @simd for i ∈ eachindex(A)
            Aᵢ = A[i]
            d = Aᵢ - mu
            s += (d * d * W[i]) * (Aᵢ == Aᵢ) # Zero if Aᵢ is NaN
        end
        return sqrt(s / w * n / (n-1))
    end
    function _nanstd(A::AbstractArray{<:AbstractFloat}, W, ::Colon)
        n = 0
        Tw = eltype(W)
        Tm = promote_type(eltype(W), eltype(A))
        w = zero(Tw)
        m = zero(Tm)
        @avx for i ∈ eachindex(A)
            Aᵢ = A[i]
            Wᵢ = W[i]
            t = Aᵢ==Aᵢ
            n += t
            w += ifelse(t, Wᵢ,  zero(Tw))
            m += ifelse(t, Wᵢ * Aᵢ, zero(Tm))
        end
        mu = m / w
        Tmu = typeof(mu)
        s = zero(Tmu)
        @avx for i ∈ eachindex(A)
            Aᵢ = A[i]
            d = Aᵢ - mu
            s += ifelse(Aᵢ==Aᵢ, d * d * W[i], zero(Tmu))
        end
        return sqrt(s / w * n / (n-1))
    end
    export nanstd


    """
    ```julia
    nanmedian(A; dims)
    ```
    Calculate the median, ignoring NaNs, of an indexable collection `A`,
    optionally along a dimension specified by `dims`.

    Also supports the `dim` keyword, which behaves identically to `dims`, but
    also drops any singleton dimensions that have been reduced over (as is the
    convention in some other languages).
    """
    nanmedian(A; dims=:, dim=:) = __nanmedian(A, dims, dim)
    __nanmedian(A, ::Colon, ::Colon) = _nanmedian(A, :)
    __nanmedian(A, region, ::Colon) = _nanmedian(A, region)
    __nanmedian(A, ::Colon, region) = _reducedims(_nanmedian(A, region), region)
    function _nanmedian(A, ::Colon)
        t = nanmask(A)
        return any(t) ? median(A[t]) : float(eltype(A))(NaN)
    end
    function _nanmedian(A, region)
        s = size(A)
        if region == 2
            t = Array{Bool}(undef, s[2])
            result = Array{float(eltype(A))}(undef, s[1], 1)
            for i=1:s[1]
                nanmask!(t, A[i,:])
                result[i] = any(t) ? median(A[i,t]) : float(eltype(A))(NaN)
            end
        elseif region == 1
            t = Array{Bool}(undef, s[1])
            result = Array{float(eltype(A))}(undef, 1, s[2])
            for i=1:s[2]
                nanmask!(t, A[:,i])
                result[i] = any(t) ? median(A[t,i]) : float(eltype(A))(NaN)
            end
        else
            result = _nanmedian(A, :)
        end
        return result
    end
    export nanmedian


    """
    ```julia
    nanmad(A; dims)
    ```
    Median absolute deviation from the median, ignoring NaNs, of an indexable
    collection `A`, optionally along a dimension specified by `dims`.
    Note that for a Normal distribution, sigma = 1.4826 * MAD.

    Also supports the `dim` keyword, which behaves identically to `dims`, but
    also drops any singleton dimensions that have been reduced over (as is the
    convention in some other languages).
    """
    nanmad(A; dims=:, dim=:) = __nanmad(A, dims, dim)
    __nanmad(A, dims, dim) = __nanmedian(abs.(A .- _nanmedian(A, dims)), dims, dim)
    __nanmad(A, ::Colon, dim) = __nanmedian(abs.(A .- _nanmedian(A, dim)), :, dim)
    export nanmad


    """
    ```julia
    nanaad(A; dims)
    ```
    Mean (average) absolute deviation from the mean, ignoring NaNs, of an
    indexable collection `A`, optionally along a dimension specified by `dims`.
    Note that for a Normal distribution, sigma = 1.253 * AAD.

    Also supports the `dim` keyword, which behaves identically to `dims`, but
    also drops any singleton dimensions that have been reduced over (as is the
    convention in some other languages).
    """
    nanaad(A; dims=:, dim=:) = __nanaad(A, dims, dim)
    __nanaad(A, dims, dim) = __nanmean(abs.(A .- _nanmean(A, dims)), dims, dim)
    __nanaad(A, ::Colon, dim) = __nanmean(abs.(A .- _nanmean(A, dim)), :, dim)
    export nanaad


## --- Normalize and nanstandardize arrays, ignoring NaNs

    """
    ```julia
    nanstandardize!(A::Array{<:AbstractFloat}; dims)
    ```
    Rescale `A` to unit variance and zero mean
    i.e. `A .= (A .- nanmean(A)) ./ nanstd(A)`
    """
    nanstandardize!(A::Array{<:AbstractFloat}; dims=:) = _nanstandardize!(A, dims)
    function _nanstandardize!(A::Array{<:AbstractFloat}, dims=:)
        μ = nanmean(A, dims=dims)
        σ = nanstd(A, dims=dims, mean=μ)
        A .-= μ
        A ./= σ
        return A
    end
    export nanstandardize!

    """
    ```julia
    nanstandardize(A; dims)
    ```
    Rescale a copy of `A` to unit variance and zero mean
    i.e. `(A .- nanmean(A)) ./ nanstd(A)`
    """
    nanstandardize(A::AbstractArray; dims=:) = _nanstandardize!(float.(A), dims)
    export nanstandardize


## -- Moving average, ignoring NaNs

    """
    ```julia
    movmean(x::AbstractVecOrMat, n::Number)
    ```
    Simple moving average of `x` in 1 or 2 dimensions, spanning `n` bins
    (or n*n in 2D), returning an array of the same size as `x`.

    For the resulting moving average to be symmetric, `n` must be an odd integer;
    if `n` is not an odd integer, the first odd integer greater than `n` will be
    used instead.
    """
    function movmean(x::AbstractVector, n::Number)
        mean_type = Base.promote_op(/, eltype(x), Int64)
        m = Array{mean_type}(undef, size(x))
        δ = ceil(Int, (n-1)/2)
        @inbounds for i ∈ eachindex(x)
            iₗ = max(i-δ, 1)
            iᵤ = min(i+δ, length(x))
            m[i] = nanmean(view(x, iₗ:iᵤ))
        end
        return m
    end
    function movmean(x::AbstractMatrix, n::Number)
        mean_type = Base.promote_op(/, eltype(x), Int64)
        m = Array{mean_type}(undef, size(x))
        δ = ceil(Int, (n-1)/2)
        𝐼 = repeat(1:size(x,1), 1, size(x,2))
        𝐽 = repeat((1:size(x,2))', size(x,1), 1)
        @inbounds for k ∈ eachindex(x)
            i = 𝐼[k]
            iₗ = max(i-δ, 1)
            iᵤ = min(i+δ, size(x,1))
            j = 𝐽[k]
            jₗ = max(j-δ, 1)
            jᵤ = min(j+δ, size(x,2))
            m[i,j] = nanmean(view(x, iₗ:iᵤ, jₗ:jᵤ))
        end
        return m
    end
    export movmean


## --- End of File
