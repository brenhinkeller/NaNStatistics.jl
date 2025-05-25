## --- Transformations of arrays with NaNs

    # Dropdims if there are dims to be dropped
    reducedims(A, dims) = A
    reducedims(A::AbstractVector, dims) = A
    reducedims(A::AbstractArray, dims) = dropdims(A; dims)

    """
    ```julia
    countnans(A)
    ```
    Return the number of elements of `A` that are `NaN`s.
    """
    function countnans(A)
        n = 0
        @inbounds @simd ivdep for i ∈ eachindex(A)
            n += A[i]!=A[i]
        end
        return n
    end
    export countnans

    """
    ```julia
    countnonnans(A)
    ```
    Return the number of elements of `A` that are not `NaN`s.
    """
    function countnotnans(A)
        n = 0
        @inbounds @simd ivdep for i ∈ eachindex(A)
            n += A[i]==A[i]
        end
        return n
    end
    export countnotnans

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
        @inbounds @simd ivdep for i ∈ eachindex(A)
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
    function zeronan!(A::AbstractArray{T}) where T
        ∅ = zero(T)
        @inbounds @simd ivdep for i ∈ eachindex(A)
            Aᵢ = A[i]
            A[i] = ifelse(Aᵢ==Aᵢ, Aᵢ, ∅)
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
    nanmax(a, b) = ifelse(a==a, ifelse(b > a, b, a), b)
    export nanmax

    """
    ```julia
    nanmin(a,b)
    ```
    As `min(a,b)`, but if either argument is `NaN`, return the other one
    """
    nanmin(a, b) = ifelse(a==a, ifelse(b < a, b, a), b)
    export nanmin

## --- Percentile statistics, excluding NaNs

    """
    ```julia
    inpctile(A, p::Number; dims)
    ```
    Return a boolean array that identifies which values of the iterable
    collection `A` fall within the central `p`th percentile, optionally along a
    dimension specified by `dims`.

    A valid percentile value must satisfy `0 <= p <= 100`.
    """
    function inpctile(A::AbstractArray{T}, p::Number) where {T}
        q₊₋ = (100 - p) / 200
        Aₜ = copyto!(Array{T,1}(undef, length(A)), A)
        return _nanquantile!(Aₜ, q₊₋, :) .< A .< _nanquantile!(Aₜ, 1-q₊₋, :)
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
        @inbounds @simd ivdep for i ∈ eachindex(A,B)
            Aᵢ, Bᵢ = A[i], B[i]
            result[i] = (Aᵢ * (Aᵢ==Aᵢ)) + (Bᵢ * (Bᵢ==Bᵢ))
        end
        return result
    end
    export nanadd

    """
    ```julia
    nanadd!(A, B)
    ```
    Add the non-NaN elements of `B` to `A`, treating `NaN`s as zeros
    """
    function nanadd!(A::Array, B::AbstractArray)
        @inbounds @simd for i ∈ eachindex(A,B)
            Aᵢ, Bᵢ = A[i], B[i]
            A[i] = (Aᵢ * (Aᵢ==Aᵢ)) + (Bᵢ * (Bᵢ==Bᵢ))
        end
        return A
    end
    export nanadd!


## --- Summary statistics of arrays with NaNs

    uinit(T::Type) = uinit(float(T))
    uinit(T::Type{<:Integer}) = typemax(T)
    uinit(::Type{Float64}) = NaN
    uinit(::Type{Float32}) = NaN32
    uinit(::Type{Float16}) = NaN16
    linit(T::Type) = linit(float(T))
    linit(T::Type{<:Integer}) = typemin(T)
    linit(::Type{Float64}) = NaN
    linit(::Type{Float32}) = NaN32
    linit(::Type{Float16}) = NaN16

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
    __nanminimum(A, ::Colon, region) = reducedims(_nanminimum(A, region), region)
    _nanminimum(A, region) = reduce(nanmin, A, dims=region, init=uinit(eltype(A)))
    _nanminimum(A, ::Colon) = reduce(nanmin, A, init=uinit(eltype(A)))
    export nanminimum


    """
    ```julia
    nanmaximum(A; dims)
    ```
    As `maximum` but ignoring `NaN`s: Find the largest non-`NaN` value of an indexable collection `A`, optionally
    along a dimension specified by `dims`.

    Also supports the `dim` keyword, which behaves identically to `dims`, but
    also drops any singleton dimensions that have been reduced over (as is the
    convention in some other languages).
    """
    nanmaximum(A; dims=:, dim=:) = __nanmaximum(A, dims, dim)
    __nanmaximum(A, ::Colon, ::Colon) = _nanmaximum(A, :)
    __nanmaximum(A, region, ::Colon) = _nanmaximum(A, region)
    __nanmaximum(A, ::Colon, region) = reducedims(_nanmaximum(A, region), region)
    _nanmaximum(A, region) = reduce(nanmax, A, dims=region, init=linit(eltype(A)))
    _nanmaximum(A, ::Colon) = reduce(nanmax, A, init=linit(eltype(A)))
    export nanmaximum

    """
    ```julia
    nanargmin(A)
    ```
    As `argmin` but ignoring `NaN`s: Find the index of the smallest non-`NaN`
    value (if any) of an indexable collection `A`
    """
    function nanargmin(x)
        imin = firstindex(x)
        @inbounds for i in eachindex(x)
            xᵢ = x[i]
            if xᵢ==xᵢ && !(xᵢ > x[imin])
                imin = i
            end
        end
        return imin
    end
    export nanargmin

    """
    ```julia
    nanargmax(A)
    ```
    As `argmax` but ignoring `NaN`s: Find the index of the largest non-`NaN`
    value (if any) of an indexable collection `A`
    """
    function nanargmax(x)
        imax = firstindex(x)
        @inbounds for i in eachindex(x)
            xᵢ = x[i]
            if xᵢ==xᵢ && !(xᵢ < x[imax])
                imax = i
            end
        end
        return imax
    end
    export nanargmax

    """
    ```julia
    nanextrema(A; dims)
    ```
    Find the extrema (maximum & minimum) of an indexable collection `A`,
    ignoring `NaN`s, optionally along a dimension specified by `dims`.

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
    ignoring `NaN`s, optionally along a dimension specified by `dims`.

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
    Ignoring `NaN`s, calculate the weighted mean of an indexable
    collection `A`, optionally along dimensions specified by `dims`.

    Also supports the `dim` keyword, which behaves identically to `dims`, but
    also drops any singleton dimensions that have been reduced over (as is the
    convention in some other languages).
    """
    nanmean(A, W; dims=:, dim=:) = __nanmean_weighted(A, W, dims, dim)
    __nanmean_weighted(A, W, ::Colon, ::Colon) = _nanmean_weighted(A, W, :)
    __nanmean_weighted(A, W, region, ::Colon) = _nanmean_weighted(A, W, region)
    __nanmean_weighted(A, W, ::Colon, region) = reducedims(_nanmean_weighted(A, W, region), region)
    function _nanmean_weighted(A, W, region)
        mask = nanmask(A)
        return sum(A.*W.*mask, dims=region) ./ sum(W.*mask, dims=region)
    end
    # Can't have NaNs if array is all Integers
    function _nanmean_weighted(A::AbstractArray{<:Integer}, W, ::Colon)
        n = zero(eltype(W))
        m = zero(promote_type(eltype(W), eltype(A)))
        @inbounds @simd ivdep for i ∈ eachindex(A,W)
            Wᵢ = W[i]
            n += Wᵢ
            m += Wᵢ * A[i]
        end
        return m / n
    end
    function _nanmean_weighted(A, W, ::Colon)
        n = ∅ₙ = zero(eltype(W))
        m = ∅ₘ = zero(promote_type(eltype(W), eltype(A)))
        @inbounds @simd ivdep for i ∈ eachindex(A,W)
            Aᵢ, Wᵢ = A[i], W[i]
            t = Aᵢ==Aᵢ
            n += ifelse(t, Wᵢ, ∅ₙ)
            m += ifelse(t, Wᵢ * Aᵢ, ∅ₘ)
        end
        return m / n
    end
    export nanmean


    """
    ```julia
    nanstd(A, W; dims)
    ```
    Calculate the weighted standard deviation, ignoring `NaN`s, of an
    indexable collection `A`, optionally along a dimension specified by `dims`.

    Also supports the `dim` keyword, which behaves identically to `dims`, but
    also drops any singleton dimensions that have been reduced over (as is the
    convention in some other languages).
    """
    nanstd(A, W; dims=:, dim=:) = __nanstd(A, W, dims, dim)
    __nanstd(A, W, ::Colon, ::Colon) = _nanstd(A, W, :)
    __nanstd(A, W, region, ::Colon) = _nanstd(A, W, region)
    __nanstd(A, W, ::Colon, region) = reducedims(_nanstd(A, W, region), region)
    function _nanstd(A, W, region)
        mask = nanmask(A)
        n = sum(mask, dims=region)
        w = sum(W.*mask, dims=region)
        s = sum(A.*W.*mask, dims=region) ./ w
        d = A .- s # Subtract mean, using broadcasting
        @inbounds @simd ivdep for i ∈ eachindex(d, W)
            dᵢ = d[i]
            d[i] = (dᵢ * dᵢ * W[i]) * mask[i]
        end
        s .= sum(d, dims=region)
        @inbounds @simd ivdep for i ∈ eachindex(s,n,w)
            s[i] = sqrt((s[i] * n[i]) / (w[i] * (n[i] - 1)))
        end
        return s
    end
    function _nanstd(A, W, ::Colon)
        @assert eachindex(A) == eachindex(W)
        n = 0
        w = zero(eltype(W))
        m = zero(promote_type(eltype(W), eltype(A)))
        @inbounds @simd ivdep for i ∈ eachindex(A,W)
            Aᵢ, Wᵢ = A[i], W[i]
            t = Aᵢ == Aᵢ
            n += t
            w += Wᵢ * t
            m += Wᵢ * Aᵢ * t
        end
        mu = m / w
        s = zero(typeof(mu))
        @inbounds @simd ivdep for i ∈ eachindex(A,W)
            Aᵢ = A[i]
            d = Aᵢ - mu
            s += (d * d * W[i]) * (Aᵢ == Aᵢ) # Zero if Aᵢ is NaN
        end
        return sqrt(s / w * n / (n-1))
    end
    export nanstd


    """
    ```julia
    nanmad(A; dims)
    ```
    Median absolute deviation from the median, ignoring `NaN`s, of an indexable
    collection `A`, optionally along a dimension specified by `dims`.
    Note that for a Normal distribution, sigma = 1.4826 * MAD.

    Also supports the `dim` keyword, which behaves identically to `dims`, but
    also drops any singleton dimensions that have been reduced over (as is the
    convention in some other languages).
    """
    nanmad(A; dims=:, dim=:) = __nanmad(A, dims, dim)
    __nanmad(A::AbstractArray{N,T}, dims, dim) where {N,T} = __nanmad!(copyto!(Array{N,T}(undef, size(A)), A), dims, dim)
    __nanmad!(A, dims, dim) = __nanmedian!(abs.(A .- _nanmedian!(A, dims)), dims, dim)
    __nanmad!(A, ::Colon, dim) = __nanmedian!(abs.(A .- _nanmedian!(A, dim)), :, dim)
    export nanmad


    """
    ```julia
    nanmad!(A; dims)
    ```
    As `nanmad` but in-place.
    """
    nanmad!(A; dims=:, dim=:) = __nanmad!(A, dims, dim)
    export nanmad!


    """
    ```julia
    nanaad(A; dims, size_threshold=NANMEAN_SIZE_THRESHOLD)
    ```
    Mean (average) absolute deviation from the mean, ignoring `NaN`s, of an
    indexable collection `A`, optionally along a dimension specified by `dims`.
    Note that for a Normal distribution, sigma = 1.253 * AAD. The
    `size_threshold` argument is supported for taking the mean, see `nanmean`
    for more information.

    Also supports the `dim` keyword, which behaves identically to `dims`, but
    also drops any singleton dimensions that have been reduced over (as is the
    convention in some other languages).
    """
    nanaad(A; dims=:, dim=:, size_threshold=NANMEAN_SIZE_THRESHOLD) = __nanaad(A, dims, dim, size_threshold)
    __nanaad(A, dims, dim, size_threshold) = __nanmean(abs.(A .- _nanmean(A, dims, size_threshold)), dims, dim, size_threshold)
    __nanaad(A, ::Colon, dim, size_threshold) = __nanmean(abs.(A .- _nanmean(A, dim, size_threshold)), :, dim, size_threshold)
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


## --- Moving sum, ignoring NaNs

    """
    ```julia
    movsum(x::AbstractVecOrMat, n::Number)
    ```
    Simple moving sum of `x` in 1 or 2 dimensions, spanning `n` bins
    (or n*n in 2D), returning an array of the same size as `x`.

    For the resulting moving sum to be symmetric, `n` must be an odd integer;
    if `n` is not an odd integer, the first odd integer greater than `n` will be
    used instead.
    """
    function movsum(x::AbstractVector, n::Number)
        m = similar(x)
        δ = ceil(Int, (n-1)/2)
        @inbounds for i ∈ eachindex(x)
            iₗ = max(i-δ, firstindex(x))
            iᵤ = min(i+δ, lastindex(x))
            m[i] = nansum(view(x, iₗ:iᵤ))
        end
        return m
    end
    function movsum(x::AbstractMatrix, n::Number)
        m = similar(x)
        δ = ceil(Int, (n-1)/2)
        𝐼 = repeat((firstindex(x,1):lastindex(x,1)), 1, size(x,2))
        𝐽 = repeat((firstindex(x,2):lastindex(x,2))', size(x,1), 1)
        @inbounds for k ∈ eachindex(𝐼,𝐽)
            i = 𝐼[k]
            iₗ = max(i-δ, firstindex(x,1))
            iᵤ = min(i+δ, lastindex(x,1))
            j = 𝐽[k]
            jₗ = max(j-δ, firstindex(x,2))
            jᵤ = min(j+δ, lastindex(x,2))
            m[i,j] = nansum(view(x, iₗ:iᵤ, jₗ:jᵤ))
        end
        return m
    end
    export movsum

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
    function movmean(x::AbstractVector{T}, n::Number) where T
        mean_type = Base.promote_op(/, T, Int64)
        m = Array{mean_type}(undef, size(x))
        δ = ceil(Int, (n-1)/2)
        @inbounds for i ∈ eachindex(x)
            iₗ = max(i-δ, firstindex(x))
            iᵤ = min(i+δ, lastindex(x))
            m[i] = nanmean(view(x, iₗ:iᵤ))
        end
        return m
    end


    """
        movmean(x::AbstractVector{T}, win::Tuple{Int, Int}=(1, 1); skip_centre=false) where {T<:Real}
    
    Compute the simple moving average of a 1-dimensional array `x` over a window defined by `win`, returning an array of the same size as `x`.
    
    # Arguments
    - `x::AbstractVector{T}`: Input array of type `T`.
    - `win::Tuple{Int64, Int64}`: Tuple defining the window size to the left and right of each element. Default is `(1, 1)`.
    - `skip_centre::Bool`: If `true`, the center element is skipped in the average calculation. Default is `false`.
    
    # Returns
    - `Vector`: An array of the same size as `x` containing the moving averages.
    
    # Example
    ```julia
    x = [1, 2, 3, 4, 5]
    win = (1, 1)
    movmean(x, win)  # returns [1.5, 2.0, 3.0, 4.0, 4.5]
    ```
    """
    function movmean(x::AbstractVector{T}, win::Tuple{Int, Int}=(1, 1); 
        skip_centre=false) where {T<:Real}
        win_left, win_right = win
        
        FT = Base.promote_op(/, T, Int64)
        z = similar(x, FT)
        ∑ = ∅ = FT(0)
        ∑w = ∅w = 0

        @inbounds @simd for i ∈ eachindex(x)
            ibeg = max(i - win_left, firstindex(x))
            iend = min(i + win_right, lastindex(x))
            ∑ = ∅
            ∑w = ∅w
            for j = ibeg:iend
                skip = skip_centre && i==j
                xᵢ = x[j]
                notnan = (xᵢ == xᵢ) && !skip
                ∑ += ifelse(notnan, xᵢ, ∅)
                ∑w += ifelse(notnan, 1, 0)
            end
            z[i] = ∑ / ∑w
        end
        z
    end

    function movmean(x::AbstractMatrix{T}, n::Number) where T
        mean_type = Base.promote_op(/, T, Int64)
        m = Array{mean_type}(undef, size(x))
        δ = ceil(Int, (n-1)/2)
        𝐼 = repeat((firstindex(x,1):lastindex(x,1)), 1, size(x,2))
        𝐽 = repeat((firstindex(x,2):lastindex(x,2))', size(x,1), 1)
        @inbounds for k ∈ eachindex(𝐼,𝐽)
            i = 𝐼[k]
            iₗ = max(i-δ, firstindex(x,1))
            iᵤ = min(i+δ, lastindex(x,1))
            j = 𝐽[k]
            jₗ = max(j-δ, firstindex(x,2))
            jᵤ = min(j+δ, lastindex(x,2))
            m[i,j] = nanmean(view(x, iₗ:iᵤ, jₗ:jᵤ))
        end
        return m
    end
    export movmean

## --- Internal helpers

# This is a copy of the fallback branch of _mapreducedim!() in
# base/reducedim.jl. This iteration style is faster than the generated functions
# for large arrays. This is a macro instead of a function to guarantee that any
# overhead from passing an inner function is avoided, instead an expression is
# interpolated in.
macro _mapreduce_impl(B, A, expr)
    quote
        indsAt = Base.safe_tail(axes($A))
        indsRt = Base.safe_tail(axes($B))
        keep, Idefault = Broadcast.shapeindexer(indsRt)

        for IA in CartesianIndices(indsAt)
            IR = Broadcast.newindex(IA, keep, Idefault)
            @inbounds for ir in axes($A, 1)
                x = A[ir, IA]
                $expr
            end
        end
    end |> esc
end

## --- End of File
