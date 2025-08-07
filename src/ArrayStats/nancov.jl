function _nancov(x::AbstractVector, y::AbstractVector, corrected::Bool, μᵪ::Number, μᵧ::Number)
    # Calculate covariance
    σᵪᵧ = ∅ = zero(Base.promote_op(*, typeof(μᵪ), typeof(μᵧ)))
    n = 0
    @inbounds @simd ivdep for i ∈ eachindex(x,y)
            δᵪ = x[i] - μᵪ
            δᵧ = y[i] - μᵧ
            δ² = δᵪ * δᵧ
            notnan = δ²==δ²
            n += notnan
            σᵪᵧ += ifelse(notnan, δ², ∅)
    end
    σᵪᵧ = σᵪᵧ / (n-corrected)
    return σᵪᵧ
end

function _nancov(x::AbstractVector, y::AbstractVector, corrected::Bool)
    # Parwise nan-means
    n = 0
    Σᵪ = ∅ᵪ = zero(eltype(x))
    Σᵧ = ∅ᵧ = zero(eltype(y))
    @inbounds @simd ivdep for i ∈ eachindex(x,y)
        xᵢ, yᵢ = x[i], y[i]
        notnan = (xᵢ==xᵢ) & (yᵢ==yᵢ)
        n += notnan
        Σᵪ += ifelse(notnan, xᵢ, ∅ᵪ)
        Σᵧ += ifelse(notnan, yᵢ, ∅ᵧ)
    end

    return _nancov(x, y, corrected, Σᵪ/n, Σᵧ/n)
end



"""
```julia
nancov(x::AbstractVector, y::AbstractVector; corrected::Bool=true)
```
Compute the covariance between the vectors `x` and `y`.
As `Statistics.cov`, but ignoring `NaN`s.

If `corrected` is `true` as is the default, _Bessel's correction_ will be applied,
such that the sum is scaled by `n-1` rather than `n`, where `n = length(x)`.
"""
function nancov(x::AbstractVector, y::AbstractVector; corrected::Bool=true)
    @assert eachindex(x) == eachindex(y)
    return _nancov(x, y, corrected)
end

"""
```julia
nancov(X::AbstractMatrix; dims::Int=1, corrected::Bool=true)
```
Compute the covariance matrix of the matrix `X`, along dimension `dims`.
As `Statistics.cov`, but ignoring `NaN`s.

If `corrected` is `true` as is the default, _Bessel's correction_ will be applied,
such that the sum is scaled by `n-1` rather than `n`, where `n = length(x)`.
"""
function nancov(X::AbstractMatrix; dims::Int=1, corrected::Bool=true)
    Tₘ = Base.promote_op(/, eltype(X), Int)
    Tₒ = Base.promote_op(*, Tₘ, Tₘ)
    n = size(X, dims)
    m = size(X, mod(dims,2)+1)
    Σ = similar(X, Tₒ, (m, m))
    # Only two dimensions are possible, so handle each manually
    if dims == 1
        # Fill covariance matrix symmetrically
        @inbounds for i = 1:m
            for j = 1:i
                σᵢⱼ = _nancov(view(X,:,i), view(X,:,j), corrected)
                Σ[i,j] = Σ[j,i] = σᵢⱼ
            end
        end
    elseif dims == 2
        # Fill covariance matrix symmetrically
        @inbounds for i = 1:m
            for j = 1:i
                σᵢⱼ = _nancov(view(X,i,:), view(X,j,:), corrected)
                Σ[i,j] = Σ[j,i] = σᵢⱼ
            end
        end
    else
        throw("Dimension not in range")
    end
    return Σ
end
export nancov



function _nancovem(x::AbstractVector, y::AbstractVector, corrected::Bool)
    # Parwise nan-means
    n = 0
    Σᵪ = ∅ᵪ = zero(eltype(x))
    Σᵧ = ∅ᵧ = zero(eltype(y))
    @inbounds @simd ivdep for i ∈ eachindex(x,y)
        xᵢ, yᵢ = x[i], y[i]
        notnan = (xᵢ==xᵢ) & (yᵢ==yᵢ)
        n += notnan
        Σᵪ += ifelse(notnan, xᵢ, ∅ᵪ)
        Σᵧ += ifelse(notnan, yᵢ, ∅ᵧ)
    end
    μᵪ = Σᵪ/n
    μᵧ = Σᵧ/n

    return _nancov(x, y, corrected, μᵪ, μᵧ)/max(n,1)
end



"""
```julia
nancovem(x::AbstractVector, y::AbstractVector; corrected::Bool=true)
```
Compute the covariance of the error of the mean of the non-nan pairs of `x` and `y`,
where covariance of the error of the mean is to covariance 
as standard error of the mean is to standard deviation.

If `corrected` is `true` as is the default, _Bessel's correction_ will be applied,
to the covariance, such that the sum is scaled by `n-1` rather than `n`, where `n = length(x)`.
"""
function nancovem(x::AbstractVector, y::AbstractVector; corrected::Bool=true)
    @assert eachindex(x) == eachindex(y)
    return _nancovem(x, y, corrected)
end

"""
```julia
nancovem(X::AbstractMatrix; dims::Int=1, corrected::Bool=true)
```
Compute the covariance-of-error-of-the-mean matrix `X`, along dimension `dims`, ignoring NaNs.
That is, a matrix composed of the covariance of the error of the mean of
the non-nan pairs of each column (dims=1) or row (dims=2) of the matrix `X`,
where covariance of the error of the mean is to covariance 
as standard error of the mean is to standard deviation.


If `corrected` is `true` as is the default, _Bessel's correction_ will be applied,
such that the sum is scaled by `n-1` rather than `n`, where `n = length(x)`.
"""
function nancovem(X::AbstractMatrix; dims::Int=1, corrected::Bool=true)
    Tₘ = Base.promote_op(/, eltype(X), Int)
    Tₒ = Base.promote_op(*, Tₘ, Tₘ)
    n = size(X, dims)
    m = size(X, mod(dims,2)+1)
    Σ = similar(X, Tₒ, (m, m))
    # Only two dimensions are possible, so handle each manually
    if dims == 1
        # Fill covariance matrix symmetrically
        @inbounds for i = 1:m
            for j = 1:i
                σᵢⱼ = _nancovem(view(X,:,i), view(X,:,j), corrected)
                Σ[i,j] = Σ[j,i] = σᵢⱼ
            end
        end
    elseif dims == 2
        # Fill covariance matrix symmetrically
        @inbounds for i = 1:m
            for j = 1:i
                σᵢⱼ = _nancovem(view(X,i,:), view(X,j,:), corrected)
                Σ[i,j] = Σ[j,i] = σᵢⱼ
            end
        end
    else
        throw("Dimension not in range")
    end
    return Σ
end
export nancovem




"""
```julia
nancor(x::AbstractVector, y::AbstractVector)
```
Compute the (Pearson's product-moment) correlation between the vectors `x` and `y`.
As `Statistics.cor`, but ignoring `NaN`s.

Equivalent to `nancov(x,y) / (nanstd(x) * nanstd(y))`.
"""
function nancor(x::AbstractVector, y::AbstractVector; corrected::Bool=true)
    @assert eachindex(x) == eachindex(y)
    return _nancor(x, y, corrected)
end
# Pair-wise nan-covariance
function _nancor(x::AbstractVector{Tx}, y::AbstractVector{Ty}, corrected::Bool) where {Tx, Ty}
    # Parwise nan-means
    n = 0
    Σᵪ = ∅ᵪ = zero(Tx)
    Σᵧ = ∅ᵧ = zero(Ty)
    @inbounds @simd ivdep for i ∈ eachindex(x,y)
        xᵢ, yᵢ = x[i], y[i]
        notnan = (xᵢ==xᵢ) & (yᵢ==yᵢ)
        n += notnan
        Σᵪ += ifelse(notnan, xᵢ, ∅ᵪ)
        Σᵧ += ifelse(notnan, yᵢ, ∅ᵧ)
    end
    μᵪ = Σᵪ/n
    μᵧ = Σᵧ/n
    n == 0 && return (∅ᵪ+∅ᵧ)/0 # Return appropriate NaN if no data

    # Pairwise nan-variances
    σ²ᵪ = ∅²ᵪ = zero(Base.promote_op(*, typeof(μᵪ), typeof(μᵪ)))
    σ²ᵧ = ∅²ᵧ = zero(Base.promote_op(*, typeof(μᵧ), typeof(μᵧ)))
    @inbounds @simd ivdep for i ∈ eachindex(x,y)
        δᵪ = x[i] - μᵪ
        δᵧ = y[i] - μᵧ
        notnan = (δᵪ==δᵪ) & (δᵧ==δᵧ)
        σ²ᵪ += ifelse(notnan, δᵪ * δᵪ, ∅²ᵪ)
        σ²ᵧ += ifelse(notnan, δᵧ * δᵧ, ∅²ᵧ)
    end
    σᵪ = sqrt(σ²ᵪ / max(n-corrected, 0))
    σᵧ = sqrt(σ²ᵧ / max(n-corrected, 0))

    # Covariance and correlation
    σᵪᵧ = _nancov(x, y, corrected, μᵪ, μᵧ)
    ρᵪᵧ = σᵪᵧ / (σᵪ * σᵧ)

    return ρᵪᵧ
end


"""
```julia
nancor(X::AbstractMatrix; dims::Int=1)
```
Compute the (Pearson's product-moment) correlation matrix of the matrix `X`,
along dimension `dims`. As `Statistics.cor`, but ignoring `NaN`s.
"""
function nancor(X::AbstractMatrix; dims::Int=1, corrected::Bool=true)
    Tₒ = Base.promote_op(/, eltype(X), eltype(X))
    n = size(X, dims)
    m = size(X, mod(dims,2)+1)
    Ρ = similar(X, Tₒ, (m, m))
    # Diagonal must be unity
    @inbounds for i = 1:m
        Ρ[i,i] = one(Tₒ)
    end
    # Only two dimensions are possible, so handle each manually
    if dims == 1
        # Fill off-diagonals symmetrically
        @inbounds for i = 1:m
            for j = 1:i
                Ρ[i,j] = Ρ[j,i] = _nancor(view(X,:,i), view(X,:,j), corrected)
            end
        end
    elseif dims == 2
        @inbounds for i = 1:m
            for j = 1:i-1
                Ρ[i,j] = Ρ[j,i] = _nancor(view(X,i,:), view(X,j,:), corrected)
            end
        end
    else
        throw("Dimension not in range")
    end
    return Ρ
end
export nancor
