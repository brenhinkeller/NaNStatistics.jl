function _nancov(x::AbstractVector, y::AbstractVector, corrected::Bool, μᵪ::Number, μᵧ::Number)
    # Calculate covariance
    σᵪᵧ = ∅ = zero(promote_type(typeof(μᵪ), typeof(μᵧ), Int))
    n = 0
    @turbo check_empty=true for i ∈ eachindex(x,y)
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
    @turbo check_empty=true for i ∈ eachindex(x,y)
        xᵢ, yᵢ = x[i], y[i]
        notnan = (xᵢ==xᵢ) & (yᵢ==yᵢ)
        n += notnan
        Σᵪ += ifelse(notnan, xᵢ, ∅ᵪ)
        Σᵧ += ifelse(notnan, yᵢ, ∅ᵧ)
    end
    μᵪ = Σᵪ/n
    μᵧ = Σᵧ/n

    return _nancov(x, y, corrected, μᵪ, μᵧ)
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
    Tₒ = Base.promote_op(/, eltype(X), Int)
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
function _nancor(x::StridedVector{T}, y::StridedVector{T}, corrected::Bool) where T<:PrimitiveNumber
    # Parwise nan-means
    n = 0
    Σᵪ = ∅ᵪ = zero(T)
    Σᵧ = ∅ᵧ = zero(T)
    @turbo check_empty=true for i ∈ eachindex(x,y)
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
    σ²ᵪ = ∅ᵪ = zero(typeof(μᵪ))
    σ²ᵧ = ∅ᵧ = zero(typeof(μᵧ))
    @turbo check_empty=true for i ∈ eachindex(x,y)
        δᵪ = x[i] - μᵪ
        δᵧ = y[i] - μᵧ
        notnan = (δᵪ==δᵪ) & (δᵧ==δᵧ)
        σ²ᵪ += ifelse(notnan, δᵪ * δᵪ, ∅ᵪ)
        σ²ᵧ += ifelse(notnan, δᵧ * δᵧ, ∅ᵧ)
    end
    σᵪ = sqrt(σ²ᵪ / max(n-corrected, 0))
    σᵧ = sqrt(σ²ᵧ / max(n-corrected, 0))

    # Covariance and correlation
    σᵪᵧ = _nancov(x, y, corrected, μᵪ, μᵧ)
    ρᵪᵧ = σᵪᵧ / (σᵪ * σᵧ)

    return ρᵪᵧ
end
function _nancor(x::AbstractVector, y::AbstractVector, corrected::Bool)
    # Parwise nan-means
    n = 0
    Σᵪ = ∅ᵪ = zero(eltype(x))
    Σᵧ = ∅ᵧ = zero(eltype(y))
    @inbounds for i ∈ eachindex(x,y)
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
    σ²ᵪ = ∅ᵪ = zero(typeof(μᵪ))
    σ²ᵧ = ∅ᵧ = zero(typeof(μᵧ))
    @inbounds for i ∈ eachindex(x,y)
        δᵪ = x[i] - μᵪ
        δᵧ = y[i] - μᵧ
        notnan = (δᵪ==δᵪ) & (δᵧ==δᵧ)
        σ²ᵪ += ifelse(notnan, δᵪ * δᵪ, ∅ᵪ)
        σ²ᵧ += ifelse(notnan, δᵧ * δᵧ, ∅ᵧ)
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
    Tₒ = Base.promote_op(/, eltype(X), Int)
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
