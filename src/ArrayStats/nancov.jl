function _nancov(x::AbstractVector, y::AbstractVector, corrected::Bool, μᵪ::Number, μᵧ::Number)
    # Calculate covariance
    σᵪᵧ = ∅ = zero(promote_type(typeof(μᵪ), typeof(μᵧ), Int))
    n = 0
    @turbo for i ∈ indices((x,y))
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
    # Check lengths
    nᵪ = length(x)
    nᵧ = length(y)
    @assert nᵪ == nᵧ

    μᵪ = _nanmean(x,:)
    μᵧ = _nanmean(y,:)
    σᵪᵧ = _nancov(x, y, corrected, μᵪ, μᵧ)
    return σᵪᵧ
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
        # Precalculate means for each column
        μ = ntuple(m) do d
            nanmean(view(X,:,d))
        end
        # Fill covariance matrix symmetrically
        @inbounds for i = 1:m
            for j = 1:i
                σᵢⱼ = _nancov(view(X,:,i), view(X,:,j), corrected, μ[i], μ[j])
                Σ[i,j] = Σ[j,i] = σᵢⱼ
            end
        end
    elseif dims == 2
        # Precalculate means for each row
        μ = ntuple(m) do d
            nanmean(view(X,d,:))
        end
        # Fill covariance matrix symmetrically
        @inbounds for i = 1:m
            for j = 1:i
                σᵢⱼ = _nancov(view(X,i,:), view(X,j,:), corrected, μ[i], μ[j])
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
    # Check lengths
    nᵪ = length(x)
    nᵧ = length(y)
    @assert nᵪ == nᵧ

    μᵪ = nanmean(x)
    μᵧ = nanmean(y)
    σᵪ = nanstd(x, mean=μᵪ, corrected=corrected)
    σᵧ = nanstd(y, mean=μᵧ, corrected=corrected)
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
        # Precalculate means and standard deviations
        μ = ntuple(m) do d
            nanmean(view(X,:,d))
        end
        σ = ntuple(m) do d
            nanstd(view(X,:,d), mean=μ[d], corrected=corrected)
        end
        # Fill off-diagonals symmetrically
        @inbounds for i = 1:m
            for j = 1:i
                σᵢⱼ = _nancov(view(X,:,i), view(X,:,j), corrected, μ[i], μ[j])
                Ρ[i,j] = Ρ[j,i] = σᵢⱼ / (σ[i] * σ[j])
            end
        end
    elseif dims == 2
        # Precalculate means and standard deviations
        μ = ntuple(m) do d
            nanmean(view(X,d,:))
        end
        σ = ntuple(m) do d
            nanstd(view(X,d,:), mean=μ[d], corrected=corrected)
        end
        @inbounds for i = 1:m
            for j = 1:i-1
                σᵢⱼ = _nancov(view(X,i,:), view(X,j,:), corrected, μ[i], μ[j])
                Ρ[i,j] = Ρ[j,i] = σᵢⱼ / (σ[i] * σ[j])
            end
        end
    else
        throw("Dimension not in range")
    end
    return Ρ
end
export nancor
