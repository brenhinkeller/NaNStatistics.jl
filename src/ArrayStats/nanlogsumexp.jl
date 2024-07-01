
"""
```julia
nanlogsumexp(A)
```
Return the logarithm of the sum of `exp.(A)` -- i.e., `log(sum(exp.(A)))`, but ignoring `NaN`s and avoiding numerical over/underflow.
As `nancumsum`, but operating on logarithms; as `nanlogsumexp`, but returning a array of cumulative sums, rather than a single value.
## Examples
```julia
```
"""
nanlogsumexp(A) = _nanlogsumexp(A, nanmaximum(A))
export nanlogsumexp

function _nanlogsumexp(A, c)
    Σ = ∅ = zero(Base.promote_op(exp, eltype(A)))
    @inbounds @simd ivdep for i ∈ eachindex(A)
        Aᵢ = exp(A[i] - c)
        Σ += ifelse(isnan(Aᵢ), ∅, Aᵢ)
    end
    return log(Σ) + c
end


"""
```julia
nanlogcumsumexp(A)
```
Return the logarithm of the cumulative sum of `exp.(A)` -- i.e., `log.(cumsum.(exp.(A)))`, but ignoring `NaN`s and avoiding numerical over/underflow.
As `nancumsum`, but operating on logarithms; as `nanlogsumexp`, but returning a array of cumulative sums, rather than a single value.
## Examples
```julia
```
"""
nanlogcumsumexp(A; reverse=false) = _nanlogcumsumexp(A, static(reverse))
export nanlogcumsumexp

function _nanlogcumsumexp(A, ::False)
    Tᵣ = Base.promote_op(exp, eltype(A))
    Σ = ∅ = zero(Tᵣ)
    lΣ = fill!(similar(A, Tᵣ), ∅)
    c = linit(eltype(A))
    i₀ = firstindex(A)
    @inbounds for i ∈ eachindex(A)
        isnan(c) && (c = A[i])
        if A[i] > c
            c = A[i]
            Σ = ∅ 
            @simd ivdep for j ∈ i₀:i
                Aᵢ = exp(A[j] - c)
                Σ += ifelse(isnan(Aᵢ), ∅, Aᵢ)
            end
            lΣ[i] = log(Σ) + c
        else
            Aᵢ = exp(A[i] - c)
            Σ += ifelse(isnan(Aᵢ), ∅, Aᵢ)
            lΣ[i] = log(Σ) + c
        end
    end
    return lΣ
end
function _nanlogcumsumexp(A, ::True)
    Tᵣ = Base.promote_op(exp, eltype(A))
    Σ = ∅ = zero(Tᵣ)
    lΣ = fill!(similar(A, Tᵣ), ∅)
    c = linit(eltype(A))
    iₗ = lastindex(A)
    @inbounds for i ∈ reverse(eachindex(A))
        isnan(c) && (c = A[i])
        if A[i] > c
            c = A[i]
            Σ = ∅ 
            @simd ivdep for j ∈ i:iₗ
                Aᵢ = exp(A[j] - c)
                Σ += ifelse(isnan(Aᵢ), ∅, Aᵢ)
            end
            lΣ[i] = log(Σ) + c
        else
            Aᵢ = exp(A[i] - c)
            Σ += ifelse(isnan(Aᵢ), ∅, Aᵢ)
            lΣ[i] = log(Σ) + c
        end
    end
    return lΣ
end