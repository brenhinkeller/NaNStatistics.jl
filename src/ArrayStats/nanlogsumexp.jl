
"""
```julia
nanlogsumexp(A)
```
Return the logarithm of the sum of `exp.(A)` -- i.e., `log(sum(exp.(A)))`, but ignoring `NaN`s and avoiding numerical over/underflow.
## Examples
```julia
julia> A = -10:-1
-10:-1

julia> nanlogsumexp(A)
-0.5413702555732887
```
"""
nanlogsumexp(A) = _nanlogsumexp(A, nanmaximum(A))
export nanlogsumexp

function _nanlogsumexp(A, c)
    isfinite(c) || return c
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
As `nanlogsumexp`, but returning a array of cumulative sums, rather than a single value.
## Examples
```julia
julia> A = -10:-1
-10:-1

julia> nanlogcumsumexp(A)
10-element Vector{Float64}:
 -10.0
 -8.686738312481777
 -7.592394035555619
 -6.559810301438804
 -5.548085604062407
 -4.543806683981877
 -3.5422371525957574
 -2.541660373520995
 -1.5414482720326212
 -0.5413702555732887
```
"""
nanlogcumsumexp(A; reverse=false) = _nanlogcumsumexp(A, _static(reverse))
export nanlogcumsumexp

function _nanlogcumsumexp(A, ::_False)
    Tᵣ = Base.promote_op(exp, eltype(A))
    Σ = zero(Tᵣ)
    lΣ = similar(A, Tᵣ)
    if !isempty(A)
        c = first(A)
        @inbounds for i ∈ eachindex(A)
            if isnan(c)
                isnan(A[i]) || (Σ += one(Tᵣ))
                lΣ[i] = c = A[i]
            else
                isnan(A[i]) || (Σ += exp(A[i] - c))
                if c < A[i] && isfinite(A[i])
                    Σ *= exp(c - A[i])
                    c = A[i]
                end
                lΣ[i] = log(Σ) + c
            end
        end
    end
    return lΣ
end
function _nanlogcumsumexp(A, ::_True)
    Tᵣ = Base.promote_op(exp, eltype(A))
    Σ = zero(Tᵣ)
    lΣ = similar(A, Tᵣ)
    if !isempty(A)
        c = last(A)
        @inbounds for i ∈ reverse(eachindex(A))
            if isnan(c)
                isnan(A[i]) || (Σ += one(Tᵣ))
                lΣ[i] = c = A[i]
            else
                isnan(A[i]) || (Σ += exp(A[i] - c))
                if c < A[i] && isfinite(A[i])
                    Σ *= exp(c - A[i])
                    c = A[i]
                end
                lΣ[i] = log(Σ) + c
            end
        end
    end
    return lΣ
end