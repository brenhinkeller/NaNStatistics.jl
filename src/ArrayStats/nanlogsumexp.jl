"""
```julia
nanlogaddexp(a, b)
```
Return the logarithm of the sum of the exponentials of `a` and `b`, i.e., 
`log(exp(a) + exp(b))`, ignoring `NaN`s and avoiding numerical over/underflow.

## Examples
julia> nanlogaddexp(1,1)
1.6931471805599454

julia> nanlogaddexp(10000, 10000)
10000.69314718056
```
"""
function nanlogaddexp(a, b)
    if isnan(a)
        b
    elseif isnan(b)
        a
    elseif isinf(a)
        a > 0 ? a : b
    elseif isinf(b)
        b > 0 ? b : a
    elseif a < b
        log1p(exp(a - b)) + b
    else
        log1p(exp(b - a)) + a
    end
end
export nanlogaddexp

"""
```julia
nanlogsumexp(A)
```
Return the logarithm of the sum of `exp.(A)` -- i.e., `log(sum(exp.(A)))`, but ignoring `NaN`s and avoiding numerical over/underflow.
## Examples
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
    lΣ = similar(A, Tᵣ)
    @assert eachindex(A) == eachindex(lΣ) # Should always be true, but just to be sure
    if !isempty(A)
        i = firstindex(A)
        lΣ[i] = lΣ₋ = A[i]
        @inbounds for i ∈ Iterators.drop(eachindex(A),1)
            lΣ[i] = lΣ₋ = nanlogaddexp(lΣ₋, A[i])
        end
    end
    return lΣ
end
function _nanlogcumsumexp(A, ::_True)
    Tᵣ = Base.promote_op(exp, eltype(A))
    lΣ = similar(A, Tᵣ)
    @assert eachindex(A) == eachindex(lΣ) # Should always be true, but just to be sure
    if !isempty(A)
        i = lastindex(A)
        lΣ[i] = lΣ₋ = A[i]
        @inbounds for i ∈ Iterators.drop(reverse(eachindex(A)),1)
            lΣ[i] = lΣ₋ = nanlogaddexp(lΣ₋, A[i])
        end
    end
    return lΣ
end