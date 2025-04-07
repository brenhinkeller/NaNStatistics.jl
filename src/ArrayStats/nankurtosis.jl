"""
```julia
nankurtosis(A; dims=:, mean=nothing)
```
Compute the kurtosis of all non-`NaN` elements in `A`, optionally over dimensions specified by `dims`.
As `StatsBase.kurtosis`, but ignoring `NaN`s.

A precomputed `mean` may optionally be provided, which results in a somewhat faster
calculation. If `corrected` is `true`, then _Bessel's correction_ is applied, such
that the sum is divided by `n-1` rather than `n`.

As an alternative to `dims`, `nankurtosis` also supports the `dim` keyword, which
behaves identically to `dims`, but also drops any singleton dimensions that have
been reduced over (as is the convention in some other languages).


## Examples
```julia
julia> using NaNStatistics

julia> A = [1 2 3 ; 4 5 6; 7 8 9]
3×3 Matrix{Int64}:
 1  2  3
 4  5  6
 7  8  9

julia> nankurtosis(A, dims=1)
1×3 Matrix{Float64}:
 -1.5  -1.5  -1.5

julia> nankurtosis(A, dims=2)
3-element Vector{Float64}:
 -1.5
 -1.5
 -1.5
```
"""
nankurtosis(A; dims=:, dim=:, mean=nothing, corrected=false) = __nankurtosis(mean, corrected, A, dims, dim)
__nankurtosis(mean, corrected, A, ::Colon, ::Colon) = _nankurtosis(mean, corrected, A, :)
__nankurtosis(mean, corrected, A, region, ::Colon) = _nankurtosis(mean, corrected, A, region)
__nankurtosis(mean, corrected, A, ::Colon, region) = reducedims(__nankurtosis(mean, corrected, A, region, :), region)
export nankurtosis

# If dims is an integer, wrap it in a tuple
_nankurtosis(μ, corrected::Bool, A, dims::Int) = _nankurtosis(μ, corrected, A, (dims,))

# If the mean isn't known, compute it
_nankurtosis(::Nothing, corrected::Bool, A, dims::Tuple) = _nankurtosis!(_nanmean(A, dims), corrected, A, dims)
# Reduce all the dims!
function _nankurtosis(::Nothing, corrected::Bool, A, ::Colon)
    T = eltype(A)
    n = 0
    Σ = ∅ = zero(T)
    @inbounds @simd ivdep for i ∈ eachindex(A)
        Aᵢ = A[i]
        notnan = Aᵢ==Aᵢ
        n += notnan
        Σ += ifelse(notnan, Aᵢ, ∅)
    end
    μ = Σ / n
    μ₂ = ∅² = zero(Base.promote_op(*, typeof(μ), typeof(μ)))
    μ₄ = ∅⁴ = zero(Base.promote_op(*, typeof(μ₂), typeof(μ₂)))
    @inbounds @simd ivdep for i ∈ eachindex(A)
        δ = A[i] - μ
        notnan = δ==δ
        δ² = δ * δ
        μ₂ += ifelse(notnan, δ², ∅²)
        μ₄ += ifelse(notnan, δ² * δ², ∅⁴)
    end
    σ = sqrt(μ₂ / max(n-corrected,0))
    return (μ₄/n)/σ^4 - 3
end
function _nankurtosis(::Nothing, corrected::Bool, A::AbstractArray{T}, ::Colon) where T<:Integer
    n = length(A)
    Σ = zero(T)
    @inbounds @simd ivdep for i ∈ eachindex(A)
        Σ += A[i]
    end
    μ = Σ / n
    μ₂ = zero(Base.promote_op(*, typeof(μ), typeof(μ)))
    μ₄ = zero(Base.promote_op(*, typeof(μ₂), typeof(μ₂)))
    @inbounds @simd ivdep for i ∈ eachindex(A)
        δ = A[i] - μ
        δ² = δ * δ
        μ₂ += δ²
        μ₄ += δ² * δ²
    end
    σ = sqrt(μ₂ / max(n-corrected,0))
    return (μ₄/n)/σ^4 - 3
end


# If the mean is known, pass it on in the appropriate form
_nankurtosis(μ, corrected::Bool, A, dims::Tuple) = _nankurtosis!(collect(μ), corrected, A, dims)
_nankurtosis(μ::AbstractArray, corrected::Bool, A, dims::Tuple) = _nankurtosis!(copy(μ), corrected, A, dims)
_nankurtosis(μ::Number, corrected::Bool, A, dims::Tuple) = _nankurtosis!([μ], corrected, A, dims)
# Reduce all the dims!
function _nankurtosis(μ::Number, corrected::Bool, A, ::Colon)
    n = 0
    μ₂ = ∅² = zero(Base.promote_op(*, typeof(μ), typeof(μ)))
    μ₄ = ∅⁴ = zero(Base.promote_op(*, typeof(μ₂), typeof(μ₂)))
    @inbounds @simd ivdep for i ∈ eachindex(A)
        δ = A[i] - μ
        notnan = δ==δ
        n += notnan
        δ² = δ * δ
        μ₂ += ifelse(notnan, δ², ∅²)
        μ₄ += ifelse(notnan, δ² * δ², ∅⁴)
    end
    σ = sqrt(μ₂ / max(n-corrected,0))
    return (μ₄/n)/σ^4 - 3
end
function _nankurtosis(μ::Number, corrected::Bool, A::AbstractArray{T}, ::Colon) where T<:Integer
    μ₂ = zero(Base.promote_op(*, typeof(μ), typeof(μ)))
    μ₄ = zero(Base.promote_op(*, typeof(μ₂), typeof(μ₂)))
    if μ==μ
        @inbounds @simd ivdep for i ∈ eachindex(A)
            δ = A[i] - μ
            δ² = δ * δ
            μ₂ += δ²
            μ₄ += δ² * δ²
        end
        n = length(A)
    else
        n = 0
    end
    σ = sqrt(μ₂ / max(n-corrected,0))
    return (μ₄/n)/σ^4 - 3
end

# Metaprogramming magic adapted from Chris Elrod example:
# Generate customized set of loops for a given ndims and a vector
# `static_dims` of dimensions to reduce over
function staticdim_nankurtosis_quote(static_dims::Vector{Int}, N::Int)
  M = length(static_dims)
  # `static_dims` now contains every dim we're taking the var over.
  Bᵥ = Expr(:call, :view, :B)
  reduct_inds = Int[]
  nonreduct_inds = Int[]
  # Firstly, build our expressions for indexing each array
  Aind = :(A[])
  Bind = :(Bᵥ[])
  inds = Vector{Symbol}(undef, N)
  for n ∈ 1:N
    ind = Symbol(:i_,n)
    inds[n] = ind
    push!(Aind.args, ind)
    if n ∈ static_dims
      push!(reduct_inds, n)
      push!(Bᵥ.args, :(firstindex(B,$n)))
    else
      push!(nonreduct_inds, n)
      push!(Bᵥ.args, :)
      push!(Bind.args, ind)
    end
  end
  firstn = first(nonreduct_inds)
  # Secondly, build up our set of loops
  block = Expr(:block)
  loops = Expr(:for, :($(inds[firstn]) = indices((A,B),$firstn)), block)
  if length(nonreduct_inds) > 1
    for n ∈ @view(nonreduct_inds[2:end])
      newblock = Expr(:block)
      push!(block.args, Expr(:for, :($(inds[n]) = indices((A,B),$n)), newblock))
      block = newblock
    end
  end
  rblock = block
  # Push more things here if you want them at the beginning of the reduction loop
  push!(rblock.args, :(μ = $Bind))
  push!(rblock.args, :(n = 0))
  push!(rblock.args, :(μ₄ = μ₂ = ∅))
  # Build the reduction loop
  for n ∈ reduct_inds
    newblock = Expr(:block)
    if n==last(reduct_inds)
      push!(block.args, Expr(:macrocall, Symbol("@simd"), :ivdep, Expr(:for, :($(inds[n]) = axes(A,$n)), newblock)))
    else
      push!(block.args, Expr(:for, :($(inds[n]) = axes(A,$n)), newblock))
    end
    block = newblock
  end
  # Push more things here if you want them in the innermost loop
  push!(block.args, :(δ = $Aind - μ))
  push!(block.args, :(notnan = δ==δ))
  push!(block.args, :(n += notnan))
  push!(block.args, :(δ² = δ * δ))
  push!(block.args, :(μ₂ += ifelse(notnan, δ², ∅)))
  push!(block.args, :(μ₄ += ifelse(notnan, δ² * δ², ∅)))

  # Push more things here if you want them at the end of the reduction loop
  push!(rblock.args, :(σ = sqrt(μ₂ * inv(max(n-corrected,0))) ))
  push!(rblock.args, :($Bind = (μ₄/n)/σ^4 - 3))

  # Put it all together
  quote
    ∅ = zero(eltype(B))
    Bᵥ = $Bᵥ
    @inbounds $loops
    return B
  end
end

# Turn non-static integers in `dims` tuple into `StaticInt`s
# so we can construct `static_dims` vector within @generated code
function branches_nankurtosis_quote(N::Int, M::Int, D)
  static_dims = Int[]
  for m ∈ 1:M
    param = D.parameters[m]
    if param <: StaticInt
      new_dim = _dim(param)::Int
      @assert new_dim ∉ static_dims
      push!(static_dims, new_dim)
    else
      t = Expr(:tuple)
      for n ∈ static_dims
        push!(t.args, :(StaticInt{$n}()))
      end
      q = Expr(:block, :(dimm = dims[$m]))
      qold = q
      ifsym = :if
      for n ∈ 1:N
        n ∈ static_dims && continue
        tc = copy(t)
        push!(tc.args, :(StaticInt{$n}()))
        qnew = Expr(ifsym, :(dimm == $n), :(return _nankurtosis!(B, corrected, A, $tc)))
        for r ∈ m+1:M
          push!(tc.args, :(dims[$r]))
        end
        push!(qold.args, qnew)
        qold = qnew
        ifsym = :elseif
      end
      push!(qold.args, Expr(:block, :(throw("Dimension `$dimm` not found."))))
      return q
    end
  end
  staticdim_nankurtosis_quote(static_dims, N)
end

# Efficient @generated in-place var
@generated function _nankurtosis!(B::AbstractArray{Tₒ,N}, corrected::Bool, A::AbstractArray{T,N}, dims::D) where {Tₒ,T,N,M,D<:Tuple{Vararg{IntOrStaticInt,M}}}
  N == M && return :(B[1] = _nankurtosis(B[1], corrected, A, :); B)
  branches_nankurtosis_quote(N, M, D)
end
