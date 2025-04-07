"""
```julia
nanskewness(A; dims=:, mean=nothing)
```
Compute the skewness of all non-`NaN` elements in `A`, optionally over dimensions specified by `dims`.
As `StatsBase.skewness`, but ignoring `NaN`s.

A precomputed `mean` may optionally be provided, which results in a somewhat faster
calculation.

As an alternative to `dims`, `nanskewness` also supports the `dim` keyword, which
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

julia> nanskewness(A, dims=1)
1×3 Matrix{Float64}:
 0.0  0.0  0.0


julia> nanskewness(A, dims=2)
3-element Vector{Float64}:
 0.0
 0.0
 0.0
```
"""
nanskewness(A; dims=:, dim=:, mean=nothing, corrected=false) = __nanskewness(mean, corrected, A, dims, dim)
__nanskewness(mean, corrected, A, ::Colon, ::Colon) = _nanskewness(mean, corrected, A, :)
__nanskewness(mean, corrected, A, region, ::Colon) = _nanskewness(mean, corrected, A, region)
__nanskewness(mean, corrected, A, ::Colon, region) = reducedims(__nanskewness(mean, corrected, A, region, :), region)
export nanskewness

# If dims is an integer, wrap it in a tuple
_nanskewness(μ, corrected::Bool, A, dims::Int) = _nanskewness(μ, corrected, A, (dims,))

# If the mean isn't known, compute it
_nanskewness(::Nothing, corrected::Bool, A, dims::Tuple) = _nanskewness!(_nanmean(A, dims), corrected, A, dims)
# Reduce all the dims!
function _nanskewness(::Nothing, corrected::Bool, A, ::Colon)
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
    μ₃ = ∅³ = zero(Base.promote_op(*, typeof(μ₂), typeof(μ)))
    @inbounds @simd ivdep for i ∈ eachindex(A)
        δ = A[i] - μ
        notnan = δ==δ
        δ² = δ * δ
        μ₂ += ifelse(notnan, δ², ∅²)
        μ₃ += ifelse(notnan, δ² * δ, ∅³)
    end
    σ = sqrt(μ₂ / max(n-corrected,0))
    return (μ₃/n)/σ^3
end
function _nanskewness(::Nothing, corrected::Bool, A::AbstractArray{T}, ::Colon) where T<:Integer
    n = length(A)
    Σ = zero(T)
    @inbounds @simd ivdep for i ∈ eachindex(A)
        Σ += A[i]
    end
    μ = Σ / n
    μ₂ = zero(Base.promote_op(*, typeof(μ), typeof(μ)))
    μ₃ = zero(Base.promote_op(*, typeof(μ₂), typeof(μ)))
    @inbounds @simd ivdep for i ∈ eachindex(A)
        δ = A[i] - μ
        δ² = δ * δ
        μ₂ += δ²
        μ₃ += δ² * δ
    end
    σ = sqrt(μ₂ / max(n-corrected,0))
    return (μ₃/n)/σ^3
end


# If the mean is known, pass it on in the appropriate form
_nanskewness(μ, corrected::Bool, A, dims::Tuple) = _nanskewness!(collect(μ), corrected, A, dims)
_nanskewness(μ::AbstractArray, corrected::Bool, A, dims::Tuple) = _nanskewness!(copy(μ), corrected, A, dims)
_nanskewness(μ::Number, corrected::Bool, A, dims::Tuple) = _nanskewness!([μ], corrected, A, dims)
# Reduce all the dims!
function _nanskewness(μ::Number, corrected::Bool, A, ::Colon)
    n = 0
    μ₃ = μ₂ = ∅ₒ = zero(typeof(μ))
    @inbounds @simd ivdep for i ∈ eachindex(A)
        δ = A[i] - μ
        notnan = δ==δ
        n += notnan
        δ² = δ * δ
        μ₂ += ifelse(notnan, δ², ∅ₒ)
        μ₃ += ifelse(notnan, δ² * δ, ∅ₒ)
    end
    σ = sqrt(μ₂ / max(n-corrected,0))
    return (μ₃/n)/σ^3
end
function _nanskewness(μ::Number, corrected::Bool, A::AbstractArray{T}, ::Colon) where T<:Integer
    μ₃ = μ₂ = zero(typeof(μ))
    if μ==μ
        @inbounds @simd ivdep for i ∈ eachindex(A)
            δ = A[i] - μ
            δ² = δ * δ
            μ₂ += δ²
            μ₃ += δ² * δ
        end
        n = length(A)
    else
        n = 0
    end
    σ = sqrt(μ₂ / max(n-corrected,0))
    return (μ₃/n)/σ^3
end

# Metaprogramming magic adapted from Chris Elrod example:
# Generate customized set of loops for a given ndims and a vector
# `static_dims` of dimensions to reduce over
function staticdim_nanskewness_quote(static_dims::Vector{Int}, N::Int)
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
  push!(rblock.args, :(μ₃ = μ₂ = ∅))
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
  push!(block.args, :(μ₃ += ifelse(notnan, δ² * δ, ∅)))

  # Push more things here if you want them at the end of the reduction loop
  push!(rblock.args, :(σ = sqrt(μ₂ * inv(max(n-corrected,0))) ))
  push!(rblock.args, :($Bind = (μ₃/n)/σ^3 ))

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
function branches_nanskewness_quote(N::Int, M::Int, D)
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
        qnew = Expr(ifsym, :(dimm == $n), :(return _nanskewness!(B, corrected, A, $tc)))
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
  staticdim_nanskewness_quote(static_dims, N)
end

# Efficient @generated in-place var
@generated function _nanskewness!(B::AbstractArray{Tₒ,N}, corrected::Bool, A::AbstractArray{T,N}, dims::D) where {Tₒ,T,N,M,D<:Tuple{Vararg{IntOrStaticInt,M}}}
  N == M && return :(B[1] = _nanskewness(B[1], corrected, A, :); B)
  branches_nanskewness_quote(N, M, D)
end
