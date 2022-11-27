"""
```julia
nanmean(A; dims)
```
Compute the mean of all non-`NaN` elements in `A`, optionally over dimensions
specified by `dims`. As `Statistics.mean`, but ignoring `NaN`s.

As an alternative to `dims`, `nanmean` also supports the `dim` keyword, which
behaves identically to `dims`, but also drops any singleton dimensions that have
been reduced over (as is the convention in some other languages).

## Examples
```julia
julia> using NaNStatistics

julia> A = [1 2; 3 4]
2×2 Matrix{Int64}:
 1  2
 3  4

julia> nanmean(A, dims=1)
1×2 Matrix{Float64}:
 2.0  3.0

julia> nanmean(A, dims=2)
2×1 Matrix{Float64}:
 1.5
 3.5
```
"""
nanmean(A; dims=:, dim=:) = __nanmean(A, dims, dim)
__nanmean(A, ::Colon, ::Colon) = _nanmean(A, :)
__nanmean(A, region, ::Colon) = _nanmean(A, region)
__nanmean(A, ::Colon, region) = reducedims(_nanmean(A, region), region)
export nanmean


# Reduce one dim
_nanmean(A, dims::Int) = _nanmean(A, (dims,))

# Reduce some dims
function _nanmean(A::AbstractArray{T,N}, dims::Tuple) where {T,N}
    sᵢ = size(A)
    sₒ = ntuple(Val(N)) do d
        ifelse(d ∈ dims, 1, sᵢ[d])
    end
    Tₒ = Base.promote_op(/, T, Int)
    B = similar(A, Tₒ, sₒ)
    _nanmean!(B, A, dims)
end

# Reduce all the dims!
function _nanmean(A::AbstractArray, ::Colon)
    Tₒ = Base.promote_op(/, eltype(A), Int)
    n = 0
    Σ = ∅ = zero(Tₒ)
    @turbo for i ∈ eachindex(A)
        Aᵢ = A[i]
        notnan = Aᵢ==Aᵢ
        n += notnan
        Σ += ifelse(notnan, A[i], ∅)
    end
    return Σ / n
end
function _nanmean(A::AbstractArray{<:Integer}, ::Colon)
    Tₒ = Base.promote_op(/, eltype(A), Int)
    Σ = zero(Tₒ)
    @turbo for i ∈ eachindex(A)
        Σ += A[i]
    end
    return Σ / length(A)
end
# Fallback method for non-arrays
function _nanmean(A, ::Colon)
    Tₒ = Base.promote_op(/, eltype(A), Int)
    n = 0
    Σ = ∅ = zero(Tₒ)
    @inbounds for i ∈ eachindex(A)
        Aᵢ = A[i]
        notnan = Aᵢ==Aᵢ
        n += notnan
        Σ += ifelse(notnan, A[i], ∅)
    end
    return Σ / n
end


# Metaprogramming magic adapted from Chris Elrod example:
# Generate customized set of loops for a given ndims and a vector
# `static_dims` of dimensions to reduce over
function staticdim_nanmean_quote(static_dims::Vector{Int}, N::Int)
  M = length(static_dims)
  # `static_dims` now contains every dim we're taking the mean over.
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
  push!(rblock.args, :(n = 0))
  push!(rblock.args, :(Σ = ∅))
  # Build the reduction loop
  for n ∈ reduct_inds
    newblock = Expr(:block)
    push!(block.args, Expr(:for, :($(inds[n]) = axes(A,$n)), newblock))
    block = newblock
  end
  # Push more things here if you want them in the innermost loop
  push!(block.args, :(Aᵢ = $Aind))
  push!(block.args, :(notnan = Aᵢ==Aᵢ))
  push!(block.args, :(n += notnan))
  push!(block.args, :(Σ += ifelse(notnan, Aᵢ, ∅)))
  # Push more things here if you want them at the end of the reduction loop
  push!(rblock.args, :($Bind = Σ * inv(n)))
  # Put it all together
  quote
    ∅ = zero(eltype(B))
    Bᵥ = $Bᵥ
    @turbo $loops
    return B
  end
end

# Chris Elrod metaprogramming magic:
# Turn non-static integers in `dims` tuple into `StaticInt`s
# so we can construct `static_dims` vector within @generated code
function branches_nanmean_quote(N::Int, M::Int, D)
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
        qnew = Expr(ifsym, :(dimm == $n), :(return _nanmean!(B, A, $tc)))
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
  staticdim_nanmean_quote(static_dims, N)
end

# Efficient @generated in-place mean
@generated function _nanmean!(B::AbstractArray{Tₒ,N}, A::AbstractArray{T,N}, dims::D) where {Tₒ,T,N,M,D<:Tuple{Vararg{IntOrStaticInt,M}}}
  N == M && return :(B[1] = _nanmean(A, :); B)
  branches_nanmean_quote(N, M, D)
end


## ---
