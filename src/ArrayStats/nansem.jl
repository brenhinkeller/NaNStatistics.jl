"""
```julia
nansem(A; dims=:, mean=nothing, corrected=true)
```
Compute the standard error of the mean of all non-`NaN` elements in `A`, optionally 
over dimensions specified by `dims`.

A precomputed `mean` may optionally be provided, which results in a somewhat faster
calculation. If `corrected` is `true`, then _Bessel's correction_ is applied, such
that the sum is divided by `n-1` rather than `n`.

As an alternative to `dims`, `nansem` also supports the `dim` keyword, which
behaves identically to `dims`, but also drops any singleton dimensions that have
been reduced over (as is the convention in some other languages).


## Examples
```julia
julia> using NaNStatistics

julia> A = [1 2; 3 4]
2×2 Matrix{Int64}:
 1  2
 3  4

julia> nansem(A, dims=1)
1×2 Matrix{Float64}:
 1.0  1.0

julia> nansem(A, dims=2)
2×1 Matrix{Float64}:
 0.5
 0.5
```
"""
nansem(A; dims=:, dim=:, mean=nothing, corrected=true) = __nansem(mean, corrected, A, dims, dim)
__nansem(mean, corrected, A, ::Colon, ::Colon) = _nansem(mean, corrected, A, :)
__nansem(mean, corrected, A, region, ::Colon) = _nansem(mean, corrected, A, region)
__nansem(mean, corrected, A, ::Colon, region) = reducedims(__nansem(mean, corrected, A, region, :), region)
export nansem

# If dims is an integer, wrap it in a tuple
_nansem(μ, corrected::Bool, A, dims::Int) = _nansem(μ, corrected, A, (dims,))

# If the mean isn't known, compute it
_nansem(::Nothing, corrected::Bool, A, dims::Tuple) = _nansem!(_nanmean(A, dims, NANMEAN_SIZE_THRESHOLD), corrected, A, dims)
# Reduce all the dims!
function _nansem(::Nothing, corrected::Bool, A, ::Colon)
  Tₒ = Base.promote_op(/, eltype(A), Int)
  n = 0
  Σ = ∅ = zero(Tₒ)
  @inbounds @simd ivdep for i ∈ eachindex(A)
      Aᵢ = A[i]
      notnan = Aᵢ==Aᵢ
      n += notnan
      Σ += ifelse(notnan, Aᵢ, ∅)
  end
  μ = Σ / n
  σ² = ∅² = zero(Base.promote_op(*,Tₒ,Tₒ))
  @inbounds @simd ivdep for i ∈ eachindex(A)
      Aᵢ = A[i]
      δ = Aᵢ - μ
      σ² += ifelse(Aᵢ==Aᵢ, δ * δ, ∅²)
  end
  return sqrt(σ² / max(n-corrected,0) / n)
end
function _nansem(::Nothing, corrected::Bool, A::AbstractArray{T}, ::Colon) where T<:Integer
    Tₒ = Base.promote_op(/, T, Int)
    n = length(A)
    Σ = zero(Tₒ)
    @inbounds @simd ivdep for i ∈ eachindex(A)
        Σ += A[i]
    end
    μ = Σ / n
    σ² = zero(Base.promote_op(*,Tₒ,Tₒ))
    @inbounds @simd ivdep for i ∈ eachindex(A)
        δ = A[i] - μ
        σ² += δ * δ
    end
    return sqrt(σ² / max(n-corrected,0) / n)
end


# If the mean is known, pass it on in the appropriate form
function _nansem(μ, corrected::Bool, A, dims::Tuple)
  μ_alloc = allocate_nanmean(A, dims)
  @assert length(μ) == length(μ_alloc) "Precomputed mean has $(length(μ)) elements, expected $(length(μ_alloc))"
  return _nansem!((μ_alloc .= μ), corrected, A, dims)
end
# Reduce all the dims!
function _nansem(μ::Number, corrected::Bool, A, ::Colon)
  n = 0
  σ² = ∅ = zero(typeof(μ))
  @inbounds @simd ivdep for i ∈ eachindex(A)
      Aᵢ = A[i]
      δ = Aᵢ - μ
      notnan = Aᵢ==Aᵢ
      n += notnan
      σ² += ifelse(notnan, δ * δ, ∅)
  end
  return sqrt(σ² / max(n-corrected, 0) / n)
end
function _nansem(μ::Number, corrected::Bool, A::AbstractArray{T}, ::Colon) where T<:Integer
    σ² = zero(typeof(μ))
    if μ==μ
        @inbounds @simd ivdep for i ∈ eachindex(A)
            δ = A[i] - μ
            σ² += δ * δ
        end
        n = length(A)
    else
        n = 0
    end
    return sqrt(σ² / max(n-corrected, 0) / n)
end


function staticdim_nansem_quote(static_dims::Vector{Int}, N::Int)
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
  # Reverse the axis lists so that we build up the loops from slow-axis to
  # fast-axis for column major arrays.
  reverse!(nonreduct_inds)
  reverse!(reduct_inds)
  
  # Secondly, build up our set of loops
  block = loops = Expr(:block)
  for n ∈ nonreduct_inds
    newblock = Expr(:block)
    push!(block.args, Expr(:for, :($(inds[n]) = axes(A,$n)), newblock))
    block = newblock
  end
  # Push more things here if you want them at the beginning of the reduction loop
  rblock = block
  push!(rblock.args, :(μ = $Bind))
  push!(rblock.args, :(n = 0))
  push!(rblock.args, :(σ² = ∅))
  # Build the reduction loop
  for n ∈ reduct_inds
    newblock = Expr(:block)
    push!(block.args, Expr(:for, :($(inds[n]) = axes(A,$n)), newblock))
    block = newblock
  end
  # Push more things here if you want them in the innermost loop
  push!(block.args, :(δ = $Aind - μ))
  push!(block.args, :(notnan = $Aind==$Aind))
  push!(block.args, :(n += notnan))
  push!(block.args, :(σ² += ifelse(notnan, δ * δ, ∅)))
  # Push more things here if you want them at the end of the reduction loop
  push!(rblock.args, :($Bind = sqrt(σ² * inv(max(n-corrected,0)) * inv(n))))

  # Put it all together
  quote
    ∅ = zero(eltype(B))
    Bᵥ = $Bᵥ
    @inbounds $loops
    return B
  end
end

function branches_nansem_quote(N::Int, M::Int, D)
  static_dims = Int[]
  for m ∈ 1:M
    param = D.parameters[m]
    if param <: _StaticInt
      new_dim = _dim(param)::Int
      @assert new_dim ∉ static_dims
      push!(static_dims, new_dim)
    else
      t = Expr(:tuple)
      for n ∈ static_dims
        push!(t.args, :(_StaticInt{$n}()))
      end
      q = Expr(:block, :(dimm = dims[$m]))
      qold = q
      ifsym = :if
      for n ∈ 1:N
        n ∈ static_dims && continue
        tc = copy(t)
        push!(tc.args, :(_StaticInt{$n}()))
        qnew = Expr(ifsym, :(dimm == $n), :(return _nansem!(B, corrected, A, $tc)))
        for r ∈ m+1:M
          push!(tc.args, :(dims[$r]))
        end
        push!(qold.args, qnew)
        qold = qnew
        ifsym = :elseif
      end
      # Else, if dimm ∉ 1:N, drop it from list and continue
      tc = copy(t)
      for r ∈ m+1:M
          push!(tc.args, :(dims[$r]))
      end
      push!(qold.args, Expr(:block, :(return _nansem!(B, corrected, A, $tc))))
      return q
    end
  end
  staticdim_nansem_quote(static_dims, N)
end

# Efficient @generated in-place var
@generated function _nansem!(B::AbstractArray{Tₒ,N}, corrected::Bool, A::AbstractArray{T,N}, dims::D) where {Tₒ,T,N,M,D<:Tuple{Vararg{_IntOrStaticInt,M}}}
  # N == M && return :(B[1] = _nansem(B[1], corrected, A, :); B)
  branches_nansem_quote(N, M, D)
end
