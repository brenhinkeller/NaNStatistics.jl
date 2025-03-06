"""
```julia
nancumsum(A; dims)
```
Calculate the sum of an indexable collection `A`, ignoring `NaN`s, optionally
along dimensions specified by `dims`.

## Examples
```julia
julia> using NaNStatistics

julia> A = [1 2; 3 4; NaN NaN]
3×2 Matrix{Float64}:
   1.0    2.0
   3.0    4.0
 NaN    NaN

julia> nancumsum(A, dims=1)
3×2 Matrix{Float64}:
 1.0  2.0
 4.0  6.0
 4.0  6.0

julia> nancumsum(vec(A))
6-element Vector{Float64}:
  1.0
  4.0
  4.0
  6.0
 10.0
 10.0
```
"""
nancumsum(A; dims=:) = _nancumsum(A, dims)
export nancumsum

nancumsum!(A; dims=:) = _nancumsum!(A, dims) # In-place variant
export nancumsum!

# Reduce one dim
_nancumsum(A, dims::Int) = _nancumsum(A, (dims,))
_nancumsum!(A, dims::Int) = _nancumsum!(A, A, (dims,))

# Reduce some dims
function _nancumsum(A::AbstractArray{T,N}, dims::Tuple) where {T,N}
    Tₒ = Base.promote_op(+, T, Int)
    B = similar(A, Tₒ)
    _nancumsum!(B, A, dims)
end

# Reduce all the dims!
function _nancumsum(A::AbstractArray, ::Colon)
    Tₒ = Base.promote_op(+, eltype(A), Int)
    B = similar(A, Tₒ)
    _nancumsum!(B, A, :)
end
@inline function _nancumsum!(B::AbstractArray{T}, A::AbstractArray, ::Colon) where T
    Σ = ∅ = zero(T)
    @inbounds for i ∈ eachindex(A)
        Aᵢ = A[i]
        notnan = Aᵢ==Aᵢ
        Σ += ifelse(notnan, Aᵢ, ∅)
        B[i] = Σ
    end
    return B
end
@inline function _nancumsum!(A::AbstractArray{T}, ::Colon) where T
  Σ = ∅ = zero(T)
  @inbounds for i ∈ eachindex(A)
      Aᵢ = A[i]
      notnan = Aᵢ==Aᵢ
      Σ += ifelse(notnan, Aᵢ, ∅)
      A[i] = Σ
  end
  return A
end
@inline function _nancumsum!(B::AbstractArray{T}, A::AbstractArray{<:Integer}, ::Colon) where T
    Σ = zero(T)
    @inbounds for i ∈ eachindex(A)
        Σ += A[i]
        B[i] = Σ
    end
    return B
end
@inline function _nancumsum!(A::AbstractArray{T}, ::Colon) where T<:Integer
  Σ = zero(T)
  @inbounds for i ∈ eachindex(A)
      Σ += A[i]
      A[i] = Σ
  end
  return A
end



# Metaprogramming magic adapted from Chris Elrod example:
# Generate customized set of loops for a given ndims and a vector
# `static_dims` of dimensions to reduce over
function staticdim_nancumsum_quote(static_dims::Vector{Int}, N::Int)
  M = length(static_dims)
  # `static_dims` now contains every dim we're taking the sum over.
  Bᵥ = Expr(:call, :view, :B)
  reduct_inds = Int[]
  nonreduct_inds = Int[]
  # Firstly, build our expressions for indexing each array
  Aind = :(A[])
  Bind = :(B[])
  inds = Vector{Symbol}(undef, N)
  for n ∈ 1:N
    ind = Symbol(:i_,n)
    inds[n] = ind
    push!(Aind.args, ind)
    push!(Bind.args, ind)
    if n ∈ static_dims
      push!(reduct_inds, n)
    else
      push!(nonreduct_inds, n)
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
  push!(block.args, :(Σ += ifelse(notnan, Aᵢ, ∅)))
  push!(block.args, :($Bind = Σ))
  # Push more things here if you want them at the end of the reduction loop
  # push!(rblock.args, :(....))
  # Put it all together
  quote
    ∅ = zero(eltype(B))
    @inbounds $loops
    return B
  end
end

# Chris Elrod metaprogramming magic:
# Turn non-static integers in `dims` tuple into `StaticInt`s
# so we can construct `static_dims` vector within @generated code
function branches_nancumsum_quote(N::Int, M::Int, D)
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
        qnew = Expr(ifsym, :(dimm == $n), :(return _nancumsum!(B, A, $tc)))
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
  staticdim_nancumsum_quote(static_dims, N)
end

# Efficient @generated in-place sum
@generated function _nancumsum!(B::AbstractArray{Tₒ,N}, A::AbstractArray{T,N}, dims::D) where {Tₒ,T,N,M,D<:Tuple{Vararg{IntOrStaticInt,M}}}
  N == M && return _nancumsum!(B, A, :)
  branches_nancumsum_quote(N, M, D)
end

function staticdim_nancumsum!_quote(static_dims::Vector{Int}, N::Int)
  reduct_inds = Int[]
  nonreduct_inds = Int[]
  # Firstly, build our expressions for indexing each array
  Aind = :(A[])
  inds = Vector{Symbol}(undef, N)
  for n ∈ 1:N
    ind = Symbol(:i_,n)
    inds[n] = ind
    push!(Aind.args, ind)
    if n ∈ static_dims
      push!(reduct_inds, n)
    else
      push!(nonreduct_inds, n)
    end
  end
  firstn = first(nonreduct_inds)
  # Secondly, build up our set of loops
  block = Expr(:block)
  loops = Expr(:for, :($(inds[firstn]) = indices(A,$firstn)), block)
  if length(nonreduct_inds) > 1
    for n ∈ @view(nonreduct_inds[2:end])
      newblock = Expr(:block)
      push!(block.args, Expr(:for, :($(inds[n]) = indices(A,$n)), newblock))
      block = newblock
    end
  end
  rblock = block
  # Push more things here if you want them at the beginning of the reduction loop
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
  push!(block.args, :(Σ += ifelse(notnan, Aᵢ, ∅)))
  push!(block.args, :($Aind = Σ))
  # Push more things here if you want them at the end of the reduction loop
  # push!(rblock.args, :(....))
  # Put it all together
  quote
    ∅ = zero(eltype(A))
    @inbounds $loops
    return A
  end
end

function branches_nancumsum!_quote(N::Int, M::Int, D)
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
        qnew = Expr(ifsym, :(dimm == $n), :(return _nancumsum!(A, $tc)))
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
  staticdim_nancumsum!_quote(static_dims, N)
end

# Efficient @generated in-place sum
@generated function _nancumsum!(A::AbstractArray{T,N}, dims::D) where {T,N,M,D<:Tuple{Vararg{IntOrStaticInt,M}}}
  N == M && return _nancumsum!(A, :)
  branches_nancumsum!_quote(N, M, D)
end
