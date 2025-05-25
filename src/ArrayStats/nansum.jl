"""
```julia
nansum(A; dims)
```
Calculate the sum of an indexable collection `A`, ignoring `NaN`s, optionally
along dimensions specified by `dims`.

Also supports the `dim` keyword, which behaves identically to `dims`, but
also drops any singleton dimensions that have been reduced over (as is the
convention in some other languages).

## Examples
```julia
julia> using NaNStatistics

julia> A = [1 2; 3 4]
2×2 Matrix{Int64}:
 1  2
 3  4

julia> nansum(A, dims=1)
1×2 Matrix{Int64}:
 4  6

julia> nansum(A, dims=2)
2×1 Matrix{Int64}:
 3
 7
```
"""
nansum(A; dims=:, dim=:) = __nansum(A, dims, dim)
__nansum(A, ::Colon, ::Colon) = _nansum(A, :)
__nansum(A, region, ::Colon) = _nansum(A, region)
__nansum(A, ::Colon, region) = reducedims(__nansum(A, region, :), region)
export nansum


# Reduce one dim
_nansum(A, dims::Int) = _nansum(A, (dims,))

# Reduce some dims
function _nansum(A::AbstractArray{T,N}, dims::Tuple) where {T,N}
    if 1 in dims
        # Use the generated-function approach for reducing over the first
        # dimension since it's faster.
        sᵢ = size(A)
        sₒ = ntuple(Val{N}()) do d
            ifelse(d ∈ dims, 1, sᵢ[d])
        end
        B = similar(A, T, sₒ)
        _nansum!(B, A, dims)
    else
        # Base.sum() is faster at everything else
        sum(x -> ifelse(isnan(x), zero(x), x), A; dims)
    end
end
function _nansum(A::AbstractArray{T,N}, dims::Tuple) where {T<:Integer,N}
    sᵢ = size(A)
    sₒ = ntuple(Val{N}()) do d
        ifelse(d ∈ dims, 1, sᵢ[d])
    end
    Tₒ = Base.promote_op(+, T, Int)
    B = similar(A, Tₒ, sₒ)
    _nansum!(B, A, dims)
end

# Reduce all the dims!
function _nansum(A, ::Colon)
    Tₒ = eltype(A)
    Σ = ∅ = zero(Tₒ)
    @inbounds @simd ivdep for i ∈ eachindex(A)
        Aᵢ = A[i]
        notnan = Aᵢ==Aᵢ
        Σ += ifelse(notnan, Aᵢ, ∅)
    end
    return Σ
end
function _nansum(A::AbstractArray{T}, ::Colon) where T<:Integer
    Tₒ = Base.promote_op(+, T, Int)
    Σ = zero(Tₒ)
    @inbounds @simd ivdep for i ∈ eachindex(A)
        Σ += A[i]
    end
    return Σ
end


# Metaprogramming magic adapted from Chris Elrod example:
# Generate customized set of loops for a given ndims and a vector
# `static_dims` of dimensions to reduce over
function staticdim_nansum_quote(static_dims::Vector{Int}, N::Int)
    M = length(static_dims)
    # `static_dims` now contains every dim we're taking the sum over.
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
    firstn = first(nonreduct_inds)
    block = Expr(:block)
    loops = Expr(:for, :($(inds[firstn]) = indices((A,B),$firstn)), block)
    if length(nonreduct_inds) > 1
        for n ∈ @view(nonreduct_inds[2:end])
            newblock = Expr(:block)
            push!(block.args, Expr(:for, :($(inds[n]) = indices((A,B),$n)), newblock))
            block = newblock
        end
    end

    # Push more things here if you want them at the beginning of the reduction loop
    rblock = block
    push!(rblock.args, :(Σ = ∅))

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
    push!(block.args, :(Aᵢ = $Aind))
    push!(block.args, :(notnan = Aᵢ==Aᵢ))
    push!(block.args, :(Σ += ifelse(notnan, Aᵢ, ∅)))

    # Push more things here if you want them at the end of the reduction loop
    push!(rblock.args, :($Bind = Σ))

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
function branches_nansum_quote(N::Int, M::Int, D)
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
                qnew = Expr(ifsym, :(dimm == $n), :(return _nansum!(B, A, $tc)))
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
    staticdim_nansum_quote(static_dims, N)
end

# Efficient @generated in-place sum
@generated function _nansum!(B::AbstractArray{Tₒ,N}, A::AbstractArray{T,N}, dims::D) where {Tₒ,T,N,M,D<:Tuple{Vararg{IntOrStaticInt,M}}}
    N == M && return :(B[1] = _nansum(A, :); B)
    branches_nansum_quote(N, M, D)
end


## ---
