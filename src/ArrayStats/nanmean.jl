# Default to 1MB
NANMEAN_SIZE_THRESHOLD::Union{Int, Symbol} = 2^20

get_size_threshold(x::Integer) = x

"""
```julia
nanmean(A; dims, size_threshold)
```
Compute the mean of all non-`NaN` elements in `A`, optionally over dimensions
specified by `dims`. As `Statistics.mean`, but ignoring `NaN`s.

As an alternative to `dims`, `nanmean` also supports the `dim` keyword, which
behaves identically to `dims`, but also drops any singleton dimensions that have
been reduced over (as is the convention in some other languages).

`nanmean` has optimized implementations for big and small arrays for reducing
over slow dimensions. Which implementation is used can be tuned by setting the
`size_threshold` keyword argument to either an integer number of bytes, or the
symbols `:L1`/`:L2`/`:L3` to use a specific cache size. The Hwloc.jl package
must be loaded to support setting the threshold by cache symbol. The default
threshold is 1MB, but if Hwloc.jl is loaded then the default will be set to
`:L1` for the L1 cache size.

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
nanmean(A; dims=:, dim=:, size_threshold=NANMEAN_SIZE_THRESHOLD) = __nanmean(A, dims, dim, size_threshold)
__nanmean(A, ::Colon, ::Colon, st) = _nanmean(A, :, st)
__nanmean(A, region, ::Colon, st) = _nanmean(A, region, st)
__nanmean(A, ::Colon, region, st) = reducedims(__nanmean(A, region, :, st), region)
export nanmean


# Reduce one dim
_nanmean(A, dims::Int, st) = _nanmean(A, (dims,), st)

# Reduce some dims
function _nanmean(A::AbstractArray{T,N}, dims::Tuple, st) where {T,N}
    sᵢ = size(A)
    sₒ = ntuple(Val{N}()) do d
        ifelse(d ∈ dims, 1, sᵢ[d])
    end
    Tₒ = Base.promote_op(/, T, Int)
    B = similar(A, Tₒ, sₒ)

    if 1 in dims || sizeof(A) < get_size_threshold(st)
        # The generated-function approach is faster for small arrays and if
        # we're reducing over the first dimension.
        _nanmean!(B, A, dims, st)
    else
        # For reducing over the slow axes of large arrays we use the mapreduce
        # approach.
        _nanmean_mapreduce!(B, A, dims)
    end
end

# Reduce all the dims!
function _nanmean(A, ::Colon, _)
    Tₒ = Base.promote_op(/, eltype(A), Int)
    n = 0
    Σ = ∅ = zero(Tₒ)
    @inbounds @simd ivdep for i ∈ eachindex(A)
        Aᵢ = A[i]
        notnan = Aᵢ==Aᵢ
        n += notnan
        Σ += ifelse(notnan, A[i], ∅)
    end
    return Σ / n
end
function _nanmean(A::AbstractArray{T}, ::Colon, _) where T<:Integer
    Tₒ = Base.promote_op(/, T, Int)
    Σ = zero(Tₒ)
    @inbounds @simd ivdep for i ∈ eachindex(A)
        Σ += A[i]
    end
    return Σ / length(A)
end

# Tight loops for the mapreduce-style implementation. This is in a separate
# function to ensure type-stability.
function _nanmean_mapreduce_impl!(B, A, counts)
    @_mapreduce_impl(B, A,
                     if !isnan(x)
                         B[ir, IR] += x
                         counts[ir, IR] += 1
                     end)

    ∅ = zero(eltype(B))
    @inbounds for i in eachindex(B)
        if iszero(counts[i])
            B[i] = ∅
        else
            B[i] /= counts[i]
        end
    end
end

function _nanmean_mapreduce!(B, A, dims::Tuple)
    # Compute the length of the dimensions we're reducing over to pick
    # the smallest valid eltype for the counts array. The downside is type
    # instability and a few more allocations, but on large arrays this
    # can improve performance by ~25%.
    to_reduce_dims = ntuple(i -> size(B)[i] == 1 ? size(A)[i] : 1, ndims(B))
    max_reduction_len = prod(to_reduce_dims)
    counts_T = if max_reduction_len < typemax(UInt8)
        UInt8
    elseif max_reduction_len < typemax(UInt16)
        UInt16
    elseif max_reduction_len < typemax(UInt32)
        UInt32
    else
        UInt64
    end
    counts = similar(B, counts_T)
    fill!(counts, zero(counts_T))

    # We need to initialize B with zeros because we start off using it to sum
    # the array.
    fill!(B, zero(eltype(B)))

    _nanmean_mapreduce_impl!(B, A, counts)

    return B
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
    push!(rblock.args, :(n = 0))
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
    push!(block.args, :(n += notnan))
    push!(block.args, :(Σ += ifelse(notnan, Aᵢ, ∅)))

    # Push more things here if you want them at the end of the reduction loop
    push!(rblock.args, :($Bind = Σ * inv(n)))

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
function branches_nanmean_quote(N::Int, M::Int, D, st)
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
                qnew = Expr(ifsym, :(dimm == $n), :(return _nanmean!(B, A, $tc, st)))
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
@generated function _nanmean!(B::AbstractArray{Tₒ,N}, A::AbstractArray{T,N}, dims::D, st) where {Tₒ,T,N,M,D<:Tuple{Vararg{IntOrStaticInt,M}}}
    N == M && return :(B[1] = _nanmean(A, :, st); B)
    branches_nanmean_quote(N, M, D, st)
end

## ---
