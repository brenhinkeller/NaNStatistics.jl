# Check for sortedness, assuming no NaNs
@inline function issortedrange(A, iₗ, iᵤ)
    @inbounds for i = iₗ+1:iᵤ
        if A[i-1] > A[i]
            return false
        end
    end
    return true
end

# Check for anti-sortedness, assuming no NaNs
@inline function isantisortedrange(A, iₗ, iᵤ)
    @inbounds for i = iₗ+1:iᵤ
        if A[i-1] < A[i]
            return false
        end
    end
    return true
end

# Reverse an array, faster than Base.reverse!
@inline function vreverse!(A, iₗ, iᵤ)
    N = (iᵤ - iₗ) + 1
    n = (N ÷ 2) - 1
    if N < 32
        @inbounds for i ∈ 0:n
            𝔦ₗ, 𝔦ᵤ = iₗ+i, iᵤ-i
            A[𝔦ₗ], A[𝔦ᵤ] = A[𝔦ᵤ], A[𝔦ₗ]
        end
    else
        @inbounds @simd ivdep for i ∈ 0:n
            𝔦ₗ = iₗ+i
            𝔦ᵤ = iᵤ-i
            l = A[𝔦ₗ]
            u = A[𝔦ᵤ]
            A[𝔦ₗ] = u
            A[𝔦ᵤ] = l
        end
    end
    return A
end

# Move all NaNs to the end of the array A
function sortnans!(A::AbstractArray, iₗ::Int=firstindex(A), iᵤ::Int=lastindex(A))
    # Return early if range is empty
    iₗ >= iᵤ && return A, iₗ, iᵤ

    # Count up NaNs
    Nₙₐₙ = 0
    @inbounds @simd ivdep for i = iₗ:iᵤ
        Nₙₐₙ += isnan(A[i])
    end
    # If none, return early
    Nₙₐₙ == 0 && return A, iₗ, iᵤ

    # Otherwise, swap all NaNs
    i = iₗ
    j = iᵤ
    N = iᵤ - iₗ
    @inbounds for n = 0:N-Nₙₐₙ
        i = iₗ + n
        if A[i] != A[i]
            while A[j] != A[j]
                j -= 1
            end
            j <= i && break
            A[i], A[j] = A[j], A[i]
            j -= 1
        end
    end
    return A, iₗ, iᵤ - Nₙₐₙ
end
function argsortnans!(I::AbstractArray, A::AbstractArray, iₗ::Int=firstindex(A), iᵤ::Int=lastindex(A))
    # Return early if range is empty
    iₗ >= iᵤ && return I, A, iₗ, iᵤ
    
    # Count up NaNs
    Nₙₐₙ = 0
    @inbounds @simd ivdep for i = iₗ:iᵤ
        Nₙₐₙ += isnan(A[i])
    end
    # If none, return early
    Nₙₐₙ == 0 && return I, A, iₗ, iᵤ

    # Otherwise, swap all NaNs
    i = iₗ
    j = iᵤ
    N = iᵤ - iₗ
    @inbounds for n = 0:N-Nₙₐₙ
        i = iₗ + n
        if A[i] != A[i]
            while A[j] != A[j]
                j -= 1
            end
            j <= i && break
            A[i], A[j] = A[j], A[i]
            I[i], I[j] = I[j], I[i]
            j -= 1
        end
    end
    return I, A, iₗ, iᵤ - Nₙₐₙ
end
# For integers, don't need to check for NaNs
sortnans!(A::AbstractArray{<:Integer}, iₗ::Int=firstindex(A), iᵤ::Int=lastindex(A)) = A, iₗ, iᵤ
argsortnans!(I::AbstractArray, A::AbstractArray{<:Integer}, iₗ::Int=firstindex(A), iᵤ::Int=lastindex(A)) = I, A, iₗ, iᵤ

# Partially sort `A` around the `k`th sorted element and return that element
@inline function quickselect!(A::AbstractArray, iₗ=firstindex(A), iᵤ=lastindex(A), k=(iₗ+iᵤ)÷2)
    Base.Sort.partialsort!(view(A, iₗ:iᵤ), k-(iₗ-1))
end

# Sort `A`, assuming no NaNs
function quicksort!(A::TA, iₗ=firstindex(A), iᵤ=lastindex(A)) where {TA<:AbstractArray}
    # If already sorted, we're done here
    issortedrange(A, iₗ, iᵤ) && return A::TA

    # Otherwise, we have to sort
    N = iᵤ - iₗ + 1
    if isantisortedrange(A, iₗ, iᵤ)
        vreverse!(A, iₗ, iᵤ)
    elseif N == 3
        # We know we are neither sorted nor antisorted, so only four possibilities remain
        iₘ = iₗ + 1
        a,b,c = A[iₗ], A[iₘ], A[iᵤ]
        if a <= b
            if a <= c
                A[iₘ], A[iᵤ] = c, b             # a ≤ c ≤ b
            else
                A[iₗ], A[iₘ], A[iᵤ] = c, a, b   # c ≤ a ≤ b
            end
        else
            if a <= c
                A[iₗ], A[iₘ] = b, a             # b ≤ a ≤ c
            else
                A[iₗ], A[iₘ], A[iᵤ] = b, c, a   # b ≤ c ≤ a
            end
        end
    else
        # Pick a pivot for partitioning
        iₚ = iₗ + (N >> 2)
        A[iₗ], A[iₚ] = A[iₚ], A[iₗ]
        pivot = A[iₗ]

        # Count up elements that must be moved to upper partition
        Nᵤ = 0
        @inbounds @simd ivdep for i = (iₗ+1):iᵤ
            Nᵤ += A[i] >= pivot
        end
        Nₗ = N - Nᵤ

        # Swap elements between upper and lower partitions
        i = iₗ
        j = iᵤ
        @inbounds for n = 1:Nₗ-1
            i = iₗ + n
            if A[i] >= pivot
                while A[j] >= pivot
                    j -= 1
                end
                j <= i && break
                A[i], A[j] = A[j], A[i]
                j -= 1
            end
        end
        # Move pivot to the top of the lower partition
        iₚ = iₗ + Nₗ - 1
        A[iₗ], A[iₚ] = A[iₚ], A[iₗ]
        # Recurse: sort both upper and lower partitions
        quicksort!(A, iₗ, iₚ)::TA
        quicksort!(A, iₚ+1, iᵤ)::TA
    end
    return A::TA
end

# Argsort: sort A and permute I to match `A`, assuming no NaNs
function argsort!(I::TI, A::TA, iₗ::Int=firstindex(A), iᵤ::Int=lastindex(A)) where {TI<:AbstractArray, TA<:AbstractArray}
    # If already sorted, we're done here
    issortedrange(A, iₗ, iᵤ) && return (I, A)::Tuple{TI, TA}
        
    # Otherwise, we have to sort
    N = iᵤ - iₗ + 1
    if isantisortedrange(A, iₗ, iᵤ)
        vreverse!(A, iₗ, iᵤ)
        vreverse!(I, iₗ, iᵤ)
    elseif N == 3
        # We know we are neither sorted nor antisorted, so only four possibilities remain
        iₘ = iₗ + 1
        a,b,c = A[iₗ], A[iₘ], A[iᵤ]
        if a <= b
            if a <= c
                A[iₘ], A[iᵤ] = c, b             # a ≤ c ≤ b
                I[iₘ], I[iᵤ] = I[iᵤ], I[iₘ]
            else
                A[iₗ], A[iₘ], A[iᵤ] = c, a, b   # c ≤ a ≤ b
                I[iₗ], I[iₘ], I[iᵤ] = I[iᵤ], I[iₗ], I[iₘ]
            end
        else
            if a <= c
                A[iₗ], A[iₘ] = b, a             # b ≤ a ≤ c
                I[iₗ], I[iₘ] = I[iₘ], I[iₗ]
            else
                A[iₗ], A[iₘ], A[iᵤ] = b, c, a   # b ≤ c ≤ a
                I[iₗ], I[iₘ], I[iᵤ] = I[iₘ], I[iᵤ], I[iₗ]
            end
        end
    else
        # Pick a pivot for partitioning
        iₚ = iₗ + (N >> 2)
        A[iₗ], A[iₚ] = A[iₚ], A[iₗ]
        I[iₗ], I[iₚ] = I[iₚ], I[iₗ]
        pivot = A[iₗ]

        # Count up elements that must be moved to upper partition
        Nᵤ = 0
        @inbounds @simd ivdep for i = (iₗ+1):iᵤ
            Nᵤ += A[i] >= pivot
        end
        Nₗ = N - Nᵤ

        # Swap elements between upper and lower partitions
        i = iₗ
        j = iᵤ
        @inbounds for n = 1:Nₗ-1
            i = iₗ + n
            if A[i] >= pivot
                while A[j] >= pivot
                    j -= 1
                end
                j <= i && break
                A[i], A[j] = A[j], A[i]
                I[i], I[j] = I[j], I[i]
                j -= 1
            end
        end
        # Move pivot to the top of the lower partition
        iₚ = iₗ + Nₗ - 1
        A[iₗ], A[iₚ] = A[iₚ], A[iₗ]
        I[iₗ], I[iₚ] = I[iₚ], I[iₗ]
        # Recurse: sort both upper and lower partitions
        argsort!(I, A, iₗ, iₚ)::Tuple{TI, TA}
        argsort!(I, A, iₚ+1, iᵤ)::Tuple{TI, TA}
    end
    return (I, A)::Tuple{TI, TA}
end

@inline function nansort!(A)
    A, iₗ, iᵤ = sortnans!(A)
    quicksort!(A, iₗ, iᵤ)
    return A
end
export nansort!
@inline function nanargsort!(I, A)
    I, A, iₗ, iᵤ = argsortnans!(I, A)
    argsort!(I, A, iₗ, iᵤ)
    return I, A
end
export nanargsort!