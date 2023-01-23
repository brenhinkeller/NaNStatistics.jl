# Move all NaNs to the end of the array A
function sortnans!(A, iₗ=firstindex(A), iᵤ=lastindex(A))
    # Return early if range is empty
    iₗ >= iᵤ && return A, iₗ, iᵤ

    # Count up NaNs
    Nₙₐₙ = 0
    @turbo for i = iₗ:iᵤ
        Nₙₐₙ += A[i] != A[i]
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
# For integers, don't need to check for NaNs
sortnans!(A::AbstractArray{<:Integer}, iₗ=firstindex(A), iᵤ=lastindex(A)) = A, iₗ, iᵤ

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
        @turbo for i ∈ 0:n
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

# Partially sort `A` around the `k`th sorted element and return that element
function quickselect!(A::AbstractArray, iₗ=firstindex(A), iᵤ=lastindex(A), k=(iₗ+iᵤ)÷2)
    # Pick a pivot for partitioning
    N = iᵤ - iₗ + 1
    A[iₗ], A[k] = A[k], A[iₗ]
    pivot = A[iₗ]

    # Count up elements that must be moved to upper partition
    Nᵤ = 0
    @turbo for i = (iₗ+1):iᵤ
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
    # Recurse: select from partition containing k
    (iₗ <= k < iₚ) && quickselect!(A, iₗ, iₚ, k)
    (iₚ < k <= iᵤ) && quickselect!(A, iₚ+1, iᵤ, k)
    return A[k]
end


# Sort `A`, assuming no NaNs
function quicksort!(A, iₗ=firstindex(A), iᵤ=lastindex(A))
    if issortedrange(A, iₗ, iᵤ)
        # If already sorted, we're done here
        return A
    end
    # Otherwise, we have to sort
    N = iᵤ - iₗ + 1
    if isantisortedrange(A, iₗ, iᵤ)
        vreverse!(A, iₗ, iᵤ)
        return A
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
        return A
    else
        # Pick a pivot for partitioning
        iₚ = iₗ + (N >> 2)
        A[iₗ], A[iₚ] = A[iₚ], A[iₗ]
        pivot = A[iₗ]

        # Count up elements that must be moved to upper partition
        Nᵤ = 0
        @turbo for i = (iₗ+1):iᵤ
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
        quicksort!(A, iₗ, iₚ)
        quicksort!(A, iₚ+1, iᵤ)
    end
end

# Sort `A`, assuming no NaNs, multithreaded
function quicksortt!(A, iₗ=firstindex(A), iᵤ=lastindex(A), level=1)
    if issortedrange(A, iₗ, iᵤ)
        # If already sorted, we're done here
        return A
    end
    # Otherwise, we have to sort
    N = iᵤ - iₗ + 1
    if isantisortedrange(A, iₗ, iᵤ)
        vreverse!(A, iₗ, iᵤ)
        return A
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
        return A
    else
        # Pick a pivot for partitioning
        iₚ = iₗ + (N >> 2)
        A[iₗ], A[iₚ] = A[iₚ], A[iₗ]
        pivot = A[iₗ]

        # Count up elements that must be moved to upper partition
        Nᵤ = 0
        @turbo for i = (iₗ+1):iᵤ
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
        if level < 7
            @sync begin
                Threads.@spawn quicksortt!(A, iₗ, iₚ, level+1)
                Threads.@spawn quicksortt!(A, iₚ+1, iᵤ, level+1)
            end
        else
            quicksort!(A, iₗ, iₚ)
            quicksort!(A, iₚ+1, iᵤ)
        end
        return A
    end
end
