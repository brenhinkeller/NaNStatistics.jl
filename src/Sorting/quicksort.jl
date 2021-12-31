# Move all NaNs to the end of the array A
function sortnans!(A, iₗ=firstindex(A), iᵤ=lastindex(A))
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

# Partially sort `A` around the `k`th sorted element and return that element
function quickselect!(A::AbstractArray, iₗ=firstindex(A), iᵤ=lastindex(A), k=(iₗ+iᵤ)÷2)
    # Pick a pivot for partitioning
    N = iᵤ - iₗ + 1
    A[iₗ], A[k] = A[k], A[iₗ]
    pivot = A[iₗ]

    # Count up elements that must be moved to upper partition
    Nᵤ = 0
    @turbo for i = iₗ:iᵤ
        Nᵤ += A[i] > pivot
    end
    Nₗ = N - Nᵤ

    # Swap elements between upper and lower partitions
    i = iₗ
    j = iᵤ
    @inbounds for n = 1:Nₗ-1
        i = iₗ + n
        if A[i] > pivot
            while A[j] > pivot
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

# Check for sortedness, assuming no NaNs
function issortedrange(A, iₗ, iᵤ)
    @inbounds for i = iₗ+1:iᵤ
        if A[i-1] > A[i]
            return false
        end
    end
    return true
end

# Sort `A`, assuming no NaNs
function quicksort!(A, iₗ=firstindex(A), iᵤ=lastindex(A))
    if issortedrange(A, iₗ, iᵤ)
        # If already sorted, we're done here
        return A
    end
    # Otherwise, we have to sort
    N = iᵤ - iₗ + 1
    if  N == 2
        # If we've gotten here, we know we're not sorted, so reverse elements
        A[iᵤ], A[iₗ] = A[iₗ], A[iᵤ]
        return A
    elseif N == 3
        # For N==3, can sort with 3 more comparisons, worst-case
        iₘ = iₗ + 1
        a,b,c = A[iₗ], A[iₘ], A[iᵤ]
        if a > b
            if b > c
                # c < b < a
                A[iₗ], A[iᵤ] = c, a
            elseif c > a
                # b < a < c
                A[iₗ], A[iₘ] = b, a
            else
                # b <= c <= a
                A[iₗ], A[iₘ], A[iᵤ] = b, c, a
            end
        else # a <= b
            if c > a
                # a < c < b
                A[iₘ], A[iᵤ] = c, b
            else
                # c <= a <= b
                A[iₗ], A[iₘ], A[iᵤ] = c, a, b
            end
        end
        return A
    else
        # Pick a pivot for partitioning
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
    if  N == 2
        # If we've gotten here, we know we're not sorted, so reverse elements
        A[iᵤ], A[iₗ] = A[iₗ], A[iᵤ]
        return A
    elseif N == 3
        # For N==3, can sort with 3 more comparisons, worst-case
        iₘ = iₗ + 1
        a,b,c = A[iₗ], A[iₘ], A[iᵤ]
        if a > b
            if b > c
                # c < b < a
                A[iₗ], A[iᵤ] = c, a
            elseif c > a
                # b < a < c
                A[iₗ], A[iₘ] = b, a
            else
                # b <= c <= a
                A[iₗ], A[iₘ], A[iᵤ] = b, c, a
            end
        else # a <= b
            if c > a
                # a < c < b
                A[iₘ], A[iᵤ] = c, b
            else
                # c <= a <= b
                A[iₗ], A[iₘ], A[iᵤ] = c, a, b
            end
        end
        return A
    else
        # Pick a pivot for partitioning
        if iᵤ-iₗ < 4096
            # Just use first element as pivot
            pivot = A[iₗ]
        else
            # Put a modicum of effort into choosing a pivot
            # This little maneuver will cost us about 32 ns
            pivot = semimedian(A, iₗ, iₗ+32)
            iₚ = findfirstinrange(A, pivot, iₗ, iₗ+32)
            A[iₗ], A[iₚ] = A[iₚ], A[iₗ]
        end

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

# Find the exact median of three elements, assuming no NaNs
function median_of_three(a,b,c)
    if a < b
        if b < c
            b
        else
            max(a,c)
        end
    else
        if b < c
            min(a,c)
        else
            b
        end
    end
end

# Find an approximate median for pivoting
function semimedian(A, iₗ=firstindex(A), iᵤ=lastindex(A))
    m = A[iₗ]
    @inbounds for i = iₗ+2:3:iᵤ
        m = median_of_three(m, A[i-1], A[i])
    end
    return m
end

# Return fist matching linear index within range
function findfirstinrange(A, target, iₗ, iᵤ)
    @inbounds for i = iₗ:iᵤ-1
        if A[i] == target
            return i
        end
    end
    return iᵤ
end
