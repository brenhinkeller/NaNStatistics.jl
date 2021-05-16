
function histcounts(x::AbstractArray, xedges::AbstractRange; T=Int64)
    N = fill(zero(T), length(xedges)-1)
    histcounts!(N, x, xedges)
    return N
end
histcounts(x, xmin::Number, xmax::Number, nbins::Integer; T=Int64) = histcounts(x, range(xmin, xmax, length=nbins+1); T=T)
export histcounts

function histcounts!(N::Array, x::AbstractArray, xedges::AbstractRange)
    # What is the size of each bin?
    nbins = length(xedges) - 1
    xmin, xmax = extrema(xedges)
    δiδx = nbins/(xmax-xmin)

    # Make sure we don't have a segfault by filling beyond the length of N
    # in the @inbounds loop below
    if length(N) < nbins
        nbins = length(N)
        @warn "length(N) < nbins; any bins beyond length(N) will not be filled"
    end

    # Loop through each element of x
    @inbounds for i ∈ eachindex(x)
        xᵢ = x[i]
        if xᵢ==xᵢ # If not a NaN
            δi = (xᵢ - xmin) * δiδx
            if 0 < δi <= nbins
                binindex = ceil(Int, δi)
                N[binindex] += 1
            end
        end
    end
    return N
end
function histcounts!(N::Array, x::AbstractArray{<:Integer}, xedges::AbstractRange)
    # What is the size of each bin?
    nbins = length(xedges) - 1
    xmin, xmax = extrema(xedges)
    δiδx = nbins/(xmax-xmin)

    # Make sure we don't have a segfault by filling beyond the length of N
    # in the @inbounds loop below
    if length(N) < nbins
        nbins = length(N)
        @warn "length(N) < nbins; any bins beyond length(N) will not be filled"
    end

    # Loop through each element of x
    @inbounds for i ∈ eachindex(x)
        xᵢ = x[i]
        δi = (xᵢ - xmin) * δiδx
        if 0 < δi <= nbins
            binindex = ceil(Int, δi)
            N[binindex] += 1
        end
    end
    return N
end
histcounts!(N, x, xmin::Number, xmax::Number, nbins::Integer; T=Int64) = histcounts!(N, x, range(xmin, xmax, length=nbins+1); T=T)
export histcounts!
