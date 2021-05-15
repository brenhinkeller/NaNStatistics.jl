
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
    binwidth = (xmax-xmin)/nbins

    # Loop through each element of x
    @inbounds for i ∈ eachindex(x)
        xᵢ = x[i]
        if xᵢ==xᵢ # If not a NaN
            binindex = Integer((xᵢ - xmin) ÷ binwidth) + 1
            if 1 <= binindex <= nbins
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
    binwidth = (xmax-xmin)/nbins

    # Loop through each element of x
    @inbounds for i ∈ eachindex(x)
        xᵢ = x[i]
        binindex = Integer((xᵢ - xmin) ÷ binwidth) + 1
        if 1 <= binindex <= nbins
            N[binindex] += 1
        end
    end
    return N
end
histcounts!(N, x, xmin::Number, xmax::Number, nbins::Integer; T=Int64) = histcounts!(N, x, range(xmin, xmax, length=nbins+1); T=T)
export histcounts!
