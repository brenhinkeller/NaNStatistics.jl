module NaNStatisticsUnitfulExt

    using NaNStatistics, Unitful

    # With Unitful quantities, the assumption that a numeric value will have the same type 
    # as that value squared does not hold, so several methods must be adapted

    # Functions to strip units if and only if neccessary
    unitless(x) = x
    unitless(x::Unitful.Quantity) = Unitful.ustrip(x)
    unitless(x::AbstractArray{<:Unitful.Quantity}) = Unitful.ustrip(x)

    # Cannot take in-place square root
    NaNStatistics.sqrt!(A::AbstractArray{<:Unitful.Quantity}) = sqrt.(A)
    
    # Mean has different units than variance and standard error
    function NaNStatistics.__nanvar(mean, corrected, A::AbstractArray{T}, region::Union{Int, Dims}, ::Colon) where {T<:Unitful.Quantity}
        u = Unitful.unit(T)
        return NaNStatistics.__nanvar(unitless(mean), corrected, unitless(A), region, :) .* u^2
    end
    function NaNStatistics.__nansem(mean, corrected, A::AbstractArray{T}, region::Union{Int, Dims}, ::Colon) where {T<:Unitful.Quantity}
        u = Unitful.unit(T)
        return NaNStatistics.__nansem(unitless(mean), corrected, unitless(A), region, :) .* u
    end
    
    # Skewness and kurtosis are unitless
    function NaNStatistics.__nanskewness(mean, corrected, A::AbstractArray{T}, region::Union{Int, Dims}, ::Colon) where {T<:Unitful.Quantity}
        return NaNStatistics.__nanskewness(unitless(mean), corrected, unitless(A), region, :)
    end
    function NaNStatistics.__nankurtosis(mean, corrected, A::AbstractArray{T}, region::Union{Int, Dims}, ::Colon) where {T<:Unitful.Quantity}
        return NaNStatistics.__nankurtosis(unitless(mean), corrected, unitless(A), region, :)
    end
end