module NaNStatisticsDimensionalDataExt

import DimensionalData as DD
using DimensionalData: AbstractDimArray, rebuild, dims, dimnum

import NaNStatistics


# Crete methods for all the reduction functions
for func in (:nanmean, :nansum, :nanstd, :nanvar,
             :nanaad, :nansem, :nanskewness, :nankurtosis,
             :nanmad, :nanmad!,
             :nanmedian, :nanmedian!,
             :nanpctile, :nanpctile!,
             :nanquantile, :nanquantile!)

    # First we define low-level private functions that handle the three different
    # reduction cases:
    # - Reduce to a scalar
    # - Reduce and keep all dimensions
    # - Reduce and drop singleton dimensions
    func_impl = Symbol('_', func)
    @eval begin
        # Reduce to a scalar
        $func_impl(A, ::Colon, ::Colon, args...; kwargs...) = NaNStatistics.$func(parent(A), args...; kwargs...)

        # Reduce and keep all dimensions
        function $func_impl(A, dims, ::Colon, args...; kwargs...)
            data = NaNStatistics.$func(parent(A), args...; dims=dimnum(A, dims), kwargs...)
            rebuild(A, data, DD.reducedims(A, dims))
        end

        # Reduce and drop singleton dimensions
        function $func_impl(A, ::Colon, dim, args...; kwargs...)
            data = NaNStatistics.$func(parent(A), args...; dim=dimnum(A, dim), kwargs...)
            rebuild(A, data, DD.otherdims(A, dim))
        end
    end

    # And then we define the public methods. Some of these need to be special-cased.
    if func ∈ (:nanpctile, :nanpctile!, :nanquantile, :nanquantile!)
        # These functions need their second argument explicitly type-annotated
        # to avoid method ambiguities.
        @eval NaNStatistics.$func(A::AbstractDimArray, n::Number; dims=:, dim=:) = $func_impl(A, dims, dim, n)
    else
        @eval NaNStatistics.$func(A::AbstractDimArray, args...; dims=:, dim=:, kwargs...) = $func_impl(A, dims, dim, args...; kwargs...)
    end

    if func === :nanstd
        # nanstd() has a method for passing weights and the weights need to be a
        # regular array to avoid internal broadcasting errors, so we define an
        # additional method for cases where both the data and weights are
        # AbstractDimArray's.
        @eval NaNStatistics.$func(A::AbstractDimArray, W::AbstractDimArray; dims=:, dim=:) = $func_impl(A, dims, dim, parent(W))
    end
end

# Instead of explicitly rebuilding we take advantage of DimArray's built-in
# support for the array interface for nancumsum() to work. The only thing we
# have to do is convert `dims` to indices so that users can pass in dimensions
# by name.
NaNStatistics.nancumsum(A::AbstractDimArray; dims=:) = NaNStatistics._nancumsum(A, dimnum(A, dims))
NaNStatistics.nancumsum!(A::AbstractDimArray; dims=:) = NaNStatistics._nancumsum(A, dimnum(A, dims))

function NaNStatistics.movmean(A::DD.AbstractDimVecOrMat, n::Number)
    data = NaNStatistics.movmean(parent(A), n)
    rebuild(A, data)
end

end
