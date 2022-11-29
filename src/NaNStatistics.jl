module NaNStatistics

    # AVX vectorziation tools
    using IfElse: ifelse
    using LoopVectorization

    using Static
    const IntOrStaticInt = Union{Integer, StaticInt}
    _dim(::Type{StaticInt{N}}) where {N} = N::Int

    include("ArrayStats/ArrayStats.jl")
    include("ArrayStats/nanmean.jl")
    include("ArrayStats/nansum.jl")
    include("ArrayStats/nancumsum.jl")
    include("ArrayStats/nanvar.jl")
    include("ArrayStats/nanstd.jl")
    include("ArrayStats/nancov.jl")
    include("Sorting/quicksort.jl")
    include("Sorting/nanmedian.jl")
    include("Sorting/nanpctile.jl")
    include("Histograms.jl")
    include("Binning.jl")

    using SnoopPrecompile
    @precompile_setup begin

        @precompile_all_calls begin
            for T in (Float64,)
                for nd in 1:3
                    A = ones(T, ntuple(i->10, nd))
                    nansum(A)
                    nanmean(A)
                    nanstd(A)
                    nanvar(A)
                    nanminimum(A)
                    nanmaximum(A)

                    for d in 1:nd
                        nansum(A, dims=d)
                        nanmean(A, dims=d)
                        nanstd(A, dims=d)
                        nanvar(A, dims=d)
                        nanminimum(A, dims=d)
                        nanmaximum(A, dims=d)
                    end
                    if nd > 1
                        for i = 2:nd
                            for j = 1:i-1
                                nansum(A, dims=(j,i))
                                nanmean(A, dims=(j,i))
                                nanstd(A, dims=(j,i))
                                nanvar(A, dims=(j,i))
                                nanminimum(A, dims=(j,i))
                                nanmaximum(A, dims=(j,i))
                            end
                        end
                    end

                end
            end
        end
    end


end
