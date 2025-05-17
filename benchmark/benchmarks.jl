using NaNStatistics
using BenchmarkTools

const SUITE = BenchmarkGroup()

# Benchmark nanmean and nansum optimizations for large arrays and ensure they
# don't effect performance on small arrays.
SUITE["nanmean"] = BenchmarkGroup()
SUITE["nansum"] = BenchmarkGroup()

big::Array{Float64, 3} = rand(1000, 1000, 10)
for i in rand(eachindex(big), 1000)
    big[i] = NaN
end

small::Array{Float64, 3} = rand(10, 10, 10)
for i in rand(eachindex(small), 50)
    small[i] = NaN
end

small_vector::Vector{Float64} = rand(10)

for i in 1:3
    SUITE["nanmean"]["big array, dim $(i)"] = @benchmarkable nanmean($big; dim=$i)
    SUITE["nanmean"]["small array, dim $(i)"] = @benchmarkable nanmean($small; dim=$i)
    SUITE["nanmean"]["small vector"] = @benchmarkable nanmean($small_vector)
    SUITE["nansum"]["big array, dim $(i)"] = @benchmarkable nansum($big; dim=$i)
    SUITE["nansum"]["small array, dim $(i)"] = @benchmarkable nansum($small; dim=$i)
    SUITE["nansum"]["small vector"] = @benchmarkable nansum($small_vector)
end
