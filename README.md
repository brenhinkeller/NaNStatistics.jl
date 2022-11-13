# NaNStatistics
[![Docs][docs-dev-img]][docs-dev-url]
[![CI][ci-img]][ci-url]
[![Coverage][codecov-img]][codecov-url]

*Because* `NaN` *is just* `missing` *with hardware support!*

Fast (often [LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl)-based) summary statistics, histograms, and binning — all ignoring `NaN`s, as if `NaN` represented missing data.

See also [JuliaSIMD/VectorizedStatistics.jl](https://github.com/JuliaSIMD/LoopVectorization.jl) for similar vectorized implementations that don't ignore` `NaN`s

### Summary statistics
Summary statistics exported by NaNStatistics are generally named the same as their normal counterparts, but with "nan" in front of the name, similar to the Matlab and NumPy conventions. Options include:
##### Reductions
* `nansum`
* `nanminimum`
* `nanmaximum`
* `nanextrema`

##### Measures of central tendency
* `nanmean` &emsp; arithmetic mean, ignoring `NaN`s
* `nanmedian` &emsp; median, ignoring `NaN`s
* `nanmedian!` &emsp; as `nanmedian` but quicksorts in-place for efficiency

##### Measures of dispersion
* `nanvar` &emsp; variance
* `nanstd` &emsp; standard deviation
* `nancov` &emsp; covariance
* `nancor` &emsp; Pearson's product-moment correlation
* `nanaad` &emsp; mean (average) absolute deviation from the mean
* `nanmad` &emsp; median absolute deviation from the median
* `nanmad!` &emsp; as `nanmad` but quicksorts in-place for efficiency
* `nanrange` &emsp; range between nanmaximum and nanminimum
* `nanpctile` &emsp; percentile
* `nanpctile!` &emsp; as `nanpctile` but quicksorts in-place for efficiency

Note that, regardless of implementation, functions involving medians or percentiles are generally significantly slower than other summary statistics, since calculating a median or percentile requires a quicksort or quickselect of the input array; if not done in-place as in `nanmedian!` and `nanpctile!` then a copy of the entire array must also be made.

These functions will generally support the same `dims` keyword argument as their normal Julia counterparts (though are most efficient when operating on an entire collection).
As an alternative to `dims`, the `dim` keyword is also supported, which behaves identially to `dims` except that it also (as is the norm in some other languages) drops any singleton dimensions that have been reduced over.
```
julia> a = rand(100000);

julia> @btime minimum($a)
  51.950 μs (0 allocations: 0 bytes)
7.630517166790085e-6

julia> using NaNStatistics

julia> @btime nanminimum($a)
  19.690 μs (0 allocations: 0 bytes)
7.630517166790085e-6

julia> a[rand(1:100000, 10000)] .= NaN;

julia> @btime nanminimum($a)
  19.663 μs (0 allocations: 0 bytes)
7.630517166790085e-6
```
### Histograms
The main 1D and 2D histogram function is `histcounts` (with an in-place variant `histcounts!`), and will, as you might expect for this package, ignore NaNs. However, it might be worth using for speed even if your data don't contain any NaNs:
```
julia> b = 10 * rand(100000);

julia> using StatsBase

julia> @btime fit(Histogram, $b, 0:1:10, closed=:right)
  2.633 ms (2 allocations: 224 bytes)
Histogram{Int64, 1, Tuple{StepRange{Int64, Int64}}}
edges:
  0:1:10
weights: [10128, 10130, 10084, 9860, 9973, 10062, 10003, 10045, 9893, 9822]
closed: right
isdensity: false

julia> using NaNStatistics

julia> @btime histcounts($b, 0:1:10)
  1.037 ms (1 allocation: 160 bytes)
10-element Vector{Int64}:
 10128
 10130
 10084
  9860
  9973
 10062
 10003
 10045
  9893
  9822
```

### Binning
NaNStatistics also provides functions that will efficiently calculate the summary statistics of a given dependent variable `y` binned by an independent variable `x`. These currently include:
* `nanbinmean` / `nanbinmean!`
* `nanbinmedian` / `nanbinmedian!`
```
julia> x = 10 * rand(100000);

julia> y = x.^2 .+ randn.();

julia> xmin, xmax, nbins = 0, 10, 10;

julia> @btime nanbinmean($x,$y,xmin,xmax,nbins)
  364.082 μs (2 allocations: 320 bytes)
10-element Vector{Float64}:
  0.3421697507351903
  2.3065542448799015
  6.322448227456871
 12.340306767007629
 20.353233411797074
 30.347815506059405
 42.31866909140384
 56.32256214256441
 72.35387230251672
 90.35682945641588
```
### Other functions
* `movmean`
A simple moving average function, which can operate in 1D or 2D, ignoring NaNs.
```
julia> A = rand(1:10, 4,4)
4×4 Matrix{Int64}:
 3  5  10  3
 4  2   5  8
 5  6   8  8
 2  6  10  6

julia> movmean(A, 3)
4×4 Matrix{Float64}:
 3.5      4.83333  5.5      6.5
 4.16667  5.33333  6.11111  7.0
 4.16667  5.33333  6.55556  7.5
 4.75     6.16667  7.33333  8.0
 ```

 * `nanstandardize` / `nanstandardize!`
 De-mean and set to unit variance


[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://brenhinkeller.github.io/NaNStatistics.jl/stable/
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://brenhinkeller.github.io/NaNStatistics.jl/dev/
[ci-img]: https://github.com/brenhinkeller/NaNStatistics.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/brenhinkeller/NaNStatistics.jl/actions?query=workflow%3ACI
[codecov-img]: https://codecov.io/gh/brenhinkeller/NaNStatistics.jl/branch/main/graph/badge.svg
[codecov-url]: http://codecov.io/github/brenhinkeller/NaNStatistics.jl?branch=main
