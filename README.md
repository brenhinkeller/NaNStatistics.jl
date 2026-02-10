# NaNStatistics
[![Docs][docs-dev-img]][docs-dev-url]
[![CI][ci-img]][ci-url]
[![CI (julia nightly)][ci-nightly-img]][ci-nightly-url]
[![Coverage][codecov-img]][codecov-url]

*Because* `NaN` *is just* `missing` *with hardware support!*

Fast summary statistics, histograms, and binning — all ignoring `NaN`s, as if `NaN` represented missing data.

See also [JuliaSIMD/VectorizedStatistics.jl](https://github.com/JuliaSIMD/VectorizedStatistics.jl) for similar vectorized implementations that don't ignore `NaN`s.

### Summary statistics
Summary statistics exported by NaNStatistics are generally named the same as their normal counterparts, but with "nan" in front of the name, similar to the Matlab and NumPy conventions. Options include:
##### Reductions
* `nansum`
* `nansum!`
* `nanminimum`
* `nanmaximum`
* `nanextrema`

##### Measures of central tendency
* `nanmean` &emsp; arithmetic mean, ignoring `NaN`s
* `nanmean!`&emsp; as `nanmean`, but writes to a given output array
* `nanmedian` &emsp; median, ignoring `NaN`s
* `nanmedian!` &emsp; as `nanmedian` but quicksorts in-place for efficiency

##### Measures of dispersion
* `nanvar` &emsp; variance
* `nanstd` &emsp; standard deviation
* `nansem` &emsp; standard error of the mean
* `nancov` &emsp; covariance
* `nancor` &emsp; Pearson's product-moment correlation
* `nanaad` &emsp; mean (average) absolute deviation from the mean
* `nanmad` &emsp; median absolute deviation from the median
* `nanmad!` &emsp; as `nanmad` but quicksorts in-place for efficiency
* `nanrange` &emsp; range between nanmaximum and nanminimum
* `nanpctile` &emsp; percentile
* `nanpctile!` &emsp; as `nanpctile` but quicksorts in-place for efficiency

##### Other summary statistics
* `nanskewness` &emsp; skewness
* `nankurtosis` &emsp; excess kurtosis

Note that, regardless of implementation, functions involving medians or percentiles are generally significantly slower than other summary statistics, since calculating a median or percentile requires a quicksort or quickselect of the input array; if not done in-place as in `nanmedian!` and `nanpctile!` then a copy of the entire array must also be made.

These functions will generally support the same `dims` keyword argument as their normal Julia counterparts (though are most efficient when operating on an entire collection).
As an alternative to `dims`, the `dim` keyword is also supported, which behaves identially to `dims` except that it also (as is the norm in some other languages) drops any singleton dimensions that have been reduced over.
```
julia> a = rand(100000);

julia> minimum(a)
9.70221275542471e-7

julia> using NaNStatistics

julia> nanminimum(a)
9.70221275542471e-7

julia> a[rand(1:100000, 10000)] .= NaN;

julia> nanminimum(a)
7.630517166790085e-6
```
### Histograms
The main 1D and 2D histogram function is `histcounts` (with an in-place variant `histcounts!`), and will, as you might expect for this package, ignore NaNs. However, it might be worth using for speed even if your data don't contain any NaNs:
```
julia> b = 10 * rand(100000);

julia> using StatsBase

julia> @btime fit(Histogram, $b, 0:1:10, closed=:right)
  526.750 μs (2 allocations: 208 bytes)
Histogram{Int64, 1, Tuple{StepRange{Int64, Int64}}}
edges:
  0:1:10
weights: [10042, 10105, 9976, 9980, 10073, 10038, 9983, 9802, 10056, 9945]
closed: right
isdensity: false

julia> using NaNStatistics

julia> @btime histcounts($b, 0:1:10)
  155.083 μs (2 allocations: 176 bytes)
10-element Vector{Int64}:
 10042
 10105
  9976
  9980
 10073
 10038
  9983
  9802
 10056
  9945
```
(Timings as of Julia v1.10.4, NaNStatistics v0.6.36, Apple M1 Max)

In addition, several functions are provided to estimate the summary statistics of a dataset from its histogram, specifically
* `histmean` &emsp; arithmetic mean
* `histvar` &emsp; variance
* `histstd` &emsp; standard deviation
* `histskewness` &emsp; skewness
* `histkurtosis` &emsp; excess kurtosis

### Binning
NaNStatistics also provides functions that will efficiently calculate the summary statistics of a given dependent variable `y` binned by an independent variable `x`. These currently include:
* `nanbinmean` / `nanbinmean!`
* `nanbinmedian` / `nanbinmedian!`
```
julia> x = 10 * rand(100000);

julia> y = x.^2 .+ randn.();

julia> xmin, xmax, nbins = 0, 10, 10;

julia> @btime nanbinmean($x,$y,xmin,xmax,nbins)
  222.542 μs (2 allocations: 288 bytes)
10-element Vector{Float64}:
  0.3482167982440996
  2.32463720126215
  6.348942343257478
 12.352990978599395
 20.34955219534221
 30.31123519946431
 42.3578375163112
 56.33841854482159
 72.23884588251572
 90.30275863080671
```
### Other functions
* `movmean` / `movmean!`
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

### Allocation functions
To use mutating functions like `nanmean!` you can call the appropriate
allocation function and get back an array that can be passed as the output
argument.

* `allocate_nanmean`
* `allocate_nansum`
* `allocate_movmean`

### DimensionalData support
Almost all functions support
[DimArrays](https://rafaqz.github.io/DimensionalData.jl/stable/dimarrays) and
will preserve array metadata like dimensions and lookups etc. The ones that
currently don't are:
- `nancor` and `nancov` (does work but won't preserve array metadata)
- `nanstandardize` and `nanstandardize!` (completely unsupported)

### Benchmarks
We maintain a few benchmarks in `benchmark/benchmarks.jl` for use with
[AirspeedVelocity.jl](https://astroautomata.com/AirspeedVelocity.jl). If you're
developing the package and want to benchmark the current state of the code vs
`main`, install AirspeedVelocity.jl and run `benchpkg`. See the
AirspeedVelocity.jl docs for more information.

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://brenhinkeller.github.io/NaNStatistics.jl/stable/
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://brenhinkeller.github.io/NaNStatistics.jl/dev/
[ci-img]: https://github.com/brenhinkeller/NaNStatistics.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/brenhinkeller/NaNStatistics.jl/actions?query=workflow%3ACI
[ci-nightly-img]:https://github.com/brenhinkeller/NaNStatistics.jl/workflows/CI%20(Julia%20nightly)/badge.svg
[ci-nightly-url]:https://github.com/brenhinkeller/NaNStatistics.jl/actions/workflows/CI-julia-nightly.yml
[codecov-img]: https://codecov.io/gh/brenhinkeller/NaNStatistics.jl/branch/main/graph/badge.svg
[codecov-url]: http://codecov.io/github/brenhinkeller/NaNStatistics.jl?branch=main
