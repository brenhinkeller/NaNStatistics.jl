# NaNStatistics
*Because* `NaN` *is just* `missing` *with hardware support!*

[![Dev][docs-dev-img]][docs-dev-url]
[![Build Status][ci-img]][ci-url]
[![codecov.io][codecov-img]][codecov-url]

Fast (often [LoopVectorization](https://github.com/JuliaSIMD/LoopVectorization.jl)-based) summary statistics, histograms, and binning – ignoring NaNs

### Summary statistics
Summary statistics exported by NaNStatistics are generally named the same as their normal counterparts, but with "nan" in front of the name, similar to the matlab and numpy conventions. Options include:
* `nansum`
* `nanmean`
* `nanmedian`
* `nanpctile`
* `nanstd` (standard deviation)
* `nanmad` (median absolute deviation from the median)
* `nanaad` (mean (average) absolute deviation from the mean)
* `nanminimum`
* `nanmaximum`
* `nanextrema`
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
The main histogram function is `histcounts` (with an in-place variant `histcounts!`), and will, as you might expect for this package, ignore NaNs. However, it might be worth using for speed even if your data don't contain any NaNs:
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

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://brenhinkeller.github.io/NaNStatistics.jl/stable/
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://brenhinkeller.github.io/NaNStatistics.jl/dev/
[ci-img]: https://github.com/brenhinkeller/NaNStatistics.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/brenhinkeller/NaNStatistics.jl/actions?query=workflow%3ACI
[codecov-img]: http://codecov.io/github/brenhinkeller/NaNStatistics.jl/coverage.svg?branch=main
[codecov-url]: http://codecov.io/github/brenhinkeller/NaNStatistics.jl?branch=main
