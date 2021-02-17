# NaNStatistics

[![Dev][docs-dev-img]][docs-dev-url]
[![Build Status][ci-img]][ci-url]
[![codecov.io][codecov-img]][codecov-url]

Fast ([LoopVectorization](https://github.com/chriselrod/LoopVectorization.jl)-based) summary statistics, histograms, and binning – ignoring NaNs

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


[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://brenhinkeller.github.io/NaNStatistics.jl/stable/
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://brenhinkeller.github.io/NaNStatistics.jl/dev/
[ci-img]: https://github.com/brenhinkeller/NaNStatistics.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/brenhinkeller/NaNStatistics.jl/actions?query=workflow%3ACI
[codecov-img]: http://codecov.io/github/brenhinkeller/NaNStatistics.jl/coverage.svg?branch=main
[codecov-url]: http://codecov.io/github/brenhinkeller/NaNStatistics.jl?branch=main
