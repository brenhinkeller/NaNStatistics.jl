module NaNStatisticsHwlocExt

import Hwloc
import NaNStatistics

const cache_sizes = Ref{@NamedTuple{L1::Int, L2::Int, L3::Int}}((; L1=0, L2=0, L3=0))

NaNStatistics.get_size_threshold(x::Symbol) = cache_sizes[][x]

function __init__()
    # Not `Hwloc.cachesize()`: it throws when a level is missing (no L3 on Apple
    # silicon) and warns when sizes differ across cores (P/E cores)
    sizes = map(f -> maximum(f(); init=0), (Hwloc.l1cache_sizes, Hwloc.l2cache_sizes, Hwloc.l3cache_sizes))
    # A missing level falls back to the largest cache below it, so `:L3` stays a
    # meaningful threshold instead of 0
    L1, L2, L3 = accumulate(max, sizes)
    cache_sizes[] = (; L1, L2, L3)
end

end
