module NaNStatisticsHwlocExt

import Hwloc
import NaNStatistics

cache_sizes::@NamedTuple{L1::Int, L2::Int, L3::Int} = (; L1=0, L2=0, L3=0)

NaNStatistics.get_size_threshold(x::Symbol) = cache_sizes[x]

function __init__()
    # Not `Hwloc.cachesize()`: it throws when a level is missing (no L3 on Apple
    # silicon) and warns when sizes differ across cores (P/E cores)
    L1, L2, L3 = (maximum(f(); init=0) for f in (Hwloc.l1cache_sizes, Hwloc.l2cache_sizes, Hwloc.l3cache_sizes))
    global cache_sizes = (; L1, L2, L3)
    L1 > 0 && (NaNStatistics.NANMEAN_SIZE_THRESHOLD = :L1)
end

end
