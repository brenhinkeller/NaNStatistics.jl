module NaNStatisticsHwlocExt

import Hwloc
import NaNStatistics

cache_sizes::@NamedTuple{L1::Int, L2::Int, L3::Int} = (; L1=0, L2=0, L3=0)

NaNStatistics.get_size_threshold(x::Symbol) = cache_sizes[x]

function __init__()
    global cache_sizes = Hwloc.cachesize()
    NaNStatistics.NANMEAN_SIZE_THRESHOLD = :L1
end

end
