x = DimArray(rand(5, 10), (; foo=1:5, bar=sort(rand(UInt8, 10))))

# `res` is a global variable with a fixed type. We assign all results to it to
# ensure that they do in fact return a DimArray rather than an Array
res::DimArray = DimArray([1], :foo)

@testset "Reduction functions" begin
    global res

    # Test `dim`
    res = nanmean(x; dim=:foo)
    @test res isa DimArray
    @test hasdim(res, (:foo, :bar)) == (false, true)
    @test res == nanmean(parent(x); dim=1)

    # Test `dims`
    res = nanmean(x; dims=:foo)
    @test hasdim(res, (:foo, :bar)) == (true, true)

    # Test that extra args are passed through
    res = nanpctile(x, 50; dims=:foo)
    @test res isa DimArray
    @test res == nanpctile(parent(x), 50; dims=1)

    # Test that extra kwargs are passed through
    res = nansem(x; dims=:foo, corrected=false)
    @test res == nansem(parent(x); dims=1, corrected=false)

    # Test nanstd()'s with weights
    res = nanstd(x, x; dims=:foo)
    @test res == nanstd(parent(x), parent(x); dims=1)

    res = zeros(Dim{:bar}(size(x, :bar)))

    # Since all the reduction functions use the same machinery we only do sanity
    # checks for the other ones.
    @test nanmean(x; dim=:foo) isa DimArray
    @test nanmean!(res, x; dim=:foo) isa DimArray
    @test nansum(x; dim=:foo) isa DimArray
    @test nansum!(res, x; dim=:foo) isa DimArray
    @test nanstd(x; dim=:foo) isa DimArray
    @test nanvar(x; dim=:foo) isa DimArray
    @test nanaad(x; dim=:foo) isa DimArray
    @test nansem(x; dim=:foo) isa DimArray
    @test nanskewness(x; dim=:foo) isa DimArray
    @test nankurtosis(x; dim=:foo) isa DimArray
    @test nanmad(x; dim=:foo) isa DimArray
    @test nanmad!(x; dim=:foo) isa DimArray
    @test nanmedian(x; dim=:foo) isa DimArray
    @test nanmedian!(x; dim=:foo) isa DimArray
    @test nanpctile(x, 50; dim=:foo) isa DimArray
    @test nanpctile!(x, 50; dim=:foo) isa DimArray
    @test nanquantile(x, 0.5; dims=:foo) isa DimArray
    @test nanquantile!(x, 0.5; dims=:foo) isa DimArray
end

@testset "Functions with manual support" begin
    global res

    res = nancumsum(x; dims=:foo)
    @test res == nancumsum(parent(x); dims=1)

    res = nancumsum!(x; dims=:foo)
    @test res == nancumsum!(parent(x); dims=1)

    res = movmean(x, 2)
    @test res == movmean(parent(x), 2)
end

# These functions magically work already \o/
@testset "Functions with implicit support" begin
    global res

    res = movmean(x[foo=1], (1, 2))
    @test res == movmean(parent(x[foo=1]), (1, 2))

    res = movsum(x, 2)
    @test res == movsum(parent(x), 2)

    res = nanextrema(x; dim=:foo)
    @test res == nanextrema(parent(x); dim=1)

    res = nanmaximum(x; dim=:foo)
    @test res == nanmaximum(parent(x); dim=1)

    res = nanminimum(x; dim=:foo)
    @test res == nanminimum(parent(x); dim=1)

    res = nanrange(x; dim=:foo)
    @test res == nanrange(parent(x); dim=1)
end

@testset "Allocation helpers" begin
    data = rand(X(5), Y(11:15))

    # Just test allocate_nanmean() since the other reduction allocators all go
    # through _allocate_reduce().
    out = NaNStatistics.allocate_nanmean(data, 1)
    @test size(out) == (1, 5)
    @test out isa DimMatrix{Float64}
    @test lookup(out, Y) == lookup(data, Y)
end
