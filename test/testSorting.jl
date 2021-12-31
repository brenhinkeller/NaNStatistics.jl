## --- Test sorting functions directly

    A = rand(100)

    # SortNaNs
    B, iₗ, iᵤ = NaNStatistics.sortnans!(copy(A))
    @test B == A
    @test (iₗ, iᵤ) == (1, 100)

    # Quicksort
    NaNStatistics.quicksort!(A)
    sort!(B)
    @test A == B

    # Multithreaded quicksort
    A = rand(100)
    B = copy(A)
    NaNStatistics.quicksortt!(A)
    sort!(B)
    @test A == B

    # Partialsort
    A = rand(101)
    m = median(A)
    NaNStatistics.quickselect!(A, 1, 101, 51)
    @test A[51] == m

    # # Vsort, Float64
    # A = rand(100)
    # B = NaNStatistics.vsort(A, multithreaded=false)
    # @test issorted(B)
    # A = rand(100)
    # B = NaNStatistics.vsort(A, multithreaded=true)
    # @test issorted(B)
    #
    # # Vsort, Int64
    # A = rand(Int, 100)
    # B = NaNStatistics.vsort(A, multithreaded=false)
    # @test issorted(B)
    # A = rand(Int, 100)
    # B = NaNStatistics.vsort(A, multithreaded=true)
    # @test issorted(B)


## --- Test nanmedian!

    @test nanmedian!(0:10) == 5
    @test nanmedian!(1:10) == 5.5

    A = rand(100)
    @test nanmedian!(copy(A)) == median(A)

    A = rand(55,82)
    @test nanmedian!(copy(A)) == median(A)
    @test nanmedian!(copy(A), dims=1) == median(A, dims=1)
    @test nanmedian!(copy(A), dims=2) == median(A, dims=2)

    A = rand(10,11,12)
    @test nanmedian!(copy(A)) == median(A)
    @test nanmedian!(copy(A), dims=1) == median(A, dims=1)
    @test nanmedian!(copy(A), dims=2) == median(A, dims=2)
    @test nanmedian!(copy(A), dims=3) == median(A, dims=3)
    @test nanmedian!(copy(A), dims=(1,2)) == median(A, dims=(1,2))
    @test nanmedian!(copy(A), dims=(2,3)) == median(A, dims=(2,3))

## --- Test nanpctile! / nanquantile!

    @test nanpctile!(0:10, 0) == 0
    @test nanpctile!(0:10, 1) ≈ 0.1
    @test nanpctile!(0:10, 100) == 10
    @test nanpctile!(0:10, 13.582) ≈ 1.3582
    @test nanpctile!(collect(1:10), 50) == 5.5

    A = rand(100)
    @test nanpctile!(copy(A), 50) == median(A)

    A = rand(55,82)
    @test nanpctile!(copy(A), 50) == median(A)
    @test nanpctile!(copy(A), 50, dims=1) == median(A, dims=1)
    @test nanpctile!(copy(A), 50, dims=2) == median(A, dims=2)

    A = rand(10,11,12)
    @test nanpctile!(copy(A), 50) == median(A)
    @test nanpctile!(copy(A), 50, dims=1) == median(A, dims=1)
    @test nanpctile!(copy(A), 50, dims=2) == median(A, dims=2)
    @test nanpctile!(copy(A), 50, dims=3) == median(A, dims=3)
    @test nanpctile!(copy(A), 50, dims=(1,2)) == median(A, dims=(1,2))
    @test nanpctile!(copy(A), 50, dims=(2,3)) == median(A, dims=(2,3))


## ---
