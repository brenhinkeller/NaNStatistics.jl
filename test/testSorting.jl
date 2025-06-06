## --- Test internal sorting functions directly

    A = rand(100)
    Ix = collect(1:100)

    # SortNaNs
    B, iₗ, iᵤ = NaNStatistics.sortnans!(copy(A))
    @test B == A
    @test (iₗ, iᵤ) == (1, 100)

    Ix, B, iₗ, iᵤ = NaNStatistics.argsortnans!(Ix, copy(A))
    @test B == A
    @test Ix == 1:100
    @test (iₗ, iᵤ) == (1, 100)

    A[1:10:100] .= NaN
    B, iₗ, iᵤ = NaNStatistics.sortnans!(copy(A))
    @test (iₗ, iᵤ) == (1, 90)
    @test sort!(B[1:90]) == sort!([A[2:10]; A[12:20]; A[22:30]; A[32:40]; A[42:50]; A[52:60]; A[62:70]; A[72:80]; A[82:90]; A[92:100];])
    @test isequal(B[91:100], fill(NaN,10))

    A = rand(100)
    A[1:10:100] .= NaN
    Ix = collect(1:100)

    Ix, B, iₗ, iᵤ = NaNStatistics.argsortnans!(Ix, copy(A))
    @test (iₗ, iᵤ) == (1, 90)
    @test sort!(B[1:90]) == sort!([A[2:10]; A[12:20]; A[22:30]; A[32:40]; A[42:50]; A[52:60]; A[62:70]; A[72:80]; A[82:90]; A[92:100];])
    @test isequal(B[91:100], fill(NaN,10))

    # Quicksort
    A = rand(100)
    NaNStatistics.quicksort!(A)
    B = sort(A)
    @test A == B

    A = rand(1_000)
    B = sort(A)
    NaNStatistics.quicksort!(A)
    @test A == B

    A = rand(1_000_000)
    B = sort(A)
    NaNStatistics.quicksort!(A)
    @test A == B

    # Quicksort of already-sorted arrays
    @test NaNStatistics.quicksort!(collect(1:100)) == 1:100
    @test NaNStatistics.quicksort!(collect(100:-1:1)) == 1:100

    # Partialsort
    A = rand(101)
    m = median(A)
    NaNStatistics.quickselect!(A, 1, 101, 51)
    @test A[51] == m

    A = rand(1_000_001)
    m = median(A)
    NaNStatistics.quickselect!(A, 1, 1_000_001, 500_001)
    @test A[500_001] == m

## --- Test `nansort!` and `nanargsort!`

   # Float64 cases
   A = rand(100)
   Ix_known = sortperm(A)
   Ix = collect(1:length(A))
   nanargsort!(Ix, A)
   @test issorted(A)
   @test Ix_known == Ix
   A[rand(1:100, 10)] .= NaN
   nanargsort!(Ix, A)
   @test issorted(A)
   A .= rand.()
   nansort!(A)
   @test issorted(A)
   A .= rand.()
   A[rand(1:100, 10)] .= NaN
   nansort!(A)
   @test issorted(A)

   A = rand(10_000)
   Ix_known = sortperm(A)
   Ix = collect(1:length(A))
   nanargsort!(Ix, A)
   @test issorted(A)
   @test Ix_known == Ix
   A .= rand.()
   nansort!(A)
   @test issorted(A)
   A .= rand.()
   A[rand(1:10_000, 1_000)] .= NaN
   nansort!(A)
   @test issorted(A)

   # Int64 cases
   A = rand(Int, 100)
   B = nansort!(copy(A))
   @test issorted(B)
   Ix_known = sortperm(A)
   Ix = collect(1:length(A))
   nanargsort!(Ix, A)
   @test issorted(A)
   @test Ix_known == Ix

## --- Test nanmedian!

    @test nanmedian!(0:10) == 5
    @test nanmedian!(1:10) == 5.5
    @test nanmedian!(fill(NaN, 10)) === NaN

    A = rand(100)
    @test nanmedian!(copy(A)) == median(A)
    A[rand(1:100, 20)] .= NaN; B = A[.!isnan.(A)];
    @test nanmedian!(A) == median(B)
    A = rand(10_000)
    @test nanmedian!(copy(A)) == median(A)
    A[rand(1:10_000, 100)] .= NaN; B = A[.!isnan.(A)]
    @test nanmedian!(A) == median(B)
    A = rand(100_000)
    @test nanmedian!(copy(A)) == median(A)
    A[rand(1:100_000, 100)] .= NaN; B = A[.!isnan.(A)]
    @test nanmedian!(A) == median(B)

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
    @test nanpctile!(fill(NaN, 10), 50) === NaN

    A = rand(100)
    @test nanpctile!(copy(A), 50) == median(A)
    A[rand(1:100, 20)] .= NaN; B = A[.!isnan.(A)];
    @test nanpctile!(A, 50) == median(B)
    A = rand(10_000)
    @test nanpctile!(copy(A), 50) == median(A)
    A[rand(1:10_000, 100)] .= NaN; B = A[.!isnan.(A)];
    @test nanpctile!(A, 50) == median(B)
    A = rand(100_000)
    @test nanpctile!(copy(A), 50) == median(A)
    A[rand(1:100_000, 100)] .= NaN; B = A[.!isnan.(A)];
    @test nanpctile!(A, 50) == median(B)

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
    # test for n₋ >= 384
    zi = collect(1:500)
    @test NaNStatistics._nanquantile!(zi, 0.999, 1) == [499.501]
    @test nanquantile(1:500, 0.999) ≈ 499.501

    @test nanquantile(0:10, 0) == 0
    @test nanquantile(0:10, 1/100) ≈ 0.1
    @test nanquantile(0:10, 1.) == 10
    @test nanquantile(0:10, 0.13582) ≈ 1.3582
    @test nanquantile(collect(1:10), 0.5) == 5.5
    @test nanquantile(fill(NaN,10), 0.5) === NaN

    A = rand(100)
    @test nanquantile(A, 0.5) == median(A)
    A[rand(1:100, 20)] .= NaN; B = A[.!isnan.(A)];
    @test nanquantile(A, 0.5) == median(B)
    A = rand(10_000)
    @test nanquantile(A, 0.5) == median(A)
    A[rand(1:10_000, 100)] .= NaN; B = A[.!isnan.(A)];
    @test nanquantile(A, 0.5) == median(B)
    A = rand(100_000)
    @test nanquantile(A, 0.5) == median(A)
    A[rand(1:100_000, 100)] .= NaN; B = A[.!isnan.(A)];
    @test nanquantile(A, 0.5) == median(B)

    A = rand(55,82)
    @test nanquantile(A, 0.5) == median(A)
    @test nanquantile(A, 0.5, dims=1) == median(A, dims=1)
    @test nanquantile(A, 0.5, dims=2) == median(A, dims=2)

    A = rand(10,11,12)
    @test nanquantile(A, 0.5) == median(A)
    @test nanquantile(A, 0.5, dims=1) == median(A, dims=1)
    @test nanquantile(A, 0.5, dims=2) == median(A, dims=2)
    @test nanquantile(A, 0.5, dims=3) == median(A, dims=3)
    @test nanquantile(A, 0.5, dims=(1,2)) == median(A, dims=(1,2))
    @test nanquantile(A, 0.5, dims=(2,3)) == median(A, dims=(2,3))



## ---
