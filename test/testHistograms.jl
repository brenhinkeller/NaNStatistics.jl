N = 10000
nbins = 10

a = rand(N)
h = histcounts(a, 0, 1, nbins)
@test isa(h, Array{Int64,1})
@test length(h) == 10
@test sum(h) == N

@test histcounts(0:99,0,100,10) == fill(10,10)
@test histcounts(0:99.,0.,100.,10) == fill(10,10)
