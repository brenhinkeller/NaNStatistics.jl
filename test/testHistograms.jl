N = 10000

a = rand(N)
h = histcounts(a, 0, 1, 10)
@test h == histcounts(a, 0:0.1:1)
@test isa(h, Array{Int64,1})
@test length(h) == 10
@test sum(h) == N

@test histcounts(0:99,0,100,10) == fill(10,10)
@test histcounts(0:99.,0.,100.,10) == fill(10,10)
