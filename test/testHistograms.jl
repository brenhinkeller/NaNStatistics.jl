N = 10000

a = rand(N)
h = histcounts(a, 0:0.1:1)
@test h == histcounts(a, 0, 1, 10)
@test isa(h, Array{Int64,1})
@test length(h) == 10
@test sum(h) == N

@test histcounts(1:100,0:10:100) == fill(10,10)
@test histcounts(1:100,0,100,10) == fill(10,10)
@test histcounts(1:100.,0.,100.,10) == fill(10,10)

# Test results and warnings when N is too small to hold results
w = "length(N) < nbins; any bins beyond length(N) will not be filled"
@test (@test_logs (:warn, w) histcounts!(zeros(5), 1:100 ,0:10:100)) == fill(10,5)
@test (@test_logs (:warn, w) histcounts!(zeros(5), 1:100. ,0:10:100)) == fill(10,5)
