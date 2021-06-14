# Test pair-wise correlation and covariance functions
x, y = rand(1000), rand(1000)

@test nancov(x,y) ≈ cov(x,y)
@test nancov(x,y, corrected=false) ≈ cov(x,y, corrected=false)
@test nancor(x,y) ≈ cor(x,y)


# Test correlation and covariance functions as applied to matrices
X = rand(100,10)

@test nancov(X) ≈ cov(X)
@test nancov(X, dims=1) ≈ cov(X, dims=1)
@test nancov(X, dims=2) ≈ cov(X, dims=2)
@test nancov(X, corrected=false) ≈ cov(X, corrected=false)

@test nancor(X) ≈ cor(X)
@test nancor(X, dims=1) ≈ cor(X, dims=1)
@test nancor(X, dims=2) ≈ cor(X, dims=2)

# Test NaN-ignoring
xn = [1,2,3,NaN]
xnn = [1,2,3]
@test nancov(xn,xn) ≈ cov(xnn,xnn) ≈ 1
@test nancor(xn,xn) ≈ cor(xnn,xnn) ≈ 1
