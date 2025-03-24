# Test pair-wise correlation and covariance functions
x, y = rand(1000), rand(1000)

@test nancov(x,y) ≈ cov(x,y)
@test nancov(x,y, corrected=false) ≈ cov(x,y, corrected=false)
@test nancovem(x,y) ≈ cov(x,y)./length(x)
@test nancovem(x,y, corrected=false) ≈ cov(x,y, corrected=false)./length(x)
@test nancor(x,y) ≈ cor(x,y)


# Test correlation and covariance functions as applied to matrices
X = rand(100,10)

@test nancov(X) ≈ cov(X)
@test nancov(X, dims=1) ≈ cov(X, dims=1)
@test nancov(X, dims=2) ≈ cov(X, dims=2)
@test nancov(X, corrected=false) ≈ cov(X, corrected=false)

@test nancovem(X) ≈ cov(X)./size(X,1)
@test nancovem(X, dims=1) ≈ cov(X, dims=1)./size(X,1)
@test nancovem(X, dims=2) ≈ cov(X, dims=2)./size(X,2)
@test nancovem(X, corrected=false) ≈ cov(X, corrected=false)./size(X,1)

@test nancor(X) ≈ cor(X)
@test nancor(X, dims=1) ≈ cor(X, dims=1)
@test nancor(X, dims=2) ≈ cor(X, dims=2)

# Test NaN-ignoring
xn = [1,2,3,NaN]
xnn = [1,2,3]
@test nancov(xn,xn) ≈ cov(xnn,xnn) ≈ 1
@test nancor(xn,xn) ≈ cor(xnn,xnn) ≈ 1

x[rand(1:length(x), 100)] .= NaN
y[rand(1:length(y), 100)] .= NaN
t = .!(isnan.(x) .| isnan.(y))
@test nancov(x,y) ≈ cov(x[t],y[t])
@test nancor(x,y) ≈ cor(x[t],y[t])
