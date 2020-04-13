library(MASS)

# homogeneous X with diag
generateX = function(n, p, s = 1){
  X = mvrnorm(n, rep(0, p), diag(p))
  X.mean = apply(X, 2, mean)
  X.sd = apply(X, 2, sd)
  X = sweep(sweep(X, 2, X.mean), 2, X.sd, "/")
  X = X*s
  return(X)
}

AR1 = function(rrho, K){
  R = matrix(, nrow = K, ncol = K)
  for (i in 1:K){
    for (j in 1:K){
      R[i,j] = rrho^(abs(i-j))
    }
  }
  return(R)
}

generateX_AR = function(n, p, r, s = 1){
  C = AR1(r, p)
  X = mvrnorm(n, rep(0, p), C)
  X.mean = apply(X, 2, mean)
  X.sd = apply(X, 2, sd)
  X = sweep(sweep(X, 2, X.mean), 2, X.sd, "/")
  X = X*s
  return(X)
}



