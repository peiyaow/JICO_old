library(quadprog)
library(nleqslv)

maxmin = function(X, hatb, xi = 0.05){
  G = ncol(hatb)
  n = nrow(X)
  hatS <- t(X) %*% X/n
  #empirical covariance matrix of X
  H <- t(hatb) %*% hatS %*% hatb + xi * diag(G)
  #assume that it is positive definite
  #(use H + xi * I, xi > 0 small, otherwise)
  A <- rbind(rep(1,G),diag(1,G))
  #constraints
  b <- c(1,rep(0,G))
  d <- rep(0,G)
  #linear term is zero
  w <- solve.QP(H,d,t(A),b, meq = 1)
  #quadratic programming solution to
  #argmin(x^T H x) such that Ax >= b and
  #first inequality is an equality
  return(w)
}

# continuum from original paper 
# stone and brooks 1990
continuum = function(X, Y, gamma = 0.5, N = 50){
  n = nrow(X)
  p = ncol(X)
  X.sd = apply(X, 2, sd)
  X = X/X.sd
  S = t(X)%*%X
  s = t(X)%*%Y
  m = min(n-1, p)
  
  e = rev(eigen(S)$values)[1:m]
  V = eigen(S)$vectors[, rev(seq(1,m))]
  theta = seq(0, 1, length.out = N)[-N]
  rho = c()
  for (i in 1:(m-1)){
    rho = c(rho, e[i]+(e[i+1]-e[i])*theta)
  }
  rho = c(rho, e[m])
  D.list = lapply(1:((m-1)*(N-1)+1), function(jj) diag(gamma*rho[jj], m) + diag((1-gamma)*e))
  
  d = t(V)%*%s
  
  Md.list = lapply(1:((m-1)*(N-1)+1), function(jj) solve(D.list[[jj]])%*%d)
  z.list = lapply(1:((m-1)*(N-1)+1), function(jj) Md.list[[jj]]/norm(Md.list[[jj]], "2"))
  c.list = lapply(1:((m-1)*(N-1)+1), function(jj) V%*%z.list[[jj]]) 
  T.vec = sapply(1:((m-1)*(N-1)+1), function(jj) (t(c.list[[jj]])%*%s)^2*(t(c.list[[jj]])%*%S%*%c.list[[jj]])^(gamma-1))
  Z = z.list[[which.max(T.vec)]]
  C = V%*%Z
  A = diag(e)%*%t(V)%*%C
  
  
  for (k in 1:10){
    Md.list = lapply(1:((m-1)*(N-1)+1), function(jj) (solve(D.list[[jj]])-solve(D.list[[jj]])%*%A%*%solve(t(A)%*%solve(D.list[[jj]])%*%A)%*%t(A)%*%solve(D.list[[jj]]))%*%d)
    z.list = lapply(1:((m-1)*(N-1)+1), function(jj) Md.list[[jj]]/norm(Md.list[[jj]], "2"))
    c.list = lapply(1:((m-1)*(N-1)+1), function(jj) V%*%z.list[[jj]])                  
    T.vec = sapply(1:((m-1)*(N-1)+1), function(jj) (t(c.list[[jj]])%*%s)^2*(t(c.list[[jj]])%*%S%*%c.list[[jj]])^(gamma-1))
    Z = cbind(Z, z.list[[which.max(T.vec)]])
    A = cbind(A, diag(e)%*%z.list[[which.max(T.vec)]])
    print(rho[which.max(T.vec)])
    print((t(z.list[[which.max(T.vec)]])%*%diag(e)%*%z.list[[which.max(T.vec)]])[1,1])
    print(max(T.vec))
    if (max(T.vec) < 1e-4){
      break
    }
  }
  
  C = (V%*%Z)[,-(k+1)]
  C= X.sd*C
  return(C)
}

# continuum exact solution
# exact rho solution from yufeng's paper
continuum.exact = function(X, Y, gamma = 0.5){
  print(gamma)
  n = nrow(X)
  p = ncol(X)
  X.mean = apply(X, 2, mean)
  X.sd = apply(X, 2, sd)
  X = sweep(X, 2, X.mean)
  X = sweep(X, 2, X.sd, FUN = "/")
  S = t(X)%*%X
  s = t(X)%*%Y
  m = min(n-1, p)
  
  e = rev(eigen(S)$values)[1:m]
  V = eigen(S)$vectors[, rev(seq(1,m))]
  d = t(V)%*%s
  
  fn = function(rho, gamma, e, V, d){
    D = diag(gamma*rho, m) + diag((1-gamma)*e)
    Md = solve(D)%*%d
    z = Md/norm(Md, "2")
    return(t(z)%*%diag(e)%*%z - rho)
  }
#  print(nleqslv(mean(e), fn, gamma = gamma, e = e, V = V, d = d))
  nleqslv.res = nleqslv(mean(e), fn, gamma = gamma, e = e, V = V, d = d, 
                        method = "Broyden", global = NULL, 
                        control = list(xtol=1e-5, ftol = 1e-5, maxit = 100))
#   print(nleqslv.res$termcd)
  rho = nleqslv.res$x
  if (nleqslv.res$termcd != 1){
    print(paste0("Warning! The value is ", as.character(fn(rho, gamma, e, V, d))))
#     print(fn(rho, gamma, e, V, d))
#     print(nleqslv.res$termcd)
  }
  D = diag(gamma*rho, m) + diag((1-gamma)*e)
  Md = solve(D)%*%d
  Z = Md/norm(Md, "2")
  C = V%*%Z
  A = diag(e)%*%t(V)%*%C
  fn = function(rho, gamma, e, V, d, A){
    D = diag(gamma*rho, m) + diag((1-gamma)*e)
    M = solve(D)-solve(D)%*%A%*%solve(t(A)%*%solve(D)%*%A)%*%t(A)%*%solve(D)
    Md = M%*%d
    z = Md/norm(Md, "2")
    return(t(z)%*%diag(e)%*%z - rho)
  }
#  print(A)
  while (1/kappa(t(A)%*%solve(D)%*%A) > 1e-5){
#    print(nleqslv(rho, fn, gamma = gamma, e = e, V = V, d = d, A = A))
    nleqslv.res = nleqslv(rho, fn, gamma = gamma, e = e, V = V, d = d, A = A, 
                          method = "Broyden", global = NULL, 
                          control = list(xtol=1e-5, ftol = 1e-5, maxit = 100))
    rho = nleqslv.res$x
#     print(nleqslv.res$termcd)
    if (nleqslv.res$termcd != 1){
      print(paste0("Warning! The value is ", as.character(fn(rho, gamma, e, V, d, A))))
#       print(fn(rho, gamma, e, V, d, A))
#       print(nleqslv.res$termcd)
    }
    D = diag(gamma*rho, m) + diag((1-gamma)*e)
    M = solve(D)-solve(D)%*%A%*%solve(t(A)%*%solve(D)%*%A)%*%t(A)%*%solve(D)
    Md = M%*%d
    z = Md/norm(Md, "2")
    Z = cbind(Z, z)
    A = cbind(A, diag(e)%*%z)
#    print(A)
  }
  C = V%*%Z
#  print(C)
  C = C[,-ncol(C)]
  C = C/X.sd
  return(as.matrix(C))
}

# cross validation for continuum.exact
cv.continuum.exact = function(X, Y, gamma = 0.5, nfolds = 10){
  C = continuum.exact(X, Y, gamma)
  n.C = ncol(C)
  # C = as.matrix(C[, rev(seq(n.C))])
  
  if (n.C > 1){
    flds = createFolds(Y, k = nfolds, list = TRUE, returnTrain = FALSE)
    mse.list = list()
    for (k in 1:nfolds){
      X.train = X[unlist(flds[-k]), ]
      X.val = X[unlist(flds[k]), ]
      Y.train = Y[unlist(flds[-k])]
      Y.val = Y[unlist(flds[k])]
      C.train = continuum.exact(X.train, Y.train, gamma)
      n.C.train = ncol(C.train)
      # C.train = as.matrix(C.train[, rev(seq(n.C.train))])
      
      XC.train = X.train%*%C.train
      nC = min(n.C, n.C.train)
      bhat.list = lapply(1:nC, function(c) lm(Y.train~XC.train[,1:c])$coefficients)
      Yhat.val.list = lapply(1:nC, function(c) X.val%*%C.train[,1:c]%*%bhat.list[[c]][-1]+bhat.list[[c]][1])
      mse.list[[k]] = sapply(1:nC, function(c) mean((Yhat.val.list[[c]]-Y.val)^2))
      if (length(mse.list[[k]]) < n.C){
        mse.list[[k]] = c(mse.list[[k]], rep(NA, n.C-length(mse.list[[k]])))
      }
    }
    mse.result = apply(do.call(rbind, mse.list), 2, function(x) mean(x, na.rm = T))
    # print(apply(do.call(rbind, mse.list), 2, function(x) sd(x, na.rm = T)))
    # print(do.call(rbind, mse.list))
    # print(mse.result)
    n.C = which.min(mse.result)
  }
  return(as.matrix(C[, 1:n.C]))
}

# continuum ridge using exact solution
continuum.ridge0 = function(X, Y, lambda, gamma = 0){
  n = nrow(X)
  p = ncol(X)
  X.mean = apply(X, 2, mean)
  X.sd = apply(X, 2, sd)
  X = sweep(X, 2, X.mean)
  X = sweep(X, 2, X.sd, FUN = "/")
  Y.mean = mean(Y)
#  Y.sd = sd(Y)
  Y = sweep(matrix(Y), 2, Y.mean)
#  Y = sweep(matrix(Y), 2, Y.sd, FUN = "/")
  
  S = t(X)%*%X
  s = t(X)%*%Y
  m = min(n-1, p)
  
  e = rev(eigen(S)$values)[1:m]
  V = eigen(S)$vectors[, rev(seq(1,m))]
  d = t(V)%*%s
  
  fn = function(rho, gamma, lambda, e, V, d){
    D = diag(gamma*rho - lambda, m) + diag((1-gamma)*e)
    Md = solve(D)%*%d
    z = Md/norm(Md, "2")
    return(t(z)%*%diag(e)%*%z + lambda - rho)
  }
  #  print(nleqslv(mean(e), fn, gamma = gamma, e = e, V = V, d = d))
  nleqslv.res = nleqslv(e[m]+lambda, fn, gamma = gamma, lambda = lambda, e = e, V = V, d = d, 
                        method = "Broyden", global = NULL, 
                        control = list(xtol=1e-5, ftol = 1e-5, maxit = 100))
  #   print(nleqslv.res$termcd)
  rho = nleqslv.res$x
  print(rho)
  if (nleqslv.res$termcd != 1){
    print(paste0("Warning! The value is ", as.character(fn(rho, gamma, lambda, e, V, d))))
    #     print(fn(rho, gamma, e, V, d))
    #     print(nleqslv.res$termcd)
  }
  D = diag(gamma*rho - lambda, m) + diag((1-gamma)*e)
  Md = solve(D)%*%d
  Z = Md/norm(Md, "2")
  C = V%*%Z
  A = diag(e)%*%t(V)%*%C
  rho0 = rho
  
  fn = function(rho, gamma, lambda, e, V, d, A){
    D = diag(gamma*rho - lambda, m) + diag((1-gamma)*e)
    M = solve(D)-solve(D)%*%A%*%solve(t(A)%*%solve(D)%*%A)%*%t(A)%*%solve(D)
    Md = M%*%d
    z = Md/norm(Md, "2")
    return(t(z)%*%diag(e)%*%z + lambda - rho)
  }
  #  print(A)
  while (1/kappa(t(A)%*%solve(D)%*%A) > 1e-7){
    #    print(nleqslv(rho, fn, gamma = gamma, e = e, V = V, d = d, A = A))
    nleqslv.res = nleqslv(rho0, fn, gamma = gamma, lambda = lambda, e = e, V = V, d = d, A = A, 
                          method = "Broyden", global = NULL, 
                          control = list(xtol=1e-5, ftol = 1e-5, maxit = 100))
    rho = nleqslv.res$x
    print(rho)
    #     print(nleqslv.res$termcd)
    if (nleqslv.res$termcd != 1){
      print(paste0("Warning! The value is ", as.character(fn(rho, gamma, lambda, e, V, d, A))))
    }
    if (abs(rho0-rho) < 0.1){
      break
    }
    rho0 = rho
    D = diag(gamma*rho - lambda, m) + diag((1-gamma)*e)
    M = solve(D)-solve(D)%*%A%*%solve(t(A)%*%solve(D)%*%A)%*%t(A)%*%solve(D)
    Md = M%*%d
    z = Md/norm(Md, "2")
    Z = cbind(Z, z)
    A = cbind(A, diag(e)%*%z)
    #    print(A)
  }
  C = V%*%Z
  #  print(C)
#  C = C[,-ncol(C)]
  C = C/X.sd
  return(as.matrix(C))
}

# main function
continuum.ridge = function(X, Y, lambda, gamma = 0){
  n = nrow(X)
  p = ncol(X)
  X.mean = apply(X, 2, mean)
  X.sd = apply(X, 2, sd)
  X = sweep(X, 2, X.mean)
  X = sweep(X, 2, X.sd, FUN = "/")
  Y.mean = mean(Y)
  Y.sd = sd(Y)
  Y = sweep(matrix(Y), 2, Y.mean)
  Y = sweep(matrix(Y), 2, Y.sd, FUN = "/")
  
  S = t(X)%*%X + lambda*diag(p)
  s = t(X)%*%Y
  m = min(n-1, p)
  
  e = rev(eigen(S)$values)[1:m]
  V = eigen(S)$vectors[, rev(seq(1,m))]
  d = t(V)%*%s
  
  fn = function(rho, gamma, e, V, d){
    D = diag(gamma*rho, m) + diag((1-gamma)*e)
    Md = solve(D)%*%d
    z = Md/norm(Md, "2")
    return(t(z)%*%diag(e)%*%z - rho)
  }
  #  print(nleqslv(mean(e), fn, gamma = gamma, e = e, V = V, d = d))
  nleqslv.res = nleqslv(e[m], fn, gamma = gamma, e = e, V = V, d = d, 
                        method = "Broyden", global = NULL, 
                        control = list(xtol=1e-5, ftol = 1e-5, maxit = 100))
  #   print(nleqslv.res$termcd)
  rho = nleqslv.res$x
  # print(rho)
  if (nleqslv.res$termcd != 1){
    print(paste0("Warning! The value is ", as.character(fn(rho, gamma, e, V, d))))
    #     print(fn(rho, gamma, e, V, d))
    #     print(nleqslv.res$termcd)
  }
  D = diag(gamma*rho, m) + diag((1-gamma)*e)
  Md = solve(D)%*%d
  Z = Md/norm(Md, "2")
  C = V%*%Z
  A = diag(e)%*%t(V)%*%C
  rho0 = rho
  
  fn = function(rho, gamma, e, V, d, A){
    D = diag(gamma*rho, m) + diag((1-gamma)*e)
    M = solve(D)-solve(D)%*%A%*%solve(t(A)%*%solve(D)%*%A)%*%t(A)%*%solve(D)
    Md = M%*%d
    z = Md/norm(Md, "2")
    return(t(z)%*%diag(e)%*%z - rho)
  }
  #  print(A)
  while (1/kappa(t(A)%*%solve(D)%*%A) > 1e-7){
    #    print(nleqslv(rho, fn, gamma = gamma, e = e, V = V, d = d, A = A))
    nleqslv.res = nleqslv(rho0, fn, gamma = gamma, e = e, V = V, d = d, A = A, 
                          method = "Broyden", global = NULL, 
                          control = list(xtol=1e-5, ftol = 1e-5, maxit = 100))
    rho = nleqslv.res$x
    # print(rho)
    #     print(nleqslv.res$termcd)
    if (nleqslv.res$termcd != 1){
      print(paste0("Warning! The value is ", as.character(fn(rho, gamma, e, V, d, A))))
    }
    if (abs(rho0-rho) < 0.1){
      break
    }
    rho0 = rho
    D = diag(gamma*rho, m) + diag((1-gamma)*e)
    M = solve(D)-solve(D)%*%A%*%solve(t(A)%*%solve(D)%*%A)%*%t(A)%*%solve(D)
    Md = M%*%d
    z = Md/norm(Md, "2")
    Z = cbind(Z, z)
    A = cbind(A, diag(e)%*%z)
    #    print(A)
  }
  C = V%*%Z
  C = C[,1:min(m, ncol(C))]
#  C = C/X.sd*Y.sd
  return(as.matrix(C))
}

# do not scale Y by Y.sd
continuum.ridge.sd = function(X, Y, lambda, gamma = 0){
  n = nrow(X)
  p = ncol(X)
  X.mean = apply(X, 2, mean)
  X.sd = apply(X, 2, sd)
  X = sweep(X, 2, X.mean)
  X = sweep(X, 2, X.sd, FUN = "/")
  Y.mean = mean(Y)
#  Y.sd = sd(Y)
  Y = sweep(matrix(Y), 2, Y.mean)
#  Y = sweep(matrix(Y), 2, Y.sd, FUN = "/")
  
  S = t(X)%*%X + lambda*diag(p)
  s = t(X)%*%Y
  m = min(n-1, p)
  
  e = rev(eigen(S)$values)[1:m]
  V = eigen(S)$vectors[, rev(seq(1,m))]
  d = t(V)%*%s
  
  fn = function(rho, gamma, e, V, d){
    D = diag(gamma*rho, m) + diag((1-gamma)*e)
    Md = solve(D)%*%d
    z = Md/norm(Md, "2")
    return(t(z)%*%diag(e)%*%z - rho)
  }
  #  print(nleqslv(mean(e), fn, gamma = gamma, e = e, V = V, d = d))
  nleqslv.res = nleqslv(e[m], fn, gamma = gamma, e = e, V = V, d = d, 
                        method = "Broyden", global = NULL, 
                        control = list(xtol=1e-5, ftol = 1e-5, maxit = 100))
  #   print(nleqslv.res$termcd)
  rho = nleqslv.res$x
  # print(rho)
  if (nleqslv.res$termcd != 1){
    print(paste0("Warning! The value is ", as.character(fn(rho, gamma, e, V, d))))
    #     print(fn(rho, gamma, e, V, d))
    #     print(nleqslv.res$termcd)
  }
  D = diag(gamma*rho, m) + diag((1-gamma)*e)
  Md = solve(D)%*%d
  Z = Md/norm(Md, "2")
  C = V%*%Z
  A = diag(e)%*%t(V)%*%C
  rho0 = rho
  
  fn = function(rho, gamma, e, V, d, A){
    D = diag(gamma*rho, m) + diag((1-gamma)*e)
    M = solve(D)-solve(D)%*%A%*%solve(t(A)%*%solve(D)%*%A)%*%t(A)%*%solve(D)
    Md = M%*%d
    z = Md/norm(Md, "2")
    return(t(z)%*%diag(e)%*%z - rho)
  }
  #  print(A)
  while (1/kappa(t(A)%*%solve(D)%*%A) > 1e-7){
    #    print(nleqslv(rho, fn, gamma = gamma, e = e, V = V, d = d, A = A))
    nleqslv.res = nleqslv(rho0, fn, gamma = gamma, e = e, V = V, d = d, A = A, 
                          method = "Broyden", global = NULL, 
                          control = list(xtol=1e-5, ftol = 1e-5, maxit = 100))
    rho = nleqslv.res$x
    # print(rho)
    #     print(nleqslv.res$termcd)
    if (nleqslv.res$termcd != 1){
      print(paste0("Warning! The value is ", as.character(fn(rho, gamma, e, V, d, A))))
    }
    if (abs(rho0-rho) < 0.1){
      break
    }
    rho0 = rho
    D = diag(gamma*rho, m) + diag((1-gamma)*e)
    M = solve(D)-solve(D)%*%A%*%solve(t(A)%*%solve(D)%*%A)%*%t(A)%*%solve(D)
    Md = M%*%d
    z = Md/norm(Md, "2")
    Z = cbind(Z, z)
    A = cbind(A, diag(e)%*%z)
    #    print(A)
  }
  C = V%*%Z[,1:m]
  #  print(C)
  #  C = C[,-ncol(C)]
#  C = C/X.sd
  C = C[,1:min(m, ncol(C))]
  return(as.matrix(C))
}

# cc output from continuum.ridge
C2beta = function(X, Y, cc, lambda){
  X.mean = apply(X, 2, mean)
  X.sd = apply(X, 2, sd)
  X = sweep(X, 2, X.mean)
  X = sweep(X, 2, X.sd, FUN = "/")
  Y.mean = mean(Y)
  Y.sd = sd(Y)
  Y = sweep(matrix(Y), 2, Y.mean)
  Y = sweep(matrix(Y), 2, Y.sd, FUN = "/")
  
  cc = matrix(cc, nrow = ncol(X))
  beta.C = t(t(solve(t(cc)%*%t(X)%*%X%*%cc + lambda*diag(ncol(cc)))%*%(t(cc)%*%t(X)%*%(Y)))%*%t(cc))
  
  beta.C = beta.C/X.sd*Y.sd
  intercept = Y.mean - X.mean%*%beta.C

  return(list(intercept = intercept, beta = beta.C, coef = matrix(c(intercept, beta.C))))
}

# fix gamma tune lambda and omega
cv.continuum.ridge = function(X, Y, lambda, gamma = 0.5, nfolds = 10){
  n = nrow(X)
  p = ncol(X)
  n.lambda = length(lambda)
  flds = createFolds(Y, k = nfolds, list = TRUE, returnTrain = FALSE)
  m = min(n-1, p)
  mse.list = list()
  for (k in 1:nfolds){
    X.train = X[unlist(flds[-k]), ]
    X.val = X[unlist(flds[k]), ]
    Y.train = Y[unlist(flds[-k])]
    Y.val = Y[unlist(flds[k])]
    C.list = lapply(lambda, function(lam) continuum.ridge(X.train, Y.train, lam, gamma = gamma)) # length lambda
    n.C.vec = sapply(C.list, ncol)
    beta.list = lapply(1:n.lambda, function(ix.lam) lapply(1:n.C.vec[ix.lam], function(ix.C) C2beta(X.train, Y.train, C.list[[ix.lam]][,1:ix.C], lambda[ix.lam])$coef))
    # print(beta.list)
    mse.list[[k]] = lapply(1:n.lambda, function(ix.lam) lapply(1:n.C.vec[ix.lam], function(ix.C) mean((Y.val - cbind(1, X.val)%*%beta.list[[ix.lam]][[ix.C]])^2)))
    mse.list[[k]] = sapply(1:n.lambda, function(ix.lam) c(unlist(mse.list[[k]][[ix.lam]]), rep(NA, (m - n.C.vec)[ix.lam])))
    }
  MSE = Reduce('+', mse.list)/nfolds
  print(MSE)
  n.C = nrow(MSE) 
  min.id = which.min(MSE)
  n.C.selected = (min.id - 1)%%n.C + 1
  min.id.lam = ceiling(min.id/n.C)
  lam.selected =lambda[min.id.lam]
  C.selected = continuum.ridge(X, Y, lambda = lam.selected, gamma = gamma)[,1:n.C.selected]
  
  return(list(C = C.selected, n.C = n.C.selected, lam = lam.selected, MSE = MSE))
}








