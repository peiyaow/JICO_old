library(quadprog)
library(nleqslv)
library(caret)

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
continuum = function(X, Y, gam = 0.5, N = 50){
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
  D.list = lapply(1:((m-1)*(N-1)+1), function(jj) diag(gam*rho[jj], m) + diag((1-gam)*e))
  
  d = t(V)%*%s
  
  Md.list = lapply(1:((m-1)*(N-1)+1), function(jj) solve(D.list[[jj]])%*%d)
  z.list = lapply(1:((m-1)*(N-1)+1), function(jj) Md.list[[jj]]/norm(Md.list[[jj]], "2"))
  c.list = lapply(1:((m-1)*(N-1)+1), function(jj) V%*%z.list[[jj]]) 
  T.vec = sapply(1:((m-1)*(N-1)+1), function(jj) (t(c.list[[jj]])%*%s)^2*(t(c.list[[jj]])%*%S%*%c.list[[jj]])^(gam-1))
  Z = z.list[[which.max(T.vec)]]
  C = V%*%Z
  A = diag(e)%*%t(V)%*%C
  
  
  for (k in 1:10){
    Md.list = lapply(1:((m-1)*(N-1)+1), function(jj) (solve(D.list[[jj]])-solve(D.list[[jj]])%*%A%*%solve(t(A)%*%solve(D.list[[jj]])%*%A)%*%t(A)%*%solve(D.list[[jj]]))%*%d)
    z.list = lapply(1:((m-1)*(N-1)+1), function(jj) Md.list[[jj]]/norm(Md.list[[jj]], "2"))
    c.list = lapply(1:((m-1)*(N-1)+1), function(jj) V%*%z.list[[jj]])                  
    T.vec = sapply(1:((m-1)*(N-1)+1), function(jj) (t(c.list[[jj]])%*%s)^2*(t(c.list[[jj]])%*%S%*%c.list[[jj]])^(gam-1))
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
continuum.exact = function(X, Y, gam = 0.5){
  print(gam)
  n = nrow(X)
  p = ncol(X)
  X.mean = apply(X, 2, mean)
  X.sd = apply(X, 2, sd)
  X = sweep(X, 2, X.mean)
  X = sweep(X, 2, X.sd, FUN = "/")
  S = t(X)%*%X
  s = t(X)%*%Y
  m = min(n-1, p)
  
  #  e = rev(eigen(S)$values)[1:m]
  e = eigen(S)$values[1:m]
  #  V = eigen(S)$vectors[, rev(seq(1,m))]
  V = eigen(S)$vectors[, seq(1,m)]
  d = t(V)%*%s
  
  fn = function(rho, gam, e, V, d){
    D = diag(gam*rho, m) + diag((1-gam)*e)
    Md = solve(D)%*%d
    z = Md/norm(Md, "2")
    return(t(z)%*%diag(e)%*%z - rho)
  }
#  print(nleqslv(mean(e), fn, gam = gam, e = e, V = V, d = d))
  nleqslv.res = nleqslv(mean(e), fn, gam = gam, e = e, V = V, d = d, 
                        method = "Broyden", global = NULL, 
                        control = list(xtol=1e-5, ftol = 1e-5, maxit = 100))
#   print(nleqslv.res$termcd)
  rho = nleqslv.res$x
  if (nleqslv.res$termcd != 1){
    print(paste0("Warning! The value is ", as.character(fn(rho, gam, e, V, d))))
#     print(fn(rho, gam, e, V, d))
#     print(nleqslv.res$termcd)
  }
  D = diag(gam*rho, m) + diag((1-gam)*e)
  Md = solve(D)%*%d
  Z = Md/norm(Md, "2")
  C = V%*%Z
  A = diag(e)%*%t(V)%*%C
  fn = function(rho, gam, e, V, d, A){
    D = diag(gam*rho, m) + diag((1-gam)*e)
    M = solve(D)-solve(D)%*%A%*%solve(t(A)%*%solve(D)%*%A)%*%t(A)%*%solve(D)
    Md = M%*%d
    z = Md/norm(Md, "2")
    return(t(z)%*%diag(e)%*%z - rho)
  }
#  print(A)
  while (1/kappa(t(A)%*%solve(D)%*%A) > 1e-5){
#    print(nleqslv(rho, fn, gam = gam, e = e, V = V, d = d, A = A))
    nleqslv.res = nleqslv(rho, fn, gam = gam, e = e, V = V, d = d, A = A, 
                          method = "Broyden", global = NULL, 
                          control = list(xtol=1e-5, ftol = 1e-5, maxit = 100))
    rho = nleqslv.res$x
#     print(nleqslv.res$termcd)
    if (nleqslv.res$termcd != 1){
      print(paste0("Warning! The value is ", as.character(fn(rho, gam, e, V, d, A))))
#       print(fn(rho, gam, e, V, d, A))
#       print(nleqslv.res$termcd)
    }
    D = diag(gam*rho, m) + diag((1-gam)*e)
    M = solve(D)-solve(D)%*%A%*%solve(t(A)%*%solve(D)%*%A)%*%t(A)%*%solve(D)
    Md = M%*%d
    z = Md/norm(Md, "2")
    Z = cbind(Z, z)
    A = cbind(A, diag(e)%*%z)
#    print(A)
  }
  C = V%*%Z
#  print(C)
  C = C[,1:min(ncol(C), m-1)]
  C = C/X.sd
  return(as.matrix(C))
}

# cross validation for continuum.exact, tune omega
cv.continuum.exact = function(X, Y, gam = 0.5, nfolds = 10){
  C = continuum.exact(X, Y, gam)
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
      C.train = continuum.exact(X.train, Y.train, gam)
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

continuum.ridge = function(X, Y, lambda, gam){
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
  
  S = t(X)%*%X
  s = t(X)%*%Y
#  m = 10
  m = min(n-1, p)

  e = eigen(S)$values[1:m]
  E = diag(e)
  V = eigen(S)$vectors[, 1:m]
  d = t(V)%*%s
  
  tau = 1
  fn = function(rho){
#    A = diag(tau^2*rho^(gam-2)*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*rho^(gam-2)*tau^2*E
    A = diag(tau^2*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*tau^2*E
    M = solve(A)
#    q = tau*rho^(gam-1)*d
    q = tau*rho*d
    Mq = M%*%q
    z = Mq/norm(Mq, "2")
    return(t(z)%*%E%*%z + n*lambda - rho)
  }
  
  nleqslv.res = nleqslv((e[1]+e[m])/2+n*lambda, fn, method = "Newton", control = list(maxit = 50))
  rho = nleqslv.res$x
  if (nleqslv.res$termcd != 1){
    print(paste0("Warning! The value is ", as.character(fn(rho))))
    #     print(fn(rho, gam, e, V, d))
    #     print(nleqslv.res$termcd)
  }
  
#  A = diag(tau^2*rho^(gam-2)*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*rho^(gam-2)*tau^2*E
  A = diag(tau^2*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*tau^2*E
  M = solve(A)
#  q = tau*rho^(gam-1)*d
  q = tau*rho*d
  Mq = M%*%q
  Z = Mq/norm(Mq, "2")
  B = E%*%Z
  rho0 = rho
#  tau = (t(Z)%*%d)[1,1]
  
  SOLVE = function(x){
    if (sum(dim(x))){
      return(solve(x))
    }else{
      return(x)
    }
  }
    
  fn = function(rho){
#    A = diag(tau^2*rho^(gam-2)*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*rho^(gam-2)*tau^2*E
    A = diag(tau^2*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*tau^2*E
    M = solve(A)- solve(A)%*%B%*%SOLVE(t(B)%*%solve(A)%*%B)%*%t(B)%*%solve(A)
#    q = tau*rho^(gam-1)*d
    q = tau*rho*d
    Mq = M%*%q
    z = Mq/norm(Mq, "2")
    return(t(z)%*%E%*%z + n*lambda - rho)
  }
  
#  while (1/kappa(t(B)%*%solve(A)%*%B) > 1e-17){
  while (rcond(t(B)%*%solve(A)%*%B) > 1e-10){
#    print(1)
#    print(rcond(t(B)%*%solve(A)%*%B))
    nleqslv.res = nleqslv(rho0, fn, method = "Newton", control = list(maxit = 50))
    rho = nleqslv.res$x
    if (nleqslv.res$termcd != 1){
      print(paste0("Warning! The value is ", as.character(fn(rho))))
#      break
    }
#    A = diag(tau^2*rho^(gam-2)*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*rho^(gam-2)*tau^2*E
    A = diag(tau^2*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*tau^2*E
#    print(2)
#    print(rcond(t(B)%*%solve(A)%*%B))
#    if (1/kappa(t(B)%*%solve(A)%*%B) < 1e-17){
    if (rcond(t(B)%*%solve(A)%*%B) < 1e-10){
      break
    }
    M = solve(A)- solve(A)%*%B%*%SOLVE(t(B)%*%solve(A)%*%B)%*%t(B)%*%solve(A)
#    q = tau*rho^(gam-1)*d
    q = tau*rho*d
    Mq = M%*%q
    z = Mq/norm(Mq, "2")
    Z = cbind(Z, z)
    B = E%*%Z
    rho0 = rho
#    tau = (t(z)%*%d)[1,1]
    if (ncol(Z) == m-1){
      break
    }
  }
  C = V%*%Z
  #  print(C)
  C = C[,1:min(ncol(C), m-1)]
  return(as.matrix(C))
}

# cc output from continuum.ridge
C2beta = function(X, Y, cc, lambda){
  n = nrow(X)
  X.mean = apply(X, 2, mean)
  X.sd = apply(X, 2, sd)
  X = sweep(X, 2, X.mean)
  X = sweep(X, 2, X.sd, FUN = "/")
  Y.mean = mean(Y)
  Y.sd = sd(Y)
  Y = sweep(matrix(Y), 2, Y.mean)
  Y = sweep(matrix(Y), 2, Y.sd, FUN = "/")
  S = t(X)%*%X
  s = t(X)%*%Y
  
  cc = matrix(cc, nrow = ncol(X))
  om = ncol(cc) #omega
  SOLVE = function(x){
    if (sum(dim(x))){
      return(solve(x))
    }else{
      return(x)
    }
  }
  beta.C = cc%*%SOLVE(t(cc)%*%S%*%cc + n*lambda*diag(om))%*%t(cc)%*%s
  beta.C = beta.C/X.sd*Y.sd
  intercept = Y.mean - X.mean%*%beta.C

  return(list(intercept = intercept, beta = beta.C, coef = matrix(c(intercept, beta.C))))
}

# fix gam tune lambda and omega
cv.continuum.ridge = function(X, Y, lambda, gam = 1, nfolds = 10){
  n = nrow(X)
  p = ncol(X)
  n.lambda = length(lambda)
  flds = createFolds(Y, k = nfolds, list = TRUE, returnTrain = FALSE)
  m = min(n-1, p)
#  m = 10
  mse.list = list()
  for (k in 1:nfolds){
    X.train = X[unlist(flds[-k]), ]
    X.val = X[unlist(flds[k]), ]
    Y.train = Y[unlist(flds[-k])]
    Y.val = Y[unlist(flds[k])]
    C.list = lapply(lambda, function(lam) continuum.ridge(X.train, Y.train, lam, gam = gam)) # length lambda
    n.C.vec = sapply(C.list, ncol)
    print(n.C.vec)
#    beta.list = lapply(1:n.lambda, function(ix.lam) lapply(1:n.C.vec[ix.lam], function(ix.C) C2beta(X.train, Y.train, C.list[[ix.lam]][,1:ix.C], lambda[ix.lam])$coef))
    beta.list = lapply(1:n.lambda, function(ix.lam) lapply(0:n.C.vec[ix.lam], function(ix.C) C2beta(X.train, Y.train, C.list[[ix.lam]][,0:ix.C], lambda[ix.lam])$coef))
#    mse.list[[k]] = lapply(1:n.lambda, function(ix.lam) lapply(1:n.C.vec[ix.lam], function(ix.C) mean((Y.val - cbind(1, X.val)%*%beta.list[[ix.lam]][[ix.C]])^2)))
    mse.list[[k]] = lapply(1:n.lambda, function(ix.lam) lapply(0:n.C.vec[ix.lam], function(ix.C) mean((Y.val - cbind(1, X.val)%*%beta.list[[ix.lam]][[ix.C+1]])^2)))
#    mse.list[[k]] = sapply(1:n.lambda, function(ix.lam) c(unlist(mse.list[[k]][[ix.lam]]), rep(NA, (m - n.C.vec)[ix.lam])))
    mse.list[[k]] = sapply(1:n.lambda, function(ix.lam) c(unlist(mse.list[[k]][[ix.lam]]), rep(NA, (m+1 - n.C.vec)[ix.lam])))
    }
  MSE = Reduce('+', mse.list)/nfolds
#  print(MSE)
  n.C = nrow(MSE) 
  min.id = which.min(MSE)
#  n.C.selected = (min.id - 1)%%n.C + 1
  n.C.selected = (min.id - 1)%%n.C
  min.id.lam = ceiling(min.id/n.C)
  lam.selected = lambda[min.id.lam]
  C = continuum.ridge(X, Y, lambda = lam.selected, gam = gam)
#  C.selected = as.matrix(C[,1:n.C.selected])
  C.selected = as.matrix(C[,0:n.C.selected])
  
  return(list(C = C.selected, C0 = C, n.C = n.C.selected, lam = lam.selected, MSE = MSE))
}

continuum.ridge.res = function(X, Y, C, lambda, gam){
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
  
  S = t(X)%*%X
  s = t(X)%*%Y
  m = min(n-1, p)
  
  e = eigen(S)$values[1:m]
  E = diag(e)
  V = eigen(S)$vectors[, seq(1,m)]
  d = t(V)%*%s
  B = E%*%t(V)%*%C
  
  SOLVE = function(x){
    if (sum(dim(x))){
      return(solve(x))
    }else{
      return(x)
    }
  }
  
  tau = 1
  fn = function(rho){
    #    A = diag(tau^2*rho^(gam-2)*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*rho^(gam-2)*tau^2*E
    A = diag(tau^2*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*tau^2*E
    M = solve(A)- solve(A)%*%B%*%SOLVE(t(B)%*%solve(A)%*%B)%*%t(B)%*%solve(A)
#    q = tau*rho^(gam-1)*d
    q = tau*rho*d
    Mq = M%*%q
    z = Mq/norm(Mq, "2")
    return(t(z)%*%E%*%z + n*lambda - rho)
  }
  
  nleqslv.res = nleqslv((e[1]+e[m])/2+n*lambda, fn)
  rho = nleqslv.res$x
  if (nleqslv.res$termcd != 1){
    print(paste0("Warning! The value is ", as.character(fn(rho))))
  }
  
  #    A = diag(tau^2*rho^(gam-2)*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*rho^(gam-2)*tau^2*E
  A = diag(tau^2*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*tau^2*E
  M = solve(A)- solve(A)%*%B%*%SOLVE(t(B)%*%solve(A)%*%B)%*%t(B)%*%solve(A)
#  q = tau*rho^(gam-1)*d
  q = tau*rho*d
  Mq = M%*%q
  Z = Mq/norm(Mq, "2")
  B = E%*%cbind(t(V)%*%C, Z)
  rho0 = rho
  
  while (rcond(t(B)%*%solve(A)%*%B) > 1e-10){
    nleqslv.res = nleqslv(rho0, fn)
    rho = nleqslv.res$x
    if (nleqslv.res$termcd != 1){
      print(paste0("Warning! The value is ", as.character(fn(rho))))
      break
    }
#    A = diag(tau^2*rho^(gam-2)*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*rho^(gam-2)*tau^2*E
    A = diag(tau^2*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*tau^2*E
    if (rcond(t(B)%*%solve(A)%*%B) < 1e-10){
      break
    }
    M = solve(A)- solve(A)%*%B%*%SOLVE(t(B)%*%solve(A)%*%B)%*%t(B)%*%solve(A)
#    q = tau*rho^(gam-1)*d
    q = tau*rho*d
    Mq = M%*%q
    z = Mq/norm(Mq, "2")
    Z = cbind(Z, z)
    B = E%*%cbind(t(V)%*%C, Z)
    rho0 = rho
    if (ncol(Z) == m-1){
      break
    }
  }
  C = V%*%Z
  #  print(C)
  C = C[,1:min(ncol(C), m-1)]
  return(as.matrix(C))
}

# fix gam tune lambda and omega
cv.continuum.ridge.res = function(X, Y, cc, lambda, gam = 1, nfolds = 10){
  n = nrow(X)
  p = ncol(X)
  n.lambda = length(lambda)
  flds = createFolds(Y, k = nfolds, list = TRUE, returnTrain = FALSE)
  m = min(n-1, p)
  #  m = 10
  mse.list = list()
  for (k in 1:nfolds){
    X.train = X[unlist(flds[-k]), ]
    X.val = X[unlist(flds[k]), ]
    Y.train = Y[unlist(flds[-k])]
    Y.val = Y[unlist(flds[k])]
    C.list = lapply(lambda, function(lam) continuum.ridge.res(X.train, Y.train, cc, lam, gam = gam)) # length lambda
    n.C.vec = sapply(C.list, ncol)
    print(n.C.vec)
    #    beta.list = lapply(1:n.lambda, function(ix.lam) lapply(1:n.C.vec[ix.lam], function(ix.C) C2beta(X.train, Y.train, C.list[[ix.lam]][,1:ix.C], lambda[ix.lam])$coef))
    beta.list = lapply(1:n.lambda, function(ix.lam) lapply(0:n.C.vec[ix.lam], function(ix.C) C2beta(X.train, Y.train, C.list[[ix.lam]][,0:ix.C], lambda[ix.lam])$coef))
    #    mse.list[[k]] = lapply(1:n.lambda, function(ix.lam) lapply(1:n.C.vec[ix.lam], function(ix.C) mean((Y.val - cbind(1, X.val)%*%beta.list[[ix.lam]][[ix.C]])^2)))
    mse.list[[k]] = lapply(1:n.lambda, function(ix.lam) lapply(0:n.C.vec[ix.lam], function(ix.C) mean((Y.val - cbind(1, X.val)%*%beta.list[[ix.lam]][[ix.C+1]])^2)))
    #    mse.list[[k]] = sapply(1:n.lambda, function(ix.lam) c(unlist(mse.list[[k]][[ix.lam]]), rep(NA, (m - n.C.vec)[ix.lam])))
    mse.list[[k]] = sapply(1:n.lambda, function(ix.lam) c(unlist(mse.list[[k]][[ix.lam]]), rep(NA, (m+1 - n.C.vec)[ix.lam])))
  }
  MSE = Reduce('+', mse.list)/nfolds
  #  print(MSE)
  n.C = nrow(MSE) 
  min.id = which.min(MSE)
  #  n.C.selected = (min.id - 1)%%n.C + 1
  n.C.selected = (min.id - 1)%%n.C
  min.id.lam = ceiling(min.id/n.C)
  lam.selected = lambda[min.id.lam]
  C = continuum.ridge.res(X, Y, cc, lambda = lam.selected, gam = gam)
  #  C.selected = as.matrix(C[,1:n.C.selected])
  C.selected = as.matrix(C[,0:n.C.selected])
  
  return(list(C = C.selected, C0 = C, n.C = n.C.selected, lam = lam.selected, MSE = MSE))
}


continuum.ridge.fix = function(X, Y, lambda, gam, om){
  #om: number of columns
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
  
  S = t(X)%*%X
  s = t(X)%*%Y
  m = min(n-1, p)
  
  e = eigen(S)$values[1:m]
  E = diag(e)
  V = eigen(S)$vectors[, 1:m]
  d = t(V)%*%s
  
  tau = 1
  fn = function(rho){
    #    A = diag(tau^2*rho^(gam-2)*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*rho^(gam-2)*tau^2*E
    A = diag(tau^2*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*tau^2*E
    M = solve(A)
    #    q = tau*rho^(gam-1)*d
    q = tau*rho*d
    Mq = M%*%q
    z = Mq/norm(Mq, "2")
    return(t(z)%*%E%*%z + n*lambda - rho)
  }
  
  nleqslv.res = nleqslv((e[1]+e[m])/2+n*lambda, fn, method = "Newton", control = list(maxit = 50))
  rho = nleqslv.res$x
  if (nleqslv.res$termcd != 1){
    print(paste0("Warning! The value is ", as.character(fn(rho))))
    #     print(fn(rho, gam, e, V, d))
    #     print(nleqslv.res$termcd)
  }
  
  #  A = diag(tau^2*rho^(gam-2)*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*rho^(gam-2)*tau^2*E
  A = diag(tau^2*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*tau^2*E
  M = solve(A)
  #  q = tau*rho^(gam-1)*d
  q = tau*rho*d
  Mq = M%*%q
  Z = Mq/norm(Mq, "2")
  B = E%*%Z
  rho0 = rho
  #  tau = (t(Z)%*%d)[1,1]
  
  SOLVE = function(x){
    if (sum(dim(x))){
      return(solve(x))
    }else{
      return(x)
    }
  }
  
  fn = function(rho){
    #    A = diag(tau^2*rho^(gam-2)*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*rho^(gam-2)*tau^2*E
    A = diag(tau^2*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*tau^2*E
    M = solve(A)- solve(A)%*%B%*%SOLVE(t(B)%*%solve(A)%*%B)%*%t(B)%*%solve(A)
    #    q = tau*rho^(gam-1)*d
    q = tau*rho*d
    Mq = M%*%q
    z = Mq/norm(Mq, "2")
    return(t(z)%*%E%*%z + n*lambda - rho)
  }
  
  while (ncol(Z) < om & rcond(t(B)%*%solve(A)%*%B) > 1e-10){
    #    print(1)
    #    print(rcond(t(B)%*%solve(A)%*%B))
    nleqslv.res = nleqslv(rho0, fn, method = "Newton", control = list(maxit = 50))
    rho = nleqslv.res$x
    if (nleqslv.res$termcd != 1){
      print(paste0("Warning! The value is ", as.character(fn(rho))))
      #      break
    }
    #    A = diag(tau^2*rho^(gam-2)*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*rho^(gam-2)*tau^2*E
    A = diag(tau^2*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*tau^2*E
    #    print(2)
    #    print(rcond(t(B)%*%solve(A)%*%B))
    if (rcond(t(B)%*%solve(A)%*%B) < 1e-10){
      break
    }
    M = solve(A)- solve(A)%*%B%*%SOLVE(t(B)%*%solve(A)%*%B)%*%t(B)%*%solve(A)
    #    q = tau*rho^(gam-1)*d
    q = tau*rho*d
    Mq = M%*%q
    z = Mq/norm(Mq, "2")
    Z = cbind(Z, z)
    B = E%*%Z
    rho0 = rho
    #    tau = (t(z)%*%d)[1,1]
  }
  C = V%*%Z
  C = C[,0:min(ncol(C), om)]
  return(as.matrix(C))
}

continuum.ridge.res.fix = function(X, Y, C, lambda, gam, om){
  #om: number of columns
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
  
  S = t(X)%*%X
  s = t(X)%*%Y
  m = min(n-1, p)
  
  e = eigen(S)$values[1:m]
  E = diag(e)
  V = eigen(S)$vectors[, seq(1,m)]
  d = t(V)%*%s
  B = E%*%t(V)%*%C
  
  SOLVE = function(x){
    if (sum(dim(x))){
      return(solve(x))
    }else{
      return(x)
    }
  }
  
  tau = 1
  fn = function(rho){
    #    A = diag(tau^2*rho^(gam-2)*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*rho^(gam-2)*tau^2*E
    A = diag(tau^2*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*tau^2*E
    M = solve(A)- solve(A)%*%B%*%SOLVE(t(B)%*%solve(A)%*%B)%*%t(B)%*%solve(A)
    #    q = tau*rho^(gam-1)*d
    q = tau*rho*d
    Mq = M%*%q
    z = Mq/norm(Mq, "2")
    return(t(z)%*%E%*%z + n*lambda - rho)
  }
  
  nleqslv.res = nleqslv((e[1]+e[m])/2+n*lambda, fn)
  rho = nleqslv.res$x
  if (nleqslv.res$termcd != 1){
    print(paste0("Warning! The value is ", as.character(fn(rho))))
  }
  
  #    A = diag(tau^2*rho^(gam-2)*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*rho^(gam-2)*tau^2*E
  A = diag(tau^2*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*tau^2*E
  M = solve(A)- solve(A)%*%B%*%SOLVE(t(B)%*%solve(A)%*%B)%*%t(B)%*%solve(A)
  #  q = tau*rho^(gam-1)*d
  q = tau*rho*d
  Mq = M%*%q
  Z = Mq/norm(Mq, "2")
  B = E%*%cbind(t(V)%*%C, Z)
  rho0 = rho
  
  while (ncol(Z) < om & rcond(t(B)%*%solve(A)%*%B) > 1e-10){
    nleqslv.res = nleqslv(rho0, fn)
    rho = nleqslv.res$x
    if (nleqslv.res$termcd != 1){
      print(paste0("Warning! The value is ", as.character(fn(rho))))
      break
    }
    #    A = diag(tau^2*rho^(gam-2)*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*rho^(gam-2)*tau^2*E
    A = diag(tau^2*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*tau^2*E
    if (rcond(t(B)%*%solve(A)%*%B) < 1e-10){
      break
    }
    M = solve(A)- solve(A)%*%B%*%SOLVE(t(B)%*%solve(A)%*%B)%*%t(B)%*%solve(A)
    #    q = tau*rho^(gam-1)*d
    q = tau*rho*d
    Mq = M%*%q
    z = Mq/norm(Mq, "2")
    Z = cbind(Z, z)
    B = E%*%cbind(t(V)%*%C, Z)
    rho0 = rho
  }
  C = V%*%Z
  #  print(C)
  C = C[,0:min(ncol(C), om)]
  return(as.matrix(C))
}

continuum.int.iter = function(X.list, Y.list, lambda, gam, rankJ, rankA, maxiter = 1000, conv = 1e-6){
  G = length(X.list)
  X = do.call(rbind, X.list)
  Y = do.call(rbind, Y.list)
  p = ncol(X)
  nrun = 0
  converged = F
  
  Y.homo = Y
  Y.homo.list = Y.list
  C = matrix(0, p, rankJ)
  Cind = lapply(1:G, function(g) matrix(0, p, rankA[g]))
  Cind_tot = do.call(cbind, Cind)
  MSE = list()
  U = list()
  W = list()
  
  while (nrun < maxiter & !converged){
    Clast = C
    Cind_tot.last = Cind_tot
    C = continuum.ridge.fix(X, Y.homo, lambda = lambda, gam = gam, om = rankJ)
    U = list.append(U, C)
    beta.C = C2beta(X, Y.homo, C, lambda = lambda)$coef
    Yhat.homo.list = lapply(1:G, function(g) cbind(1, X.list[[g]])%*%beta.C)
    MSE = list.append(MSE, sapply(1:G, function(g) mean((Y.homo.list[[g]] - Yhat.homo.list[[g]])^2)))
    Y.heter.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.homo.list[[g]])
    Cind = lapply(1:G, function(g) continuum.ridge.res.fix(X.list[[g]], Y.heter.list[[g]], C, lambda = lambda, gam = gam, om = rankA[g]))
    W = list.append(W, Cind)
    Cind_tot = do.call(cbind, Cind)
    beta.Cind = lapply(1:G, function(g) C2beta(X.list[[g]], Y.heter.list[[g]], Cind[[g]], lambda)$coef)
    Yhat.heter.list = lapply(1:G, function(g) cbind(1, X.list[[g]])%*%beta.Cind[[g]])
    MSE = list.append(MSE, sapply(1:G, function(g) mean((Y.heter.list[[g]]-Yhat.heter.list[[g]])^2)))
    Y.homo.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.heter.list[[g]])
    Y.homo = do.call(rbind, Y.homo.list)
    print(norm(Clast - C, type = "f"))
    print(norm(Cind_tot.last - Cind_tot, type = "f"))
    if (norm(Clast - C, type = "f") <= conv & norm(Cind_tot.last - Cind_tot, type = "f") <= conv) {
      converged <- T
    }
    nrun = nrun + 1
  }
  if (converged) {
    cat(paste("Algorithm converged after ", nrun, 
              " iterations.\n"))
  }
  else {
    cat(paste("Algorithm did not converge after ", 
              nrun, " iterations.\n"))
  }
  MSE = do.call(rbind, MSE)
  return(list(C = C, Cind = Cind, beta.C = beta.C, beta.Cind = beta.Cind, MSE = MSE, U = U, W = W))
}









