# main function
# continuum ridge using exact solution
# input S is not regularized; lambda show up in function to be optimized
continuum.ridge.S0 = function(X, Y, lambda, gamma = 0){
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
  
  #  e = rev(eigen(S)$values)[1:m]
  e = eigen(S)$values[1:m]
  #  V = eigen(S)$vectors[, rev(seq(1,m))]
  V = eigen(S)$vectors[, seq(1,m)]
  d = t(V)%*%s
  
  fn = function(rho, gamma, lambda, e, V, d){
    lambda0 = gamma + lambda*(1-gamma)/rho
    #    D = diag(lambda0*rho, m) + diag((1-gamma)*e)
    D = diag(lambda0*rho, m) + diag((1-gamma)*e)
    Md = solve(D)%*%d
    z = Md/norm(Md, "2")
    return(t(z)%*%diag(e)%*%z + lambda - rho)
  }
  #  print(nleqslv(mean(e), fn, gamma = gamma, e = e, V = V, d = d))
  nleqslv.res = nleqslv((e[1]+e[m])/2, fn, gamma = gamma, lambda = n*lambda, e = e, V = V, d = d, 
                        method = "Broyden", global = NULL, 
                        control = list(xtol=1e-5, ftol = 1e-5, maxit = 100))
  #   print(nleqslv.res$termcd)
  rho = nleqslv.res$x
  if (nleqslv.res$termcd != 1){
    print(paste0("Warning! The value is ", as.character(fn(rho, gamma, lambda, e, V, d))))
    #     print(fn(rho, gamma, e, V, d))
    #     print(nleqslv.res$termcd)
  }
  lambda0 = gamma + n*lambda*(1-gamma)/rho
  #  D = diag(lambda0*rho, m) + diag((1-gamma)*e)
  D = diag(lambda0*rho, m) + diag((1-gamma)*e)
  Md = solve(D)%*%d
  Z = Md/norm(Md, "2")
  C = V%*%Z
  A = diag(e)%*%t(V)%*%C
  rho0 = rho
  
  fn = function(rho, gamma, lambda, e, V, d, A){
    lambda0 = gamma + lambda*(1-gamma)/rho
    #    D = diag(lambda0*rho, m) + diag((1-gamma)*e)
    D = diag(lambda0*rho, m) + diag((1-gamma)*e)
    M = solve(D)-solve(D)%*%A%*%solve(t(A)%*%solve(D)%*%A)%*%t(A)%*%solve(D)
    Md = M%*%d
    z = Md/norm(Md, "2")
    return(t(z)%*%diag(e)%*%z + lambda - rho)
  }
  #  print(A)
  while (1/kappa(t(A)%*%solve(D)%*%A) > 1e-7){
    #    print(nleqslv(rho, fn, gamma = gamma, e = e, V = V, d = d, A = A))
    nleqslv.res = nleqslv(rho0, fn, gamma = gamma, lambda = n*lambda, e = e, V = V, d = d, A = A, 
                          method = "Broyden", global = NULL, 
                          control = list(xtol=1e-5, ftol = 1e-5, maxit = 100))
    rho = nleqslv.res$x
    #    print(rho)
    #     print(nleqslv.res$termcd)
    if (nleqslv.res$termcd != 1){
      print(paste0("Warning! The value is ", as.character(fn(rho, gamma, lambda, e, V, d, A))))
    }
    if (abs(rho0-rho) < 0.1){
      break
    }
    rho0 = rho
    lambda0 = gamma + n*lambda*(1-gamma)/rho
    #    D = diag(lambda0*rho, m) + diag((1-gamma)*e)
    D = diag(lambda0*rho, m) + diag((1-gamma)*e)
    M = solve(D)-solve(D)%*%A%*%solve(t(A)%*%solve(D)%*%A)%*%t(A)%*%solve(D)
    Md = M%*%d
    z = Md/norm(Md, "2")
    Z = cbind(Z, z)
    A = cbind(A, diag(e)%*%z)
    #    print(A)
    if (ncol(Z) == m-1){
      break
    }
  }
  C = V%*%Z
  #  print(C)
  C = C[,1:min(ncol(C), m-1)]
  return(as.matrix(C))
}


# main function
# continuum ridge using exact solution
# input S is not regularized; lambda show up in function to be optimized
continuum.ridge.S0.rev = function(X, Y, lambda, gamma = 0){
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
  
  e = rev(eigen(S)$values)[1:m]
  # e = eigen(S)$values[1:m]
  V = eigen(S)$vectors[, rev(seq(1,m))]
  # V = eigen(S)$vectors[, seq(1,m)]
  d = t(V)%*%s
  
  fn = function(rho, gamma, lambda, e, V, d){
    lambda0 = gamma + lambda*(1-gamma)/rho
    #    D = diag(lambda0*rho, m) + diag((1-gamma)*e)
    D = diag(lambda0*rho, m) + diag((1-gamma)*e)
    Md = solve(D)%*%d
    z = Md/norm(Md, "2")
    return(t(z)%*%diag(e)%*%z + lambda - rho)
  }
  #  print(nleqslv(mean(e), fn, gamma = gamma, e = e, V = V, d = d))
  nleqslv.res = nleqslv((e[1]+e[m])/2, fn, gamma = gamma, lambda = n*lambda, e = e, V = V, d = d, 
                        method = "Broyden", global = NULL, 
                        control = list(xtol=1e-5, ftol = 1e-5, maxit = 100))
  #   print(nleqslv.res$termcd)
  rho = nleqslv.res$x
  if (nleqslv.res$termcd != 1){
    print(paste0("Warning! The value is ", as.character(fn(rho, gamma, lambda, e, V, d))))
    #     print(fn(rho, gamma, e, V, d))
    #     print(nleqslv.res$termcd)
  }
  lambda0 = gamma + n*lambda*(1-gamma)/rho
  #  D = diag(lambda0*rho, m) + diag((1-gamma)*e)
  D = diag(lambda0*rho, m) + diag((1-gamma)*e)
  Md = solve(D)%*%d
  Z = Md/norm(Md, "2")
  C = V%*%Z
  A = diag(e)%*%t(V)%*%C
  rho0 = rho
  
  fn = function(rho, gamma, lambda, e, V, d, A){
    lambda0 = gamma + lambda*(1-gamma)/rho
    #    D = diag(lambda0*rho, m) + diag((1-gamma)*e)
    D = diag(lambda0*rho, m) + diag((1-gamma)*e)
    M = solve(D)-solve(D)%*%A%*%solve(t(A)%*%solve(D)%*%A)%*%t(A)%*%solve(D)
    Md = M%*%d
    z = Md/norm(Md, "2")
    return(t(z)%*%diag(e)%*%z + lambda - rho)
  }
  #  print(A)
  while (1/kappa(t(A)%*%solve(D)%*%A) > 1e-7){
    #    print(nleqslv(rho, fn, gamma = gamma, e = e, V = V, d = d, A = A))
    nleqslv.res = nleqslv(rho0, fn, gamma = gamma, lambda = n*lambda, e = e, V = V, d = d, A = A, 
                          method = "Broyden", global = NULL, 
                          control = list(xtol=1e-5, ftol = 1e-5, maxit = 100))
    rho = nleqslv.res$x
    #    print(rho)
    #     print(nleqslv.res$termcd)
    if (nleqslv.res$termcd != 1){
      print(paste0("Warning! The value is ", as.character(fn(rho, gamma, lambda, e, V, d, A))))
    }
    if (abs(rho0-rho) < 0.1){
      break
    }
    rho0 = rho
    lambda0 = gamma + n*lambda*(1-gamma)/rho
    #    D = diag(lambda0*rho, m) + diag((1-gamma)*e)
    D = diag(lambda0*rho, m) + diag((1-gamma)*e)
    M = solve(D)-solve(D)%*%A%*%solve(t(A)%*%solve(D)%*%A)%*%t(A)%*%solve(D)
    Md = M%*%d
    z = Md/norm(Md, "2")
    Z = cbind(Z, z)
    A = cbind(A, diag(e)%*%z)
    #    print(A)
    if (ncol(Z) == m-1){
      break
    }
  }
  C = V%*%Z
  #  print(C)
  C = C[,1:min(ncol(C), m-1)]
  return(as.matrix(C))
}


# input S is directly regularized
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
  
  S = t(X)%*%X + n*lambda*diag(p)
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
  C = C[,1:min(m-1, ncol(C))]
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
  
  S = t(X)%*%X + n*lambda*diag(p)
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
  #  print(C)
  #  C = C[,-ncol(C)]
  C = C[,1:min(m-1, ncol(C))]
  return(as.matrix(C))
}