library(quadprog)
library(nleqslv)
library(Matrix)
library(rlist)

SOLVE = function(x){
  if (sum(dim(x))){
    return(solve(x))
  }else{
    return(x)
  }
}

DIAG = function(e){
  if (length(e) > 1){
    return(diag(e))
  }else{
    return(matrix(e))
  }
}

# C output from continuum.ridge
C2beta = function(X, Y, C, lambda){
  n = nrow(X)
  X.mean = apply(X, 2, mean)
  Y.mean = mean(Y)
  S = t(X)%*%X
  s = t(X)%*%Y
  
  C = matrix(C, nrow = ncol(X))
  om = ncol(C) #omega
  
  alpha.C = SOLVE(t(C)%*%S%*%C + n*lambda*diag(om))%*%t(C)%*%s
  beta.C = C%*%alpha.C
  intercept = Y.mean - X.mean%*%beta.C
  
  return(list(intercept = intercept, beta = beta.C, alpha = alpha.C, coef = matrix(c(intercept, beta.C))))
}

# decomposeX = function(X.list, U, Wind, centerValues.X, scaleValues.X){
#   G = length(X.list)
#   n = ncol(X.list[[1]])
#   p.list = lapply(1:G, function(g) nrow(X.list[[g]]))
#   p = do.call("+", p.list)
#   r = ncol(U)
#   r.list = lapply(1:G, function(g) ncol(Wind[[g]]))
#   rr = do.call("+", r.list)
#   
#   X.list = lapply(1:G, function(g) sweep(X.list[[g]], 1, centerValues.X[[g]])/scaleValues.X[[g]])
#   X = do.call(rbind, X.list)
#   W = bdiag(Wind)
#   A = SOLVE(t(U)%*%U)%*%t(U)
#   B = SOLVE(t(W)%*%W)%*%t(W)
#   SS = SOLVE(diag(rr) - B%*%U%*%A%*%W)%*%B%*%(diag(p) - U%*%A)%*%X
#   Sind = list()
#   for (g in 1:G){
#     Sind[[g]] = t(matrix(SS[1:r.list[[g]],], ncol = n))
#     SS = matrix(SS[-(1:r.list[[g]]),], ncol = n)
#   }
#   S = as.matrix(t(SOLVE(diag(r) - A%*%W%*%B%*%U)%*%A%*%(diag(p) - W%*%B)%*%X))
#   J = U%*%t(S)
#   I = lapply(1:G, function(g) Wind[[g]]%*%t(Sind[[g]]))
#   return(list(J = J, I = I, S = S, Sind = Sind))
# }
# 
# decomposeX = function(X.list, C, Cind, centerValues.X, scaleValues.X){
#   G = length(X.list)
#   n = ncol(X.list[[1]])
#   p.list = lapply(1:G, function(g) nrow(X.list[[g]]))
#   p = do.call("+", p.list)
#   
#   X.list = lapply(1:G, function(g) sweep(X.list[[g]], 1, centerValues.X[[g]])/scaleValues.X[[g]])
#   X = do.call(rbind, X.list)
#   J = X%*%C%*%SOLVE(t(C)%*%C)%*%t(C)
#   I = X - J
#   Iind = list()
#   for (g in 1:G){
#     Iind[[g]] = I[1:p.list[[g]],]%*%Cind[[g]]%*%SOLVE(t(Cind[[g]])%*%Cind[[g]])%*%t(Cind[[g]])
#     I = I[-(1:p.list[[g]]),]
#   }
#   return(list(J = J, I = Iind))
# }


continuum.multigroup.iter = function(X.list, Y.list, lambda, gam, rankJ, rankA, maxiter = 1000, conv = 1e-6, 
                                     center.X = TRUE, scale.X = TRUE, center.Y = TRUE, scale.Y = TRUE, orthIndiv = TRUE){
  G = length(X.list)
  centerValues.X <- list()
  scaleValues.X <- list()
  centerValues.Y <- list()
  scaleValues.Y <- list()
  n = c()
  for (g in 1:G){
    n[g] = nrow(X.list[[g]])
  }
  N = sum(n)
  p = ncol(X.list[[1]])
  
  for (g in 1:G){
    # X
    if (center.X){
      centerValues.X[[g]] = apply(X.list[[g]], 2, mean)
    }else{
      centerValues.X[[g]] = rep(0, p)
    }
    if (scale.X){
      scaleValues.X[[g]] = norm(X.list[[g]], type = "f") #*sqrt(N*p)
    }else{
      scaleValues.X[[g]] = 1
    }
    X.list[[g]] = sweep(X.list[[g]], 2, centerValues.X[[g]])
    X.list[[g]] = X.list[[g]]/scaleValues.X[[g]]
    
    # Y
    if (center.Y){
      centerValues.Y[[g]] = mean(Y.list[[g]])
    }else{
      centerValues.Y[[g]] = 0
    }
    if (scale.Y){
      scaleValues.Y[[g]] = norm(Y.list[[g]], type = "f") #*sqrt(N)
    }else{
      scaleValues.Y[[g]] = 1
    }
    Y.list[[g]] = sweep(matrix(Y.list[[g]]), 2, centerValues.Y[[g]])
    Y.list[[g]] = sweep(matrix(Y.list[[g]]), 2, scaleValues.Y[[g]], FUN = "/")
  }
  X = do.call(rbind, X.list)
  Y = do.call(rbind, Y.list)
  
  nrun = 0
  converged = F
  
  Y.homo = Y
  Y.homo.list = Y.list
  X.heter.list = lapply(1:G, function(g) matrix(0, nrow = nrow(X.list[[g]]), ncol = p))
#  C = matrix(0, p, rankJ)
  Cind = lapply(1:G, function(g) matrix(0, p, rankA[g]))
#  Cind_tot = do.call(cbind, Cind)
  R = X
  R[,] = 0
  r = Y
  r[,] = 0
  MSEcurrent = 0
  MSE = list()
  U = list()
  W = list()
  ct.homo = matrix(0, nrow = rankJ, ncol = 1)
  ct.heter = lapply(1:G, function(g) matrix(0, nrow = rankA[g], ncol = 1))
  
  while (nrun < maxiter & !converged){
    # initialization
    rlast = r
    Rlast = R
    MSElast = MSEcurrent
    ct.homo.last = ct.homo
    ct.heter.last = ct.heter

#    Clast = C
#    Cind_tot.last = Cind_tot
    X.homo.list = lapply(1:G, function(g) X.list[[g]] - X.heter.list[[g]])
    X.homo = do.call(rbind, X.homo.list)
    
    # joint 
    ml.homo = continuum.ridge.fix(X.homo, Y.homo, G, lambda = lambda, gam = gam, om = rankJ)
    C = ml.homo$C
    U = list.append(U, C)
    
    XC = X.homo%*%C
    ct.homo = diag(t(XC)%*%XC)^gam*cor(XC, Y.homo)
    
    beta.C = C2beta(X.homo, Y.homo, C, lambda = lambda)$beta
    Yhat.homo.list = lapply(1:G, function(g) X.homo.list[[g]]%*%beta.C)
    Y.heter.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.homo.list[[g]])
    X.homo.list = lapply(1:G, function(g) X.homo.list[[g]]%*%C%*%SOLVE(t(C)%*%C)%*%t(C))
    X.heter.list = lapply(1:G, function(g) X.list[[g]] - X.homo.list[[g]])
    
    # individual
    temp = X.heter.list
    for (g in 1:G){
      temp[[g]] <- temp[[g]]%*%(diag(p) - C%*%SOLVE(t(C)%*%C)%*%t(C))
      
      # orthogonalization
      if (orthIndiv){
        if (nrun > 0){
          for (j in (1:G)[-g]){
            temp[[g]] <- temp[[g]]%*%(diag(p) - Cind[[j]]%*%SOLVE(t(Cind[[j]])%*%Cind[[j]])%*%t(Cind[[j]]))
          }
        }
      }
      
      if (rankA[g]){
        ml.heter = continuum.ridge.fix(temp[[g]], Y.heter.list[[g]], 1, lambda = lambda, gam = gam, om = rankA[g])
        Cind[[g]] = ml.heter$C
      }
      
#      ml.heter = continuum.ridge.fix(temp[[g]], Y.heter.list[[g]], 1, lambda = lambda, gam = gam, om = rankA[g])
#      Cind[[g]] = ml.heter$C
    }
    W = list.append(W, Cind)
    
    XC.list = lapply(1:G, function(g) X.heter.list[[g]]%*%Cind[[g]])
    ct.heter = lapply(1:G, function(g) diag(t(XC.list[[g]])%*%XC.list[[g]])^gam*cor(XC.list[[g]], Y.heter.list[[g]]))

    #    Cind_tot = do.call(cbind, Cind)
    beta.Cind = lapply(1:G, function(g) C2beta(X.heter.list[[g]], Y.heter.list[[g]], Cind[[g]], lambda)$beta)
    Yhat.heter.list = lapply(1:G, function(g) X.heter.list[[g]]%*%beta.Cind[[g]])
    
    # MSE
    r.list = lapply(1:G, function(g) Y.heter.list[[g]] - Yhat.heter.list[[g]])
    r = do.call(rbind, r.list)
    MSEcurrent = sum(sapply(1:G, function(g) sum(r.list[[g]]^2)))/N
    MSE = list.append(MSE, mean(r.list[[g]]^2))
#    print(MSEcurrent)
    
    Y.homo.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.heter.list[[g]])
    Y.homo = do.call(rbind, Y.homo.list)
    X.heter.list = lapply(1:G, function(g) X.heter.list[[g]]%*%Cind[[g]]%*%SOLVE(t(Cind[[g]])%*%Cind[[g]])%*%t(Cind[[g]]))
    
    if (orthIndiv){
      temp = X.heter.list
      if (nrun == 0){
        for (g in 1:G){
          for (j in (1:G)[-g]) {
            temp[[g]] <- temp[[g]]%*%(diag(p) - Cind[[j]]%*%SOLVE(t(Cind[[j]])%*%Cind[[j]])%*%t(Cind[[j]]))
          }
        }
        for (g in 1:G){
          if (rankA[g]){
            ml.heter = continuum.ridge.fix(temp[[g]], Y.heter.list[[g]], 1, lambda = lambda, gam = gam, om = rankA[g])
            Cind[[g]] = ml.heter$C
          }
        }
#        ml.heter = lapply(1:G, function(g) continuum.ridge.fix(temp[[g]], Y.heter.list[[g]], 1, lambda = lambda, gam = gam, om = rankA[g]))
#        Cind = lapply(1:G, function(g) ml.heter[[g]]$C)
      }
    }
    
    
    R.list = lapply(1:G, function(g) X.list[[g]] - X.homo.list[[g]] - X.heter.list[[g]])
    R = do.call(rbind, R.list)

    # if (norm(Clast - C, type = "f") <= conv & norm(Cind_tot.last - Cind_tot, type = "f") <= conv) {
    #  converged <- T
    # }
    
#    print(abs(MSEcurrent - MSElast))
    # if (abs(MSEcurrent - MSElast) <= conv){
    #   converged <- T
    # }

    if (gam > 10){
      if (!nrun%%10){
        print(norm(Rlast - R, type = "f"))
      }
      if (norm(Rlast - R, type = "f") <= conv){
        converged <- T
      }
    }else{
      # CT1 = norm(ct.homo.last/ct.homo - 1, type = "f")
      # CT2 = sapply(1:G, function(g) norm(ct.heter.last[[g]]/ct.heter[[g]] - 1, type = "f"))
      # if (!nrun%%10){
      #   print(CT1)
      #   print(CT2)
      # }
      # if (CT1 <= conv & sum(CT2) <= conv) {
      #  converged <- T
      # }
      
      if (!nrun%%10){
        print(norm(rlast - r, type = "f"))
      }
      if (norm(rlast - r, type = "f") <= conv){
        converged <- T
      }
    }
    
    # if (!nrun%%10){
    #   print(norm(rlast - r, type = "f"))
    # }
    # if (norm(rlast - r, type = "f") <= conv){
    #   converged <- T
    # }
    
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
  beta.C0 = beta.C
  beta.Cind0 = beta.Cind
  beta.C = list()
  beta.Cind = list()
  intercept = list()
  for (g in 1:G){
    beta.C[[g]] = beta.C0/scaleValues.X[[g]]*scaleValues.Y[[g]]
    beta.Cind[[g]] = beta.Cind0[[g]]/scaleValues.X[[g]]*scaleValues.Y[[g]]
    intercept[[g]] = centerValues.Y[[g]] - t(beta.C[[g]])%*%centerValues.X[[g]] - t(beta.Cind[[g]])%*%centerValues.X[[g]]
  }
  return(list(C = C, Cind = Cind, intercept = intercept, beta.C = beta.C, beta.Cind = beta.Cind, MSE = MSE, U = U, W = W, converged = converged,
              J = X.homo.list, I = X.heter.list, nrun = nrun))
}

continuum.multisource.iter = function(X.list, Y, lambda, gam, rankJ, rankA, maxiter = 1000, conv = 1e-6, 
                                      center.X = TRUE, scale.X = TRUE, center.Y = TRUE, scale.Y = TRUE){
  G = length(X.list)
  # preprocess X
  centerValues.X <- list()
  scaleValues.X <- list()
  p = c()
  for (g in 1:G){
    p[g] = nrow(X.list[[g]])
    if (center.X){
      centerValues.X[[g]] = apply(X.list[[g]], 1, mean)
    }else{
      centerValues.X[[g]] = rep(0, p[g])
    }
    if (scale.X){
      scaleValues.X[[g]] = norm(X.list[[g]], type = "f")
    }else{
      scaleValues.X[[g]] = 1
    }
    X.list[[g]] = sweep(X.list[[g]], 1, centerValues.X[[g]])
    X.list[[g]] = X.list[[g]]/scaleValues.X[[g]]
  }
  
  # preprocess Y
  if (center.Y){
    centerValues.Y = mean(Y)
  }else{
    centerValues.Y = 0
  }
  if (scale.Y){
    scaleValues.Y = norm(matrix(Y), type = "f")
  }else{
    scaleValues.Y = 1
  } 
  Y = sweep(matrix(Y), 2, centerValues.Y)
  Y = sweep(matrix(Y), 2, scaleValues.Y, FUN = "/")
  
  X = do.call(rbind, X.list)
  n = ncol(X)
  P = sum(p)
  
  nrun = 0
  converged = F
  
  Y.homo = Y
  X.heter.list = lapply(1:G, function(g) matrix(0, nrow = p[g], ncol = n))
  X.heter = do.call(rbind, X.heter.list)
#  V = matrix(0, P, rankJ)
#  Vind = lapply(1:G, function(g) matrix(0, p[g], rankA[g]))
  C = matrix(0, P, rankJ)
  Cind = lapply(1:G, function(g) matrix(0, p[g], rankA[g]))
  ml.heter = list()
  R = X
  R[,] = 0
  MSE = list()
  U = list()
  W = list()
  S = list()
  Sind = list()
  mse = 0
  
  while (nrun < maxiter & !converged){
    mselast = mse
    Rlast = R
    X.homo.list = lapply(1:G, function(g) X.list[[g]] - X.heter.list[[g]])
    X.homo = do.call(rbind, X.homo.list)
    ml.homo = continuum.ridge.fix(t(X.homo), Y.homo, G, lambda = lambda, gam = gam, om = rankJ, vertical = FALSE)
    C = ml.homo$C
    a = ml.homo$a
    S = list.append(S, C)
    U. = X.homo%*%C%*%SOLVE(t(C)%*%C)
    U = list.append(U, U.)
#    ml.jive$S[[nrun+1]]/C
    beta.C = C2beta(t(X.homo)%*%X.homo, Y.homo, a, lambda = lambda)$beta
    Yhat.homo = t(X.homo)%*%X.homo%*%beta.C
    Y.heter = Y - Yhat.homo
    X.homo = U.%*%t(C)
#    X.homo = U.%*%t(V)
    X.heter.list = list()
    X.homo.list = list()
    for (g in 1:G){
      X.heter.list[[g]] = X.list[[g]] - X.homo[1:p[g],]
      X.homo.list[[g]] = X.homo[1:p[g],]
      X.homo = X.homo[-(1:p[g]),]
    }
    for (g in 1:G){
      X.heter.list[[g]] <- X.heter.list[[g]]%*%(diag(n) - C%*%SOLVE(t(C)%*%C)%*%t(C))
#      X.heter.list[[g]] <- X.heter.list[[g]]%*%(diag(n) - V%*%SOLVE(t(V)%*%V)%*%t(V))
      if (nrun > 0){
        for (j in (1:G)[-g]){
          X.heter.list[[g]] <- X.heter.list[[g]]%*%(diag(n) - Cind[[j]]%*%SOLVE(t(Cind[[j]])%*%Cind[[j]])%*%t(Cind[[j]]))
#          X.heter.list[[g]] <- X.heter.list[[g]]%*%(diag(n) - Vind[[j]]%*%SOLVE(t(Vind[[j]])%*%Vind[[j]])%*%t(Vind[[j]]))
        }
      }
      ml.heter[[g]] = continuum.ridge.fix(t(X.heter.list[[g]]), Y.heter, 1, lambda = lambda, gam = gam, om = rankA[g], 
                                          vertical = FALSE)
      Cind[[g]] = ml.heter[[g]]$C
#      Vind[[g]] = ml.heter[[g]]$V[,1:rankA[g]]
    }
#    t(Cind[[1]])%*%Cind[[2]]
    aind = lapply(1:G, function(g) ml.heter[[g]]$a)
    beta.Cind = lapply(1:G, function(g) C2beta(t(X.heter.list[[g]])%*%X.heter.list[[g]], Y.heter, aind[[g]], lambda)$beta)
    Yhat.heter.list = lapply(1:G, function(g) t(X.heter.list[[g]])%*%X.heter.list[[g]]%*%beta.Cind[[g]])
    Yhat.heter = do.call("+", Yhat.heter.list)
    mse = mean((Y.heter - Yhat.heter)^2)
    Wind = lapply(1:G, function(g) X.heter.list[[g]]%*%Cind[[g]]%*%SOLVE(t(Cind[[g]])%*%Cind[[g]]))
#    Wind = lapply(1:G, function(g) X.heter.list[[g]]%*%Vind[[g]]%*%SOLVE(t(Vind[[g]])%*%Vind[[g]]))
    X.heter.list = lapply(1:G, function(g) Wind[[g]]%*%t(Cind[[g]]))
#    X.heter.list = lapply(1:G, function(g) Wind[[g]]%*%t(Vind[[g]]))
    if (nrun == 0){
      for (g in 1:G){
        for (j in (1:G)[-g]) {
          X.heter.list[[g]] <- X.heter.list[[g]] %*%(diag(n) - Cind[[j]]%*%SOLVE(t(Cind[[j]])%*%Cind[[j]])%*%t(Cind[[j]]))
#          X.heter.list[[g]] <- X.heter.list[[g]] %*%(diag(n) - Vind[[j]]%*%SOLVE(t(Vind[[j]])%*%Vind[[j]])%*%t(Vind[[j]]))
        }
      }
      ml.heter = lapply(1:G, function(g) continuum.ridge.fix(t(X.heter.list[[g]]), Y.heter, 1, lambda = lambda, gam = gam, om = rankA[g],
                                                   vertical = FALSE))
      Cind = lapply(1:G, function(g) ml.heter[[g]]$C)
#      Vind = lapply(1:G, function(g) ml.heter[[g]]$V[,1:rankA[g]])
    }
#    ml.jive$T[[nrun+1]][[1]]/Cind[[1]]
#    ml.jive$T[[nrun+1]][[2]]/Cind[[2]]
    Sind = list.append(Sind, Cind)
    W = list.append(W, Wind)
    MSE = list.append(MSE, mse)
    R.list = lapply(1:G, function(g) X.list[[g]] - X.homo.list[[g]] - X.heter.list[[g]])
    R = do.call(rbind, R.list)
    Y.homo = Y - Yhat.heter
#    print(mse)
    print(norm(Rlast - R, type = "f"))
    if (norm(Rlast - R, type = "f") <= conv) {
      converged <- T
    }
    
    # print(mse)
    # if ((mse - mselast) <= conv) {
    #   converged <- T
    # }
    nrun = nrun + 1
  }
  if (converged) {
    cat(paste("Algorithm converged after ", nrun, 
              " iterations.\n"))
  }else{
    cat(paste("Algorithm did not converge after ", 
              nrun, " iterations.\n"))
  }
  MSE = do.call(rbind, MSE)
  beta.C = beta.C*scaleValues.Y
  beta.Cind = lapply(1:G, function(g) beta.Cind[[g]]*scaleValues.Y)
  return(list(C = C, Cind = Cind, beta.C = beta.C, beta.Cind = beta.Cind, 
              centerValues.X = centerValues.X, centerValues.Y = centerValues.Y,
              scaleValues.X = scaleValues.X, scaleValues.Y = scaleValues.Y,
              intercept = centerValues.Y, 
              MSE = MSE, S = S, T = Sind, U = U, W = W, converged = converged,
              J = do.call(rbind, X.homo.list), I = X.heter.list, nrun = nrun))
}


continuum.2step.v1 = function(X.list, Y.list, lambda, gam, rankJ, rankA, scale.X = T){
  G = length(X.list)
  centerValues.X <- list()
  scaleValues.X <- list()
  centerValues.Y <- list()
  scaleValues.Y <- list()
  
  for (g in 1:G){
    # X
    centerValues.X[[g]] = apply(X.list[[g]], 2, mean)
    X.list[[g]] = sweep(X.list[[g]], 2, centerValues.X[[g]])
    if (scale.X){
      scaleValues.X[[g]] = norm(X.list[[g]], type = "f")
    }else{
      scaleValues.X[[g]] = 1
    }
    X.list[[g]] = X.list[[g]]/scaleValues.X[[g]]
    
    # Y
    centerValues.Y[[g]] = mean(Y.list[[g]])
    Y.list[[g]] = sweep(matrix(Y.list[[g]]), 2, centerValues.Y[[g]])
    scaleValues.Y[[g]] = sqrt(mean(Y.list[[g]]^2))
    Y.list[[g]] = sweep(matrix(Y.list[[g]]), 2, scaleValues.Y[[g]], FUN = "/")
  }
  X = do.call(rbind, X.list)
  Y = do.call(rbind, Y.list)
  p = ncol(X)
  
  ml.homo = continuum.ridge.fix(X, Y, G, lambda = lambda, gam = gam, om = rankJ)
  C = ml.homo$C
  beta.C = C2beta(X, Y, C, lambda = lambda)$beta
  Yhat.homo.list = lapply(1:G, function(g) X.list[[g]]%*%beta.C)
  Y.heter.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.homo.list[[g]])
  X.homo.list = lapply(1:G, function(g) X.list[[g]]%*%C%*%SOLVE(t(C)%*%C)%*%t(C))
  X.heter.list = lapply(1:G, function(g) X.list[[g]] - X.homo.list[[g]])
  ml.heter = lapply(1:G, function(g) continuum.ridge.fix(X.heter.list[[g]], Y.heter.list[[g]], G = 1, lambda = lambda, gam = gam, om = rankA[g]))
  Cind = lapply(1:G, function(g) ml.heter[[g]]$C)
  beta.Cind = lapply(1:G, function(g) C2beta(X.heter.list[[g]], Y.heter.list[[g]], Cind[[g]], lambda)$beta)
  Yhat.heter.list = lapply(1:G, function(g) X.heter.list[[g]]%*%beta.Cind[[g]])
  MSE = sapply(1:G, function(g) mean((Y.heter.list[[g]] - Yhat.heter.list[[g]])^2))
  
  beta.C0 = beta.C
  beta.Cind0 = beta.Cind
  beta.C = list()
  beta.Cind = list()
  intercept = list()
  for (g in 1:G){
    beta.C[[g]] = beta.C0/scaleValues.X[[g]]*scaleValues.Y[[g]]
    beta.Cind[[g]] = beta.Cind0[[g]]/scaleValues.X[[g]]*scaleValues.Y[[g]]
    intercept[[g]] = centerValues.Y[[g]] - t(beta.C[[g]])%*%centerValues.X[[g]] - t(beta.Cind[[g]])%*%centerValues.X[[g]]
  }
  return(list(C = C, Cind = Cind, intercept = intercept, beta.C = beta.C, beta.Cind = beta.Cind, MSE = MSE))
}

continuum.2step.v2 = function(X.list, Y.list, lambda, gam, rankJ, rankA, scale.X = T){
  G = length(X.list)
  centerValues.X <- list()
  scaleValues.X <- list()
  centerValues.Y <- list()
  scaleValues.Y <- list()
  
  for (g in 1:G){
    centerValues.X[[g]] = apply(X.list[[g]], 2, mean)
    X.list[[g]] = sweep(X.list[[g]], 2, centerValues.X[[g]])
    if (scale.X){
      scaleValues.X[[g]] = apply(X.list[[g]], 2, function(X) sqrt(mean(X^2))) 
      X.list[[g]] = sweep(X.list[[g]], 2, scaleValues.X[[g]], FUN = "/")
    }else{
      scaleValues.X[[g]] = apply(X.list[[g]], 2, function(X) 1)
    }
    centerValues.Y[[g]] = mean(Y.list[[g]])
    Y.list[[g]] = sweep(matrix(Y.list[[g]]), 2, centerValues.Y[[g]])
    scaleValues.Y[[g]] = sqrt(mean(Y.list[[g]]^2))
    Y.list[[g]] = sweep(matrix(Y.list[[g]]), 2, scaleValues.Y[[g]], FUN = "/")
  }
  X = do.call(rbind, X.list)
  Y = do.call(rbind, Y.list)
  p = ncol(X)
  
  ml.homo = continuum.ridge.fix(X, Y, G, lambda = lambda, gam = gam, om = rankJ)
  C = ml.homo$C
  beta.C = C2beta(X, Y, C, lambda = lambda)$beta
  Yhat.homo.list = lapply(1:G, function(g) X.list[[g]]%*%beta.C)
  Y.heter.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.homo.list[[g]])
  X.homo.list = lapply(1:G, function(g) X.list[[g]]%*%C%*%SOLVE(t(C)%*%C)%*%t(C))
  X.heter.list = lapply(1:G, function(g) X.list[[g]] - X.homo.list[[g]])
  U.homo.list = lapply(1:G, function(g) svd(X.homo.list[[g]], nu = rankJ, nv = rankJ)$u)
  ml.heter = lapply(1:G, function(g) if (is.null(U.homo.list[[g]])) {
    continuum.ridge.fix(X.heter.list[[g]], Y.heter.list[[g]], G = 1, lambda = lambda, gam = gam, om = rankA[g])
  }else{
    continuum.ridge.res.fix(X.heter.list[[g]], Y.heter.list[[g]], U.homo.list[[g]], lambda = lambda, gam = gam, om = rankA[g])
  })
  Cind = lapply(1:G, function(g) ml.heter[[g]]$C)
#  t(X.homo.list[[g]])%*%X.heter.list[[g]]%*%Cind[[g]]
#  t(C)%*%Cind[[g]]
  beta.Cind = lapply(1:G, function(g) C2beta(X.heter.list[[g]], Y.heter.list[[g]], Cind[[g]], lambda)$beta)
  Yhat.heter.list = lapply(1:G, function(g) X.heter.list[[g]]%*%beta.Cind[[g]])
  MSE = sapply(1:G, function(g) mean((Y.heter.list[[g]] - Yhat.heter.list[[g]])^2))
#  print(MSE)
#  X.heter.list = lapply(1:G, function(g) X.heter.list[[g]]%*%Cind[[g]]%*%solve(t(Cind[[g]])%*%Cind[[g]])%*%t(Cind[[g]]))

  beta.C0 = beta.C
  beta.Cind0 = beta.Cind
  beta.C = list()
  beta.Cind = list()
  intercept = list()
  for (g in 1:G){
    beta.C[[g]] = beta.C0/scaleValues.X[[g]]*scaleValues.Y[[g]]
    beta.Cind[[g]] = beta.Cind0[[g]]/scaleValues.X[[g]]*scaleValues.Y[[g]]
    intercept[[g]] = centerValues.Y[[g]] - t(beta.C[[g]])%*%centerValues.X[[g]] - t(beta.Cind[[g]])%*%centerValues.X[[g]]
  }
  return(list(C = C, Cind = Cind, intercept = intercept, beta.C = beta.C, beta.Cind = beta.Cind, MSE = MSE))
}

continuum.ridge.fix = function(X, Y, G, lambda, gam, om, vertical = TRUE){
  #om: number of columns
  n = nrow(X)
  p = ncol(X)
  
#  scaleValues.Y = sqrt(mean(Y^2))
#  Y = sweep(matrix(Y), 2, scaleValues.Y, FUN = "/")
  
  if (vertical){
    S = t(X)%*%X
    s = t(X)%*%Y
    e = eigen(S)$values
    m = length(e[e > 1e-6])
    e = e[1:m]
#    m = min(n-G, p)
#    e = eigen(S)$values[1:m]
    E = DIAG(e)
    V = eigen(S)$vectors[, 1:m]
    d = t(V)%*%s
  }else{
    S = X%*%t(X)
#    m = min(n-G, p)
    e = eigen(S)$values
    m = length(e[e > 1e-6])
    e = e[1:m]
#    e = eigen(S)$values[1:m]
    E = DIAG(e)
    V = eigen(S)$vectors[, 1:m]
    d = E^(1/2)%*%t(V)%*%Y
  }

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
  a = V%*%solve(E)^(1/2)%*%Z
  return(list(C = as.matrix(C), a = a, V = as.matrix(V), Z = Z, E = E))
}


continuum.ridge.res.fix = function(X, Y, Uhomo, lambda, gam, om){
  #om: number of columns
  n = nrow(X)
  p = ncol(X)
  
  scaleValues.Y = sqrt(mean(Y^2))
  Y = sweep(matrix(Y), 2, scaleValues.Y, FUN = "/")
  
  S = t(X)%*%X
  s = t(X)%*%Y
  m = min(n-1, p)
  
  SVD = svd(X, nu = m, nv = m)
  e = SVD$d[1:m]^2
  E = diag(SVD$d[1:m])
  V = SVD$v
  U = SVD$u
  d = t(V)%*%s
  B = E%*%t(U)%*%Uhomo
  
  #  e = eigen(S)$values[1:m]
  #  E = diag(e)
  #  V = eigen(S)$vectors[, seq(1,m)]
  #  d = t(V)%*%s
  #  B = E%*%t(V)%*%C
  
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
  B = E%*%cbind(t(U)%*%Uhomo, Z)
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
  return(list(C = as.matrix(C), V = as.matrix(V)))
}

continuum.ridge.fixm = function(X, Y, G, lambda, gam, om, m, vertical = TRUE){
  #om: number of columns
  n = nrow(X)
  p = ncol(X)
  
  #  scaleValues.Y = sqrt(mean(Y^2))
  #  Y = sweep(matrix(Y), 2, scaleValues.Y, FUN = "/")
  
  if (vertical){
    S = t(X)%*%X
    s = t(X)%*%Y
    e = eigen(S)$values
    #    m = length(e[e > 1e-12])
    e = e[1:m]
    #    m = min(n-G, p)
    #    e = eigen(S)$values[1:m]
    E = DIAG(e)
    V = eigen(S)$vectors[, 1:m]
    d = t(V)%*%s
  }else{
    S = X%*%t(X)
    #    m = min(n-G, p)
    e = eigen(S)$values
    #    m = length(e[e > 1e-12])
    e = e[1:m]
    #    e = eigen(S)$values[1:m]
    E = DIAG(e)
    V = eigen(S)$vectors[, 1:m]
    d = E^(1/2)%*%t(V)%*%Y
  }
  
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
  a = V%*%solve(E)^(1/2)%*%Z
  return(list(C = as.matrix(C), a = a, V = as.matrix(V), Z = Z, E = E))
}



