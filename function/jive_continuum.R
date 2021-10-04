library(quadprog)
library(nleqslv)
library(Matrix)
library(rlist)
library(MASS)

SOLVE = function(x){
  if (sum(dim(x))){
    return(ginv(x))
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

decomposeX.iter = function(X.list, C, Cind, centerValues.X, scaleValues.X, maxiter = 1000, conv = 1e-6){
  G = length(X.list)
  p = ncol(X.list[[1]])
  
  X.list = lapply(1:G, function(g) sweep(X.list[[g]], 2, centerValues.X[[g]])/scaleValues.X[[g]])
  X = do.call(rbind, X.list)
  
  X.heter.list = lapply(1:G, function(g) matrix(0, nrow = nrow(X.list[[g]]), ncol = p))
  R = X
  nrun = 0
  converged = F
  
  while (nrun < maxiter & !converged){
    Rlast = R
    X.homo.list = lapply(1:G, function(g) X.list[[g]] - X.heter.list[[g]])
    X.homo.list = lapply(1:G, function(g) X.homo.list[[g]]%*%C%*%SOLVE(t(C)%*%C)%*%t(C))
    X.heter.list = lapply(1:G, function(g) X.list[[g]] - X.homo.list[[g]])
    X.heter.list = lapply(1:G, function(g) 
      X.heter.list[[g]]%*%Cind[[g]]%*%SOLVE(t(Cind[[g]])%*%Cind[[g]])%*%t(Cind[[g]]))
    R.list = lapply(1:G, function(g) X.list[[g]] - X.homo.list[[g]] - X.heter.list[[g]])
    R = do.call(rbind, R.list)
    #    if (!nrun%%10){
    print(norm(Rlast - R, type = "f"))
    #    }
    if (norm(Rlast - R, type = "f") <= conv){
      converged <- T
    }
    nrun = nrun + 1
  }
  
  return(list(J = X.homo.list, I = X.heter.list, nrun = nrun, converged = converged, R = R))
  
}

continuum = function(X, Y, G, lambda, gam, om, 
                     U_old=matrix(,nrow=nrow(X), ncol=0), D_old=matrix(,nrow=0, ncol=0), V_old=matrix(,nrow=0, ncol=0), Z_old=matrix(,nrow=0, ncol=0), 
                     vertical = TRUE, verbose = FALSE){
  n = nrow(X)
  p = ncol(X)
  
  if (vertical){
    s = t(X)%*%Y
    # svd.X = svd(X)
    # d = svd.X$d
    # V = svd.X$v
    # U = svd.X$u
    # m = rankMatrix(X)[1]
    # e = (d^2)[1:m]
    # D = DIAG(d[1:m])
    # E = DIAG(e)
    # V = V[,1:m]
    # U = U[,1:m]
    UDVZ = initialize.UDVZ(X)
    D = UDVZ$D
    E = UDVZ$E
    V = UDVZ$V
    U = UDVZ$U
    e = UDVZ$e
    m = UDVZ$m
    d = t(V)%*%s
    E2 = D%*%t(U)%*%U_old%*%D_old
  }else{
    svd.tX = svd(t(X))
    d = svd.tX$d
    V = svd.tX$v
    m = length(d[d > 1e-5])
    e = (d^2)[1:m]
    E = DIAG(e)
    V = V[,1:m]
    d = E^(1/2)%*%t(V)%*%Y
  }
  
  tau = 1
  # find a good initial value for rho
  
  #E_all = cbind(E2, E)
  #Z_all = as.matrix(bdiag(Z_old, Z))
  B = E2%*%Z_old
  fn = function(rho){
    A = diag(tau^2*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*tau^2*E
    M = solve(A)- solve(A)%*%B%*%SOLVE(t(B)%*%solve(A)%*%B)%*%t(B)%*%solve(A)
    q = tau*rho*d
    Mq = M%*%q
    z = Mq/norm(Mq, "2")
    return(t(z)%*%E%*%z + n*lambda - rho)
  }
  nleqslv.res = nleqslv(e[1]+n*lambda, fn, method = "Newton", global = "none", control = list(maxit = 150))
  rho = nleqslv.res$x
  if (nleqslv.res$termcd != 1 && verbose){
    print(paste0("Warning! The value is ", as.character(fn(rho))))
    print(nleqslv.res$termcd)
  }
  
  A = diag(tau^2*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*tau^2*E
  M = solve(A)- solve(A)%*%B%*%SOLVE(t(B)%*%solve(A)%*%B)%*%t(B)%*%solve(A)
  q = tau*rho*d
  Mq = M%*%q
  Z = Mq/norm(Mq, "2")
  
  E_all = cbind(E2, E)
  Z_all = as.matrix(bdiag(Z_old, Z))
  B = E_all%*%Z_all
  rho0 = rho
  
  while (ncol(Z) < om){
    nleqslv.res = nleqslv(rho0, fn, method = "Newton", global = "none", control = list(maxit = 150))
    rho = nleqslv.res$x
    if (nleqslv.res$termcd != 1 && verbose){
      print(paste0("Warning! The value is ", as.character(fn(rho))))
      print(nleqslv.res$termcd)
    }
    
    A = diag(tau^2*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*tau^2*E
    M = solve(A)- solve(A)%*%B%*%SOLVE(t(B)%*%solve(A)%*%B)%*%t(B)%*%solve(A)
    q = tau*rho*d
    Mq = M%*%q
    z = Mq/norm(Mq, "2")
    Z = cbind(Z, z)
    Z_all = as.matrix(bdiag(Z_old, Z))
    B = E_all%*%Z_all
    rho0 = rho
  }
  C = V%*%Z
  C = C[,0:min(ncol(C), om)]
  a = V%*%ginv(E)^(1/2)%*%Z
  return(list(C = as.matrix(C), a = a, V = as.matrix(V), Z = Z, E = E, D = D, U = as.matrix(U)))
}

initialize.UDVZ = function(X){
  svd.X = svd(X)
  d = svd.X$d
  V = svd.X$v
  U = svd.X$u
  m = rankMatrix(X)[1]
  
  if (m > 0){
    e = (d^2)[1:m]
    D = DIAG(d[1:m])
    E = DIAG(e)
    V = V[,1:m]
    U = U[,1:m]
  }else{
    U = matrix(,nrow=nrow(X), ncol=m)
    D = matrix(,nrow=m, ncol=m)
    E = matrix(,nrow=m, ncol=m)
    V = matrix(,nrow=ncol(X), ncol=m)
    e=0
  }
  Z = matrix(,nrow=m, ncol=0)
  return(list(U=as.matrix(U), D=D, V=as.matrix(V), Z=Z, E=E, e=e, m=m))
}

continuum.multigroup.iter = function(X.list, Y.list, lambda, gam, rankJ, rankA, maxiter = 1000, conv = 1e-7, 
                                     center.X = TRUE, scale.X = TRUE, center.Y = TRUE, scale.Y = TRUE, orthIndiv = F,
                                     I.initial = NULL, sd = 0){
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
  
  n_cumsum = c(0,cumsum(n))
  index.list = lapply(1:G, function(i) (n_cumsum[i]+1):n_cumsum[i+1] )
  
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
  
  if (gam == 0){
    rankJ = min(1, rankJ)
    if (rankJ){
      rankA = rep(0, G)
    }else{
      rankA = sapply(rankA, function(r) min(1, r))
    }
  }
  
  nrun = 0
  converged = F
  
  Y.homo = Y
  Y.homo.list = Y.list
  
  # C = matrix(0, p, rankJ)
  # Cind = lapply(1:G, function(g) matrix(0, p, rankA[g]))
  # P = matrix(0, p, p)
  C = matrix(rnorm(p*rankJ, 0, sd), p, rankJ)
  Cind = lapply(1:G, function(g) matrix(rnorm(p*rankA[g], 0, sd), p, rankA[g]))
  Cind_tot = do.call(cbind, Cind)
  P = Cind_tot%*%SOLVE(t(Cind_tot)%*%Cind_tot)%*%t(Cind_tot)
  
  if (is.null(I.initial)){
    # X.heter.list = lapply(1:G, function(g) matrix(0, nrow = nrow(X.list[[g]]), ncol = p))
    X.heter = X%*%P
    X.heter.list = lapply(index.list, function(ixs) as.matrix(X.heter[ixs,]))
    # X.heter.list = lapply(1:G, function(g) X.list[[g]]%*%Cind[[g]]%*%SOLVE(t(Cind[[g]])%*%Cind[[g]])%*%t(Cind[[g]]))
    X.homo.list = lapply(1:G, function(g) X.list[[g]] - X.heter.list[[g]])#%*%C%*%SOLVE(t(C)%*%C)%*%t(C))X.homo.list = lapply(1:G, function(g) X.list[[g]])#%*%C%*%SOLVE(t(C)%*%C)%*%t(C))
  }else{
    X.heter.list = I.initial
    X.homo.list = lapply(1:G, function(g) X.list[[g]] - X.heter.list[[g]])
  }
  
  R = X
  R[,] = 0
  
  r = Y
  r[,] = 0
  
  U = list()
  W = list()
  
  ct.homo = matrix(0, nrow = rankJ, ncol = 1)
  ct.heter = lapply(1:G, function(g) matrix(0, nrow = rankA[g], ncol = 1))
  
  UDVZ.heter.list = lapply(1:G, function(g) initialize.UDVZ(X.heter.list[[g]]))
  U.heter.list = lapply(1:G, function(g) UDVZ.heter.list[[g]]$U)
  Z.heter.list = lapply(1:G, function(g) UDVZ.heter.list[[g]]$Z)
  V.heter.list = lapply(1:G, function(g) UDVZ.heter.list[[g]]$V)
  D.heter.list = lapply(1:G, function(g) UDVZ.heter.list[[g]]$D)
  
  U.heter = as.matrix(do.call(bdiag, U.heter.list))
  Z.heter = as.matrix(do.call(bdiag, Z.heter.list))
  V.heter = do.call(rbind, V.heter.list)
  D.heter = as.matrix(do.call(bdiag, D.heter.list))
  
  # U.heter = matrix(,nrow=N, ncol=0)
  # D.heter = matrix(,nrow=0, ncol=0)
  # V.heter = matrix(,nrow=p, ncol=0)
  # Z.heter = matrix(,nrow=0, ncol=0)
  
  X.homo = do.call(rbind, X.homo.list)
  UDVZ.homo = initialize.UDVZ(X.homo)
  U.homo = UDVZ.homo$U
  Z.homo = UDVZ.homo$Z
  V.homo = UDVZ.homo$V
  D.homo = UDVZ.homo$D
  U.homo.list = lapply(index.list, function(ixs) as.matrix(U.homo[ixs,]))
  
  # U.homo.list = lapply(1:G, function(g) matrix(0, nrow = n[g], ncol = 0))
  # D.homo = matrix(,nrow=0, ncol=0)
  # V.homo = matrix(,nrow=p, ncol=0)
  # Z.homo = matrix(,nrow=0, ncol=0)
  
  while (nrun < maxiter & !converged){
    # initialization
    rlast = r
    Rlast = R
    
    ct.homo.last = ct.homo
    ct.heter.last = ct.heter
    
    X.homo.list = lapply(1:G, function(g) X.list[[g]] - X.heter.list[[g]])
    X.homo = do.call(rbind, X.homo.list)
    
    # joint
    if (rankJ){
      if (gam == 1e10){
        ml.homo = svd(X.homo%*%(diag(p) - P), nu = rankJ, nv = rankJ)
        C = ml.homo$v
      }else{
        # ml.homo = continuum.ridge.fix(X.homo%*%(diag(p) - P), Y.homo, G, lambda = lambda, gam = gam, om = rankJ)
        ml.homo = continuum(X.homo%*%(diag(p) - P), Y.homo, G, lambda = lambda, gam = gam, om = rankJ,
                            U_old=U.heter, D_old=D.heter, V_old=V.heter, Z_old=Z.heter)
        C = ml.homo$C
        
        U.homo = ml.homo$U
        Z.homo = ml.homo$Z
        V.homo = ml.homo$V
        D.homo = ml.homo$D
        U.homo.list = lapply(index.list, function(ixs) as.matrix(U.homo[ixs,]))
      }
    }
    U = list.append(U, C)
    
    XC = X.homo%*%C
    ct.homo = (diag(t(XC)%*%XC)^(gam-1))*(t(XC)%*%Y.homo)^2
    
    beta.C = C2beta(X.homo, Y.homo, C, lambda = lambda)$beta
    Yhat.homo.list = lapply(1:G, function(g) X.homo.list[[g]]%*%beta.C)
    Y.heter.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.homo.list[[g]])
    X.homo.list = lapply(1:G, function(g) X.homo.list[[g]]%*%C%*%SOLVE(t(C)%*%C)%*%t(C))
    X.heter.list = lapply(1:G, function(g) X.list[[g]] - X.homo.list[[g]])
    
    # individual
    temp = X.heter.list
    U.heter.list = list()
    Z.heter.list = list()
    V.heter.list = list()
    D.heter.list = list()
    for (g in 1:G){
      tempC = C
      # orthogonalization
      if (orthIndiv){
        if (nrun > 0){
          for (j in (1:G)[-g]){
            tempC = cbind(tempC, Cind[[j]])
          }
        }
      }
      temp[[g]] <- temp[[g]]%*%(diag(p) - tempC%*%SOLVE(t(tempC)%*%tempC)%*%t(tempC))
      # U.heter.list[[g]] = matrix(,nrow=n[g], ncol=0)
      # D.heter.list[[g]] = matrix(,nrow=0, ncol=0)
      # V.heter.list[[g]] = matrix(,nrow=p, ncol=0)
      # Z.heter.list[[g]] = matrix(,nrow=0, ncol=0)
#       if (rankA[g]){
      if (gam == 1e10){
        ml.heter = svd(temp[[g]], nu = rankA[g], nv = rankA[g])
        Cind[[g]] = ml.heter$v
      }else{
        # ml.heter = continuum.ridge.fix(temp[[g]], Y.heter.list[[g]], 1, lambda = lambda, gam = gam, om = rankA[g])
        ml.heter = continuum(temp[[g]], Y.heter.list[[g]], 1, lambda = lambda, gam = gam, om = rankA[g], U_old=U.homo.list[[g]], D_old=D.homo, V_old=V.homo, Z_old=Z.homo)
        Cind[[g]] = ml.heter$C
        
        U.heter.list[[g]] = ml.heter$U
        Z.heter.list[[g]] = ml.heter$Z
        V.heter.list[[g]] = ml.heter$V
        D.heter.list[[g]] = ml.heter$D
        # print(dim(ml.heter$U))
#       }
      }
    }
    
    U.heter = as.matrix(do.call(bdiag, U.heter.list))
    Z.heter = as.matrix(do.call(bdiag, Z.heter.list))
    V.heter = do.call(rbind, V.heter.list)
    D.heter = as.matrix(do.call(bdiag, D.heter.list))
    
    W = list.append(W, Cind)
    Cind_tot = do.call(cbind, Cind)
    P = Cind_tot%*%SOLVE(t(Cind_tot)%*%Cind_tot)%*%t(Cind_tot)
    
    XC.list = lapply(1:G, function(g) X.heter.list[[g]]%*%Cind[[g]])

    for (g in 1:G){
      print(t(X.homo.list[[g]]%*%C)%*%temp[[g]]%*%Cind[[g]])
    }
    ct.heter = lapply(1:G, function(g) 
      diag((t(XC.list[[g]])%*%XC.list[[g]])^(gam-1))*(t(XC.list[[g]])%*%Y.heter.list[[g]])^2)
    
    beta.Cind = lapply(1:G, function(g) C2beta(X.heter.list[[g]], Y.heter.list[[g]], Cind[[g]], lambda)$beta)
    Yhat.heter.list = lapply(1:G, function(g) X.heter.list[[g]]%*%beta.Cind[[g]])
    
    r.list = lapply(1:G, function(g) Y.heter.list[[g]] - Yhat.heter.list[[g]])
    r = do.call(rbind, r.list)
    
    Y.homo.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.heter.list[[g]])
    Y.homo = do.call(rbind, Y.homo.list)
    X.heter.list = lapply(1:G, function(g) X.heter.list[[g]]%*%Cind[[g]]%*%SOLVE(t(Cind[[g]])%*%Cind[[g]])%*%t(Cind[[g]]))
    
    R.list = lapply(1:G, function(g) X.list[[g]] - X.homo.list[[g]] - X.heter.list[[g]])
    R = do.call(rbind, R.list)
    
    if (gam > 1e5){
      if (!nrun%%10){
        print(norm(R, type = "f"))
      }
      if (norm(Rlast - R, type = "f") <= conv){
        converged <- T
      }
    }else{
      CT1 = (ct.homo - ct.homo.last)/ct.homo.last
      CT2 = lapply(1:G, function(g) as.vector((ct.heter[[g]]-ct.heter.last[[g]])/ct.heter.last[[g]]))
      CT = c(CT1, do.call(c, CT2))
      if (!nrun%%10){
        print(c(ct.homo, do.call(c, ct.heter)))
      }
      if (max(abs(CT), na.rm = T) <= conv){
        converged <- T
      }
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
  return(list(C = C, Cind = Cind, 
              U = U, W = W, 
              centerValues.X = centerValues.X, scaleValues.X = scaleValues.X, 
              centerValues.Y = centerValues.Y, scaleValues.Y = scaleValues.Y,
              intercept = intercept, beta.C = beta.C, beta.Cind = beta.Cind, 
              beta = beta.C0, beta_i = beta.Cind0,
              J = X.homo.list, I = X.heter.list, 
              R = R, r = r,
              converged = converged, nrun = nrun))
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

continuum.ridge.fix = function(X, Y, G, lambda, gam, om, vertical = TRUE, verbose = FALSE){
  #om: number of columns
  n = nrow(X)
  p = ncol(X)
  # scaleValues.Y = sqrt(mean(Y^2))
  # Y = sweep(matrix(Y), 2, scaleValues.Y, FUN = "/")
  
  if (vertical){
    # S = t(X)%*%X
    s = t(X)%*%Y
    
    svd.X = svd(X)
    d = svd.X$d
    V = svd.X$v
    m = rankMatrix(X)[1]
#    m = length(d[d > 1e-5])
    e = (d^2)[1:m]
    E = DIAG(e)
    V = V[,1:m]
    d = t(V)%*%s

    # e = eigen(S)$values
    # m = length(e[e > 1e-10])
    # e = e[1:m]
    # E = DIAG(e)
    # V = eigen(S)$vectors[, 1:m]
    # d = t(V)%*%s
  }else{
    # S = X%*%t(X)
    svd.tX = svd(t(X))
    d = svd.tX$d
    V = svd.tX$v
    m = length(d[d > 1e-5])
    e = (d^2)[1:m]
    E = DIAG(e)
    V = V[,1:m]
    d = E^(1/2)%*%t(V)%*%Y

    # e = eigen(S)$values
    # m = length(e[e > 1e-6])
    # e = e[1:m]
    # E = DIAG(e)
    # V = eigen(S)$vectors[, 1:m]
    # d = E^(1/2)%*%t(V)%*%Y
  }

#  tau = sqrt(2/(e[1]+e[m]))
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
  
  nleqslv.res = nleqslv(e[1]+n*lambda, fn, method = "Newton", global = "none", control = list(maxit = 150))
  # if (gam > 1){
  #   nleqslv.res = nleqslv(e[1]+n*lambda, fn, method = "Newton", global = "none", control = list(maxit = 150))
  # }else{
  #   nleqslv.res = nleqslv((e[1]+e[m])/2+n*lambda, fn, method = "Newton", global = "none", control = list(maxit = 150))
  # }
  
  rho = nleqslv.res$x
  if (nleqslv.res$termcd != 1 && verbose){
    print(paste0("Warning! The value is ", as.character(fn(rho))))
    #     print(fn(rho, gam, e, V, d))
    print(nleqslv.res$termcd)
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
  
#  while (ncol(Z) < om & rcond(t(B)%*%solve(A)%*%B) > 1e-10){
  while (ncol(Z) < om){
    # print(1)
    # print(rcond(t(B)%*%solve(A)%*%B))
    nleqslv.res = nleqslv(rho0, fn, method = "Newton", global = "none", control = list(maxit = 150))
    rho = nleqslv.res$x
    if (nleqslv.res$termcd != 1 && verbose){
      print(paste0("Warning! The value is ", as.character(fn(rho))))
      print(nleqslv.res$termcd)
    }
    
    A = diag(tau^2*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*tau^2*E
    
    # print(2)
    # print(rcond(t(B)%*%solve(A)%*%B))
    
    # if (rcond(t(B)%*%solve(A)%*%B) < 1e-10){
    #   break
    # }
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
  a = V%*%ginv(E)^(1/2)%*%Z
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
  a = V%*%ginv(E)^(1/2)%*%Z
  return(list(C = as.matrix(C), a = a, V = as.matrix(V), Z = Z, E = E))
}



