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



cv.continnum.iter = function(X.list, Y.list, lambda = 0, parameter.set, nfolds = 10, maxiter = 100,
                             center.X = TRUE, scale.X = TRUE, center.Y = TRUE, scale.Y = TRUE, orthIndiv = T, 
                             plot = F, criteria = c("min", "1se"), sd = 0){
  G = length(X.list)
  flds.list = lapply(1:G, function(g) createFolds(Y.list[[g]], k = nfolds, list = TRUE, returnTrain = FALSE))
  MSE.list = list()
  for (k in 1:nfolds){
    X.train.list = lapply(1:G, function(g) X.list[[g]][unlist(flds.list[[g]][-k]), ])
    X.val.list = lapply(1:G, function(g) X.list[[g]][unlist(flds.list[[g]][k]), ])
    Y.train.list = lapply(1:G, function(g) matrix(Y.list[[g]][unlist(flds.list[[g]][-k])]))
    Y.val.list = lapply(1:G, function(g) matrix(Y.list[[g]][unlist(flds.list[[g]][k])]))
    Y.val = do.call(rbind, Y.val.list)
    
    ml.list = lapply(parameter.set, function(parameter) 
      continuum.multigroup.iter(X.train.list, Y.train.list, lambda = lambda, maxiter = maxiter,     
                                gam = parameter$gam, rankJ = parameter$rankJ, rankA = parameter$rankA, 
                                center.X = center.X, scale.X = scale.X, center.Y = center.Y, scale.Y = scale.Y, orthIndiv = orthIndiv,
                                sd = sd))
    Yhat.list = lapply(ml.list, function(ml) 
      do.call(rbind, lapply(1:G, function(g) as.numeric(ml$intercept[[g]]) + X.val.list[[g]]%*%ml$beta.C[[g]] + X.val.list[[g]]%*%ml$beta.Cind[[g]])))
    MSE.list[[k]] = sapply(Yhat.list, function(Yhat) mean((Y.val - Yhat)^2))
  }
  MSE = do.call(rbind, MSE.list)
  rMSE = apply(MSE, 2, function(x) mean(sqrt(x)))
  if (plot){
    plot(rMSE)
  }
  if (criteria == "1se"){
    absBest = min(rMSE)
    MSEsd = apply(MSE, 2, function(x) sd(sqrt(x)))/sqrt(nfolds)
    # ix = min(which((rMSE - MSEsd) < absBest))
    # parameter = parameter.set[[ix]]
    
    rMSE.mtx = matrix(rMSE, ncol = L+1, byrow = T)
    absBest.ix = which.min(rMSE.mtx)
    #    print(absBest.ix)
    #    absBest.row.ix = (absBest.ix-1)%%(nrow(rMSE.mtx))+1
    col.ix = ceiling(absBest.ix/nrow(rMSE.mtx))
    #    print(col.ix)
    row.ix = min(which((rMSE.mtx[, col.ix] - MSEsd[absBest.ix]) < absBest))
    #    print(row.ix)
    ix = col.ix + (row.ix-1)*(L+1)
    parameter = parameter.set[[ix]]
  }
  if (criteria == "min"){
    ix = which.min(rMSE)
    parameter = parameter.set[[ix]]
  }
  return(list(rMSE = rMSE, MSE = MSE, ix = ix, parameter = parameter))
}
