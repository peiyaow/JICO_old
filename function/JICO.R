library(quadprog)
library(nleqslv)
library(Matrix)
library(rlist)
library(MASS)

parameter.set.G_2 = function(maxrankA, maxrankJ, gamma){
  # generate set of hyperparameters when there are G = 2 groups. 
  
  # maxrankA: the maximum rank for individual component
  # maxrankJ: the maximum rank for joint component
  # gamma: the fixed gamma parameter
  
  parameter.set <- list()
  for(rankA1 in 0:maxrankA)
    for(rankA2 in 0:maxrankA)
      for (rankJ in 0:maxrankJ){
        parameter.set = list.append(parameter.set, list(rankA = c(rankA1, rankA2), rankJ = rankJ, gam = gamma))
      }
  return(parameter.set)
}

parameter.set.G_3 = function(maxrankA, maxrankJ, gamma){
  # generate set of hyperparameters when there are G = 3 groups. 
  
  # maxrankA: the maximum rank for individual component
  # maxrankJ: the maximum rank for joint component
  # gamma: the fixed gamma parameter
  
  parameter.set <- list()
  for(rankA1 in 0:maxrankA)
    for(rankA2 in 0:maxrankA)
      for (rankA3 in 0:maxrankA){
        for (rankJ in 0:maxrankJ){
          parameter.set = list.append(parameter.set, list(rankA = c(rankA1, rankA2, rankA3), rankJ = rankJ, gam = gamma))
        }
      }
  return(parameter.set)
}

parameter.set.rankA_eq = function(G, maxrankA, maxrankJ, gamma.list){
  # generate set of hyperparameters when the individual ranks are the same
  
  # G: number of groups
  # maxrankA: the maximum rank for individual component
  # maxrankJ: the maximum rank for joint component
  # gamma.list: the list of candidate gamma parameters to be tuned
  
  parameter.set <- list() 
  for (gam in gamma.list)
    for(rankA in 0:maxrankA)
      for (rankJ in 0:maxrankJ){
        parameter.set = list.append(parameter.set, list(rankA = rep(rankA, G), rankJ = rankJ, gam = gam))
      }
  return(parameter.set)
}

SOLVE = function(x){
  # algorithm to inverse a matrix x
  
  if (sum(dim(x))){
    return(ginv(x))
  }else{
    return(x)
  }
}

DIAG = function(e){
  # given diagonal element, return a matrix 
  
  if (length(e) > 1){
    return(diag(e))
  }else{
    return(matrix(e))
  }
}

C2beta = function(X, Y, C, lambda){
  # Compute the coefficients from the continuum regression (CR) algorithm
  
  # C: the weight matrix computed from CR algorithm
  # X, Y: the input data pairs for CR algorithm
  # lambda: regularization parameter if L2 penalization is used for CR. JICO uses zero as default
  
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

continuum = function(X, Y, lambda, gam, om, 
                     U_old=matrix(,nrow=nrow(X), ncol=0), D_old=matrix(,nrow=0, ncol=0), V_old=matrix(,nrow=0, ncol=0), Z_old=matrix(,nrow=0, ncol=0), 
                     verbose = FALSE){
  # Implement the main CR algorithm used in JICO
  
  # X, Y: the input data pairs for continuum regression algorithm
  # lambda: regularization parameter if L2 penalization is used for CR. JICO uses 0 as default
  # gam: gamma parameter in CR algorithm. Set gam=0 for OLS model, gam=0.5 for PLS model, gam >= 1e10 for PCR model
  # om: the desired number of weight vectors to achieve in the CR algorithm, i.e. the predefined rank of joint or individual componenet
  # U_old, D_old, V_old, Z_old: the given inputs U, D, V, Z from the previous JICO iteration. Detailed notation can be found in JICO paper supplement Web Appendix B.
  # verbose: if it's desired to print out intermediate outputs

  # ----------- initialization ----------- #
  n = nrow(X)
  p = ncol(X)
  
  s = t(X)%*%Y
  UDVZ = initialize.UDVZ(X)
  D = UDVZ$D
  E = UDVZ$E
  V = UDVZ$V
  U = UDVZ$U
  e = UDVZ$e
  m = UDVZ$m
  d = t(V)%*%s
  E2 = D%*%t(U)%*%U_old%*%D_old
  
  # default to 1
  tau = 1
  
  # ----------- find a good initial value for rho ----------- #
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
  
  # ----------- iteration on the first step ----------- #
  A = diag(tau^2*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*tau^2*E
  M = solve(A)- solve(A)%*%B%*%SOLVE(t(B)%*%solve(A)%*%B)%*%t(B)%*%solve(A)
  q = tau*rho*d
  Mq = M%*%q
  Z = Mq/norm(Mq, "2")
  
  E_all = cbind(E2, E)
  Z_all = as.matrix(bdiag(Z_old, Z))
  B = E_all%*%Z_all
  rho0 = rho
  
  # ----------- iteration on the following steps if om > 1 ----------- #
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
  
  # ----------- compute final results ----------- 
  C = V%*%Z
  C = C[, 0:min(ncol(C), om)]
  a = V%*%ginv(E)^(1/2)%*%Z
  
  return(list(C = as.matrix(C), a = a, V = as.matrix(V), Z = Z, E = E, D = D, U = as.matrix(U)))
}

initialize.UDVZ = function(X){
  # helper function to compute the SVD results from a given matrix X 
  
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


continuum.multigroup.iter = function(X.list, Y.list, lambda = 0, gam, rankJ, rankA, maxiter = 1000, conv = 1e-7, 
                                     center.X = TRUE, scale.X = TRUE, center.Y = TRUE, scale.Y = TRUE, orthIndiv = F,
                                     I.initial = NULL, sd = 0){
  # main function to iteratively solve the JICO algorithm
  
  # X.list: the list of observational data in each group
  # Y.list: the list of response in each group
  # lambda: regularization parameter if L2 penalization is used for CR. JICO uses 0 as default
  # gam: gamma parameter in CR algorithm. Set gam=0 for OLS model, gam=0.5 for PLS model, gam >= 1e10 for PCR model
  # rankJ: the rank for the joint component
  # rankA: the ranks for individual components
  # maxiter: the maximum of iterations to conduct before convergence
  # conv: the tolerance level for convergence
  # center.X: if centralization is to be performed on X
  # scale.X: if scaling is performed on X
  # center.Y: if centralization is to be performed on Y
  # scale.Y: if scaling is performed on Y
  # orthIndiv: if we impose the orthogonality constraint on individual components
  # I.initial: the initial value for individual components
  # sd: standard deviation used to generate random initial values for individual weight vectors
  
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
  
  # when gam == 0, it is equivalent to OLS. Either rankJ or rankA needs to be 0 and at most to be 1
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
  
  C = matrix(rnorm(p*rankJ, 0, sd), p, rankJ)
  Cind = lapply(1:G, function(g) matrix(rnorm(p*rankA[g], 0, sd), p, rankA[g]))
  Cind_tot = do.call(cbind, Cind)
  P = Cind_tot%*%SOLVE(t(Cind_tot)%*%Cind_tot)%*%t(Cind_tot)
  
  if (is.null(I.initial)){
    X.heter = X%*%P
    X.heter.list = lapply(index.list, function(ixs) as.matrix(X.heter[ixs,]))
    X.homo.list = lapply(1:G, function(g) X.list[[g]] - X.heter.list[[g]])
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
  
  X.homo = do.call(rbind, X.homo.list)
  UDVZ.homo = initialize.UDVZ(X.homo)
  U.homo = UDVZ.homo$U
  Z.homo = UDVZ.homo$Z
  V.homo = UDVZ.homo$V
  D.homo = UDVZ.homo$D
  U.homo.list = lapply(index.list, function(ixs) as.matrix(U.homo[ixs,]))
  
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
        ml.homo = continuum(X.homo%*%(diag(p) - P), Y.homo, lambda = lambda, gam = gam, om = rankJ,
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
      
      if (gam == 1e10){
        ml.heter = svd(temp[[g]], nu = rankA[g], nv = rankA[g])
        Cind[[g]] = ml.heter$v
      }else{
        ml.heter = continuum(temp[[g]], Y.heter.list[[g]], lambda = lambda, gam = gam, om = rankA[g], U_old=U.homo.list[[g]], D_old=D.homo, V_old=V.homo, Z_old=Z.homo)
        Cind[[g]] = ml.heter$C
        
        U.heter.list[[g]] = ml.heter$U
        Z.heter.list[[g]] = ml.heter$Z
        V.heter.list[[g]] = ml.heter$V
        D.heter.list[[g]] = ml.heter$D
        
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
    
    ct.heter = lapply(1:G, function(g) 
      diag((t(XC.list[[g]])%*%XC.list[[g]])^(gam-1))*(t(XC.list[[g]])%*%Y.heter.list[[g]])^2)
    
    beta.Cind = lapply(1:G, function(g) C2beta(X.heter.list[[g]], Y.heter.list[[g]], Cind[[g]], lambda)$beta)
    Yhat.heter.list = lapply(1:G, function(g) X.heter.list[[g]]%*%beta.Cind[[g]])
    
    r.list = lapply(1:G, function(g) Y.heter.list[[g]] - Yhat.heter.list[[g]])
    r = do.call(rbind, r.list)
    
    Y.homo.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.heter.list[[g]])
    Y.homo = do.call(rbind, Y.homo.list)
    X.heter.list = lapply(1:G, function(g) X.heter.list[[g]]%*%Cind[[g]]%*%SOLVE(t(Cind[[g]])%*%Cind[[g]])%*%t(Cind[[g]]))
    
    
    # compute residuals
    R.list = lapply(1:G, function(g) X.list[[g]] - X.homo.list[[g]] - X.heter.list[[g]])
    R = do.call(rbind, R.list)
    
    if (gam > 1e5){
      if (norm(Rlast - R, type = "f") <= conv){
        converged <- T
      }
    }else{
      CT1 = (ct.homo - ct.homo.last)/ct.homo.last
      CT2 = lapply(1:G, function(g) as.vector((ct.heter[[g]]-ct.heter.last[[g]])/ct.heter.last[[g]]))
      CT = c(CT1, do.call(c, CT2))
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

createFolds <- function(strat_id, k) {
  # function to create folds for cross validation
  
  if(k > length(strat_id)) {
    k <- length(strat_id)
  }	
  perm <- sample(length(strat_id))
  strat_id <- factor(strat_id, levels=sample(unique(strat_id)))
  
  strat_order <- order(strat_id[perm])
  
  num_item_per_fold_ceil <- ceiling(length(strat_id) / k)
  
  fold_ids <- rep(seq_len(k), times= num_item_per_fold_ceil)
  fold_ids <- fold_ids[seq_along(strat_id)]
  
  folds <- split(perm[strat_order], fold_ids)
  names(folds) <- paste0("Fold", seq_len(k))	
  return(folds)
}

cv.continnum.iter = function(X.list, Y.list, lambda = 0, parameter.set, nfolds = 10, maxiter = 100,
                             center.X = TRUE, scale.X = TRUE, center.Y = TRUE, scale.Y = TRUE, orthIndiv = T, 
                             plot = F, criteria = c("min", "1se"), sd = 0){
  # function to perform cross validation to select the best parameter
  
  # X.list: the list of observational data in each group
  # Y.list: the list of response in each group
  # lambda: regularization parameter if L2 penalization is used for CR. JICO uses 0 as default
  # parameter.set: the set of parameters to be tuned on
  # nfolds: number of folds to create to perform CV
  # maxiter: the maximum of iterations to conduct for JICO before convergence
  # center.X: if centralization is to be performed on X
  # scale.X: if scaling is performed on X
  # center.Y: if centralization is to be performed on Y
  # scale.Y: if scaling is performed on Y
  # orthIndiv: if we impose the orthogonality constraint on individual components
  # plot: if we want to plot the rMSE vs different parameters
  # criteria: criteria for selecting the best parameter. 
  #           use "min" to choose the parameter giving the best performance, 
  #           use "1se" to choose the simplest model that gives performance within 1se from the best one
  # sd: standard deviation used to generate random initial values for individual weight vectors
  
  G = length(X.list)
  flds.list = lapply(1:G, function(g) createFolds(Y.list[[g]], k = nfolds))
  MSE.list = list()
  for (k in 1:nfolds){
    X.train.list = lapply(1:G, function(g) X.list[[g]][unlist(flds.list[[g]][-k]), ])
    X.val.list = lapply(1:G, function(g) X.list[[g]][unlist(flds.list[[g]][k]), ])
    Y.train.list = lapply(1:G, function(g) matrix(Y.list[[g]][unlist(flds.list[[g]][-k])]))
    Y.val.list = lapply(1:G, function(g) matrix(Y.list[[g]][unlist(flds.list[[g]][k])]))
    Y.val = do.call(rbind, Y.val.list)
    
    ml.list = lapply(parameter.set, function(parameter) 
      continuum.multigroup.iter(X.train.list, Y.train.list, maxiter = maxiter,     
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
    
    gam_list = c()
    for (para in parameter.set){
      gam_list = c(gam_list, para$gam)
    }
    num_gam = length(unique(gam_list))
    
    rMSE.mtx = matrix(rMSE, ncol = num_gam) # num of cols = number of gams 
    absBest.ix = which.min(rMSE.mtx)
    col.ix = ceiling(absBest.ix/nrow(rMSE.mtx))
    row.ix = min(which((rMSE.mtx[, col.ix] - MSEsd[absBest.ix]) < absBest))
    ix = col.ix + (row.ix-1)*(num_gam)
    parameter = parameter.set[[ix]]
  }
  if (criteria == "min"){
    ix = which.min(rMSE)
    parameter = parameter.set[[ix]]
  }
  return(list(rMSE = rMSE, MSE = MSE, ix = ix, parameter = parameter))
}
