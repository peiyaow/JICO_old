library(caret)
continuum.step1 = function(X.list, Y.list, lambda = 0, gam = 1, m, center.X = TRUE, scale.X = TRUE, center.Y = TRUE, scale.Y = TRUE, tune = TRUE){
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
  
  YOrigin.list = Y.list
  XOrigin.list = X.list
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
  
  ml = continuum.ridge.fix(X = X, Y = Y, lambda = lambda, gam = gam, om = m)
  C = ml$C
  Y.heter.list = NULL
  X.heter.list = NULL
  YOriginRes.list = NULL
  if (tune){
    beta.list = lapply(0:m, function(mm) C2beta(X = X, Y = Y, C = C[,0:mm], lambda = lambda)$beta)
    
    betaOrigin.list = list()
    intercept.list = list()
    for (i in 1:(m+1)){
      betaOrigin.list[[i]] = list()
      intercept.list[[i]] = list()
      for (g in 1:G){
        betaOrigin.list[[i]][[g]] = beta.list[[i]]/scaleValues.X[[g]]*scaleValues.Y[[g]]
        intercept.list[[i]][[g]] = centerValues.Y[[g]] - t(betaOrigin.list[[i]][[g]])%*%centerValues.X[[g]]
      }
    }
  }else{ # only compute one beta
    beta.C = C2beta(X = X, Y = Y, C = C, lambda = lambda)$beta
    Yhat.homo.list = lapply(1:G, function(g) X.list[[g]]%*%beta.C)
    Y.heter.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.homo.list[[g]])
    X.homo.list = lapply(1:G, function(g) X.list[[g]]%*%C%*%SOLVE(t(C)%*%C)%*%t(C))
    X.heter.list = lapply(1:G, function(g) X.list[[g]] - X.homo.list[[g]])
    
    betaOrigin.list = list()
    intercept.list = list()
    for (g in 1:G){
      betaOrigin.list[[g]] = beta.C/scaleValues.X[[g]]*scaleValues.Y[[g]]
      intercept.list[[g]] = centerValues.Y[[g]] - t(betaOrigin.list[[g]])%*%centerValues.X[[g]]
    }
    
    YOriginRes.list = lapply(1:G, function(g) YOrigin.list[[g]] - as.numeric(intercept.list[[g]]) - XOrigin.list[[g]]%*%betaOrigin.list[[g]])
  }
  return(list(beta = betaOrigin.list, intercept = intercept.list, C = C, 
              X.heter = X.heter.list, Y.heter = Y.heter.list, 
              XOrigin = XOrigin.list, YOriginRes = YOriginRes.list,
              centerValues.X = centerValues.X, centerValues.Y = centerValues.Y, 
              scaleValues.X = scaleValues.X, scaleValues.Y = scaleValues.Y))
}

cv.continuum.step1 = function(X.list, Y.list, lambda = 0, gam = 1, nfolds = 10, m = 5, 
                              center.X = TRUE, scale.X = TRUE, center.Y = TRUE, scale.Y = TRUE, plot = F, criteria = c("min", "1se")){
  G = length(X.list)
#  set.seed(111)
  flds.list = lapply(1:G, function(g) createFolds(Y.list[[g]], k = nfolds, list = TRUE, returnTrain = FALSE))
  MSE.list = list()
  for (k in 1:nfolds){
    X.train.list = lapply(1:G, function(g) X.list[[g]][unlist(flds.list[[g]][-k]), ])
    X.val.list = lapply(1:G, function(g) X.list[[g]][unlist(flds.list[[g]][k]), ])
    Y.train.list = lapply(1:G, function(g) matrix(Y.list[[g]][unlist(flds.list[[g]][-k])]))
    Y.val.list = lapply(1:G, function(g) matrix(Y.list[[g]][unlist(flds.list[[g]][k])]))
    Y.val = do.call(rbind, Y.val.list)
    
    ml = continuum.step1(X.train.list, Y.train.list, lambda = lambda, gam = gam, m = m, 
                         center.X = center.X, scale.X = scale.X, center.Y = center.Y, scale.Y = scale.Y)
    Yhat.list = lapply(1:(m+1), function(mm) 
      do.call(rbind, lapply(1:G, function(g) as.numeric(ml$intercept[[mm]][[g]]) + X.val.list[[g]]%*%ml$beta[[mm]][[g]])))
    MSE.list[[k]] = sapply(1:(m+1), function(mm) mean((Y.val - Yhat.list[[mm]])^2))
  }
  
  MSE = do.call(rbind, MSE.list)
  rMSE = apply(MSE, 2, function(x) mean(sqrt(x)))
  if (plot){
    plot(rMSE)
  }
  
  if (criteria =="1se"){
    absBest = min(rMSE)
    MSEsd = apply(MSE, 2, function(x) sd(sqrt(x)))/sqrt(nfolds)
    rankJ = min(which((rMSE - MSEsd) < absBest)) - 1
  }
  if (criteria =="min"){
    rankJ = which.min(rMSE) - 1
  }
  
  if (criteria == "scree"){
    rMSEratio = rMSE[-(m+1)]/rMSE[-1]
    rankJ = which.max(rMSEratio)
    if (rMSEratio[rankJ] < 1){
      rankJ = 0
    }
    if (plot){
      plot(rMSEratio)
    }
  }
  ml = continuum.step1(X.list, Y.list, lambda = lambda, gam = gam, m = rankJ, 
                       center.X = center.X, scale.X = scale.X, center.Y = center.Y, scale.Y = scale.Y, tune = F)
  return(list(rMSE = rMSE, rankJ = rankJ, ml = ml, MSE = MSE))
}

continuum.step2 = function(Xres, Yres, C, lambda = 0, gam = 1, m = 5, centerValues.X, scaleValues.X, scaleValues.Y, tune = T){
  p = ncol(Xres)
  ml = continuum.ridge.fix(X = Xres%*%(diag(p) - C%*%SOLVE(t(C)%*%C)%*%t(C)), Y = Yres, lambda = lambda, gam = gam, om = m)
  C = ml$C
  
  if (tune){
    beta.list = lapply(0:m, function(mm) C2beta(X = Xres, Y = Yres, C = C[,0:mm], lambda = lambda)$beta)
    betaOrigin = list()
    intercept = list()
    for (i in 1:(m+1)){
      betaOrigin[[i]] = beta.list[[i]]/scaleValues.X*scaleValues.Y
      intercept[[i]] = -t(betaOrigin[[i]])%*%centerValues.X
    }
  }else{
    beta.C = C2beta(X = Xres, Y = Yres, C = C, lambda = lambda)$beta
    betaOrigin = beta.C/scaleValues.X*scaleValues.Y
    intercept = -t(betaOrigin)%*%centerValues.X
  }
  return(list(beta = betaOrigin, intercept = intercept))
}

cv.continuum.step2 = function(Xres, Yres, C, XOrigin, YOriginRes, lambda = 0, gam = 1, nfolds = 10, m = 5, 
                              centerValues.X, scaleValues.X, scaleValues.Y, plot = F, criteria = c("min", "1se")){
  flds = createFolds(Yres, k = nfolds, list = TRUE, returnTrain = FALSE)
  MSE.list = list()
  for (k in 1:nfolds){
    Xres.train = Xres[unlist(flds[-k]), ]
    Xres.val = Xres[unlist(flds[k]), ]
    Yres.train = matrix(Yres[unlist(flds[-k])])
    Yres.val = matrix(Yres[unlist(flds[k])])
    
    XOrigin.val = XOrigin[unlist(flds[k]), ]
    YOriginRes.val = matrix(YOriginRes[unlist(flds[k])])
    ml = continuum.step2(Xres.train, Yres.train, C, lambda = lambda, gam = gam, m = m, centerValues.X, scaleValues.X, scaleValues.Y)
    Yhat.list = lapply(1:(m+1), function(mm) as.numeric(ml$intercept[[mm]]) + XOrigin.val%*%ml$beta[[mm]])
    MSE.list[[k]] = sapply(1:(m+1), function(mm) mean((YOriginRes.val - Yhat.list[[mm]])^2))
  }
  
  MSE = do.call(rbind, MSE.list)
  rMSE = apply(MSE, 2, function(x) mean(sqrt(x)))
  if (plot){
    plot(rMSE)
  }
  if (criteria == "1se"){
    absBest = min(rMSE)
    MSEsd = apply(MSE, 2, function(x) sd(sqrt(x)))/sqrt(nfolds)
    rankA = min(which((rMSE - MSEsd) < absBest)) - 1
  }
  if (criteria == "min"){
    rankA = which.min(rMSE) - 1
  }
  
  if (criteria == "scree"){
    if (rMSE[1]/rMSE[2] < 1){
      rankA = 0
    }else{
      rankA = which.max(rMSE[-(m+1)]/rMSE[-1])
    }
  }
  ml = continuum.step2(Xres, Yres, C, lambda = lambda, gam = gam, m = rankA, centerValues.X, scaleValues.X, scaleValues.Y, tune = F)
  
  return(list(rMSE = rMSE, rankA = rankA, ml = ml))
}


cv.continnum.2step.separate = function(X.list, Y.list, lambda = 0, gam = 1, nfolds = 10, m1 = 10, m2 = 5,
                              center.X = TRUE, scale.X = TRUE, center.Y = TRUE, scale.Y = TRUE, 
                              CT1 = "min", CT2 = "min", plot = F){
  # fix gam: separate tune rankJ and rankA
  G = length(X.list)
  ml.step1 = cv.continuum.step1(X.list, Y.list, lambda = lambda, gam = gam, nfolds = nfolds, m = m1, 
                                center.X = center.X, scale.X = scale.X, center.Y = center.Y, scale.Y = scale.Y,
                                plot = plot, criteria = CT1)
  rankJ = ml.step1$rankJ
  
  ml.step2 = lapply(1:G, function(g) cv.continuum.step2(ml.step1$ml$X.heter[[g]], ml.step1$ml$Y.heter[[g]], ml.step1$ml$C,
                                                        ml.step1$ml$XOrigin[[g]], ml.step1$ml$YOriginRes[[g]],
                                                        lambda = lambda, gam = gam, nfolds = nfolds, m = m2,
                                                        ml.step1$ml$centerValues.X[[g]], ml.step1$ml$scaleValues.X[[g]],
                                                        ml.step1$ml$scaleValues.Y[[g]],
                                                        plot = plot, criteria = CT2))
  rankA = sapply(1:G, function(g) ml.step2[[g]]$rankA)
  
  intercept = lapply(1:G, function(g) ml.step1$ml$intercept[[g]] + ml.step2[[g]]$ml$intercept)
  beta.C = lapply(1:G, function(g) ml.step1$ml$beta[[g]])
  beta.Cind = lapply(1:G, function(g) ml.step2[[g]]$ml$beta)
  return(list(intercept = intercept, beta.C = beta.C, beta.Cind = beta.Cind, ml.step1 = ml.step1, ml.step2 = ml.step2, parameter = list(gam = gam, rankJ = rankJ, rankA = rankA)))
}

cv.continnum.2step = function(X.list, Y.list, lambda = 0, parameter.set, nfolds = 10,
                              center.X = TRUE, scale.X = TRUE, center.Y = TRUE, scale.Y = TRUE, plot = F, criteria = c("min", "1se")){
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
      continuum.2step(X.train.list, Y.train.list, lambda = lambda,      
                      gam = parameter$gam, rankJ = parameter$rankJ, rankA = parameter$rankA, 
                      center.X = center.X, scale.X = scale.X, center.Y = center.Y, scale.Y = scale.Y))
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


continuum.2step = function(X.list, Y.list, lambda = 0, gam = 1, rankJ, rankA, 
                           center.X = TRUE, scale.X = TRUE, center.Y = TRUE, scale.Y = TRUE){
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
  
  if (gam == 0){
    rankJ = min(1, rankJ)
    if (rankJ){
      rankA = rep(0, G)
    }else{
      rankA = sapply(rankA, function(r) min(1, r))
    }
  }
  
  ml.homo = continuum.ridge.fix(X = X, Y = Y, lambda = lambda, gam = gam, om = rankJ)
  C = ml.homo$C
  beta.C = C2beta(X, Y, C, lambda = lambda)$beta
  Yhat.homo.list = lapply(1:G, function(g) X.list[[g]]%*%beta.C)
  Y.heter.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.homo.list[[g]])
  X.homo.list = lapply(1:G, function(g) X.list[[g]]%*%C%*%SOLVE(t(C)%*%C)%*%t(C))
  X.heter.list = lapply(1:G, function(g) X.list[[g]] - X.homo.list[[g]])
  ml.heter = lapply(1:G, function(g) 
    continuum.ridge.fix(X = X.heter.list[[g]]%*%(diag(p) - C%*%SOLVE(t(C)%*%C)%*%t(C)), Y = Y.heter.list[[g]], 
                        lambda = lambda, gam = gam, om = rankA[g]))
  Cind = lapply(1:G, function(g) ml.heter[[g]]$C)
  beta.Cind = lapply(1:G, function(g) C2beta(X.heter.list[[g]], Y.heter.list[[g]], Cind[[g]], lambda)$beta)
  Yhat.heter.list = lapply(1:G, function(g) X.heter.list[[g]]%*%beta.Cind[[g]])
  
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
  return(list(C = C, Cind = Cind, intercept = intercept, beta.C = beta.C, beta.Cind = beta.Cind))
}



cv.continnum.iter = function(X.list, Y.list, lambda = 0, parameter.set, nfolds = 10, maxiter = 100,
                             center.X = TRUE, scale.X = TRUE, center.Y = TRUE, scale.Y = TRUE, orthIndiv = T, 
                             plot = F, criteria = c("min", "1se")){
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
                      center.X = center.X, scale.X = scale.X, center.Y = center.Y, scale.Y = scale.Y, orthIndiv = orthIndiv))
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


# continuum = function(X, Y, lambda = 0, gam = 1, m, center.X = TRUE, scale.X = TRUE, center.Y = TRUE, scale.Y = TRUE, tune = TRUE){
#   n = nrow(X)
#   p = ncol(X)
#   
#   if (center.X){
#     centerValues.X = apply(X, 2, mean)
#   }else{
#     centerValues.X = rep(0, p)
#   }
#   if (scale.X){
#     scaleValues.X = norm(X, type = "f") 
#   }else{
#     scaleValues.X = 1
#   }
#   X = sweep(X, 2, centerValues.X)
#   X = X/scaleValues.X
#   
#   if (center.Y){
#     centerValues.Y = mean(Y)
#   }else{
#     centerValues.Y = 0
#   }
#   if (scale.Y){
#     scaleValues.Y = norm(Y, type = "f") #*sqrt(N)
#   }else{
#     scaleValues.Y = 1
#   }
#   
#   Y = sweep(matrix(Y), 2, centerValues.Y)
#   Y = sweep(matrix(Y), 2, scaleValues.Y, FUN = "/")
#   
#   ml = continuum.ridge.fix(X = X, Y = Y, lambda = lambda, gam = gam, om = m)
#   C = ml$C
#   if (tune){
#     beta.list = lapply(0:m, function(mm) C2beta(X = X, Y = Y, C = C[,0:mm], lambda = lambda)$beta)
#     
#     betaOrigin = list()
#     intercept = list()
#     for (i in 1:(m+1)){
#       betaOrigin[[i]] = beta.list[[i]]/scaleValues.X*scaleValues.Y
#       intercept[[i]] = centerValues.Y - t(betaOrigin[[i]])%*%centerValues.X
#     }
#   }else{ # only compute one beta
#     beta.C = C2beta(X = X, Y = Y, C = C, lambda = lambda)$beta
#     # Yhat.homo.list = lapply(1:G, function(g) X.list[[g]]%*%beta.C)
#     # Y.heter.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.homo.list[[g]])
#     # X.homo.list = lapply(1:G, function(g) X.list[[g]]%*%C%*%SOLVE(t(C)%*%C)%*%t(C))
#     # X.heter.list = lapply(1:G, function(g) X.list[[g]] - X.homo.list[[g]])
#     betaOrigin = beta.C/scaleValues.X*scaleValues.Y
#     intercept = centerValues.Y - t(betaOrigin)%*%centerValues.X
#   }
#   return(list(beta = betaOrigin, intercept = intercept, C = C, 
#               centerValues.X = centerValues.X, centerValues.Y = centerValues.Y, 
#               scaleValues.X = scaleValues.X, scaleValues.Y = scaleValues.Y))
# }


cv.continuum = function(X, Y, lambda = 0, gam = 1, nfolds = 10, m = 5, plot = F, criteria = "min",
                              center.X = TRUE, scale.X = TRUE, center.Y = TRUE, scale.Y = TRUE){
#  set.seed(111)
  flds = createFolds(Y, k = nfolds, list = TRUE, returnTrain = FALSE)
  MSE.list = list()
  for (k in 1:nfolds){
    X.train = X[unlist(flds[-k]), ]
    X.val = X[unlist(flds[k]), ]
    Y.train = matrix(Y[unlist(flds[-k])])
    Y.val = matrix(Y[unlist(flds[k])])
    
    ml = continuum(X.train, Y.train, lambda = lambda, gam = gam, m = m, 
                   center.X = center.X, scale.X = scale.X, center.Y = center.Y, scale.Y = scale.Y)
    Yhat.list = lapply(1:(m+1), function(mm) as.numeric(ml$intercept[[mm]]) + X.val%*%ml$beta[[mm]])
    MSE.list[[k]] = sapply(1:(m+1), function(mm) mean((Y.val - Yhat.list[[mm]])^2))
  }
  
  MSE = do.call(rbind, MSE.list)
  rMSE = apply(MSE, 2, function(x) mean(sqrt(x)))
  if (plot){
    plot(rMSE)
  }
  if (criteria == "1se"){
    absBest = min(rMSE)
    MSEsd = apply(MSE, 2, function(x) sd(sqrt(x)))/sqrt(nfolds)
    rankA = min(which((rMSE - MSEsd) < absBest)) - 1
  }
  if (criteria == "min"){
    rankA = which.min(rMSE) - 1
  }
  
  if (criteria == "scree"){
    if (rMSE[1]/rMSE[2] < 1){
      rankA = 0
    }else{
      rankA = which.max(rMSE[-(m+1)]/rMSE[-1])
    }
  }
  ml = continuum(X, Y, lambda = lambda, gam = gam, m = rankA, 
                 center.X = center.X, scale.X = scale.X, center.Y = center.Y, scale.Y = scale.Y,
                 tune = F)
  
  return(list(rMSE = rMSE, rankA = rankA, ml = ml, MSE = MSE))
}








