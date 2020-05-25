library(pls)
mycvplsr = function(X.list, Y.list, scale.X = TRUE){
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
  ml.pls = plsr(Y~X, validation = "CV")
  ncomp.pls = which.min(RMSEP(ml.pls)$val[1,,])-1
  if (ncomp.pls){
    Yhat.pls.list = lapply(1:G, function(g) predict(ml.pls, newdata = X.list[[g]], ncomp = ncomp.pls)[,,1])
  }else{
    Yhat.pls.list = lapply(1:G, function(g) rep(ml.pls$Ymeans, n))
  }
  Y.heter.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.pls.list[[g]])
  
  ml.pls.list = lapply(1:G, function(g) plsr(Y.heter.list[[g]] ~ X.list[[g]], validation = "CV"))
  ncomp.pls.list = lapply(1:G, function(g) which.min(RMSEP(ml.pls.list[[g]])$val[1,,])-1)
  
  if (ncomp.pls){
    beta.C0 = matrix(ml.pls$coefficients[,,ncomp.pls], nrow = p)
  }else{
    beta.C0 = matrix(0, nrow = p)
  }
  
  beta.Cind0 = lapply(1:G, function(g) if(ncomp.pls.list[[g]]) {
    matrix(ml.pls.list[[g]]$coefficients[,,ncomp.pls.list[[g]]], nrow = p)
    }else{
      matrix(0, nrow = p)
    })
  beta.C = list()
  beta.Cind = list()
  intercept = list()
  for (g in 1:G){
    beta.C[[g]] = beta.C0/scaleValues.X[[g]]*scaleValues.Y[[g]]
    beta.Cind[[g]] = beta.Cind0[[g]]/scaleValues.X[[g]]*scaleValues.Y[[g]]
    intercept[[g]] = centerValues.Y[[g]] - t(beta.C[[g]])%*%centerValues.X[[g]] - t(beta.Cind[[g]])%*%centerValues.X[[g]]
  }
  return(list(intercept = intercept, beta.C = beta.C, beta.Cind = beta.Cind, rankJ = ncomp.pls, rankA = do.call(c, ncomp.pls.list)))
}

myfixplsr = function(X.list, Y.list, rankJ, rankA, scale.X = TRUE){
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
  ml.pls = plsr(Y~X, validation = "CV")
  ncomp.pls = rankJ
  if (ncomp.pls){
    Yhat.pls.list = lapply(1:G, function(g) predict(ml.pls, newdata = X.list[[g]], ncomp = ncomp.pls)[,,1])
  }else{
    Yhat.pls.list = lapply(1:G, function(g) rep(ml.pls$Ymeans, n))
  }
  Y.heter.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.pls.list[[g]])
  
  ml.pls.list = lapply(1:G, function(g) plsr(Y.heter.list[[g]] ~ X.list[[g]], validation = "CV"))
  ncomp.pls.list = lapply(1:G, function(g) rankA[g])
  
  if (ncomp.pls){
    beta.C0 = matrix(ml.pls$coefficients[,,ncomp.pls], nrow = p)
  }else{
    beta.C0 = matrix(0, nrow = p)
  }
  
  beta.Cind0 = lapply(1:G, function(g) if(ncomp.pls.list[[g]]) {
    matrix(ml.pls.list[[g]]$coefficients[,,ncomp.pls.list[[g]]], nrow = p)
  }else{
    matrix(0, nrow = p)
  })
  beta.C = list()
  beta.Cind = list()
  intercept = list()
  for (g in 1:G){
    beta.C[[g]] = beta.C0/scaleValues.X[[g]]*scaleValues.Y[[g]]
    beta.Cind[[g]] = beta.Cind0[[g]]/scaleValues.X[[g]]*scaleValues.Y[[g]]
    intercept[[g]] = centerValues.Y[[g]] - t(beta.C[[g]])%*%centerValues.X[[g]] - t(beta.Cind[[g]])%*%centerValues.X[[g]]
  }
  return(list(intercept = intercept, beta.C = beta.C, beta.Cind = beta.Cind, rankJ = ncomp.pls, rankA = do.call(c, ncomp.pls.list)))
}

mycvpcr = function(X.list, Y.list, scale.X = TRUE){
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
  
  ml.pcr = pcr(Y~X, validation = "CV", scale = F)
  ncomp.pcr = which.min(RMSEP(ml.pcr)$val[1,,])-1
  
  if (ncomp.pcr){
    Yhat.pcr.list = lapply(1:G, function(g) predict(ml.pcr, newdata = X.list[[g]], ncomp = ncomp.pcr)[,,1])
  }else{
    Yhat.pcr.list = lapply(1:G, function(g) rep(ml.pcr$Ymeans, n))
  }
  Y.heter.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.pcr.list[[g]])
  
  ml.pcr.list = lapply(1:G, function(g) pcr(Y.heter.list[[g]] ~ X.list[[g]], validation = "CV"))
  ncomp.pcr.list = lapply(1:G, function(g) which.min(RMSEP(ml.pcr.list[[g]])$val[1,,])-1)
  
  if (ncomp.pcr){
    beta.C0 = matrix(ml.pcr$coefficients[,,ncomp.pcr], nrow = p)
  }else{
    beta.C0 = matrix(0, nrow = p)
  }
  
  beta.Cind0 = lapply(1:G, function(g) if(ncomp.pcr.list[[g]]) {
    matrix(ml.pcr.list[[g]]$coefficients[,,ncomp.pcr.list[[g]]], nrow = p)
  }else{
    matrix(0, nrow = p)
  })
  beta.C = list()
  beta.Cind = list()
  intercept = list()
  for (g in 1:G){
    beta.C[[g]] = beta.C0/scaleValues.X[[g]]*scaleValues.Y[[g]]
    beta.Cind[[g]] = beta.Cind0[[g]]/scaleValues.X[[g]]*scaleValues.Y[[g]]
    intercept[[g]] = centerValues.Y[[g]] - t(beta.C[[g]])%*%centerValues.X[[g]] - t(beta.Cind[[g]])%*%centerValues.X[[g]]
  }
  return(list(intercept = intercept, beta.C = beta.C, beta.Cind = beta.Cind, rankJ = ncomp.pcr, rankA = do.call(c, ncomp.pcr.list)))
  
}

  