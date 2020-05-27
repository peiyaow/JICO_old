test = function(X.list, Y.list, center.X = TRUE, scale.X = TRUE, center.Y = TRUE, scale.Y = TRUE){
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
  
  ml.pls = plsr(Y~X, validation = "CV", scale = F, center = F)
  ncomp.pls = selectNcomp(ml.pls, method = "randomization", plot = F)
  print(which.min(RMSEP(ml.pls)$val[1,,])-1)
  print(ncomp.pls)
}

test(X.list, Y.list)

ml.step1 = cv.continuum.step1(X.list, Y.list, lambda = 0, gam = 1, nfolds = 10, m = 10, scale.X = T, center.X = T, plot = T)
rankJ.1se = ml.step1$rankJ.1se
rankJ.min = ml.step1$rankJ.min
rankJ.1se
rankJ.min

# ml.step2 = lapply(1:G, function(g) cv.continuum.step2(ml.step1$ml.1se$X.heter[[g]], ml.step1$ml.1se$Y.heter[[g]], ml.step1$ml.1se$C, 
#                                                       ml.step1$ml.1se$XOrigin[[g]], ml.step1$ml.1se$YOriginRes[[g]], 
#                                                       lambda = 0, gam = 1, nfolds = 10, m = 5, 
#                                                       ml.step1$ml.1se$centerValues.X[[g]], ml.step1$ml.1se$scaleValues.X[[g]], 
#                                                       ml.step1$ml.1se$scaleValues.Y[[g]], 
#                                                       plot = T))

ml.step2 = lapply(1:G, function(g) cv.continuum.step2(ml.step1$ml.min$X.heter[[g]], ml.step1$ml.min$Y.heter[[g]], ml.step1$ml.min$C, 
                                                      ml.step1$ml.min$XOrigin[[g]], ml.step1$ml.min$YOriginRes[[g]], 
                                                      lambda = 0, gam = 1, nfolds = 10, m = 5, 
                                                      ml.step1$ml.min$centerValues.X[[g]], ml.step1$ml.min$scaleValues.X[[g]], 
                                                      ml.step1$ml.min$scaleValues.Y[[g]], 
                                                      plot = T))

rankA.1se = sapply(1:G, function(g) ml.step2[[g]]$rankA.1se)
rankA.min = sapply(1:G, function(g) ml.step2[[g]]$rankA.min)

rankA.1se
rankA.min








MSE = do.call(rbind, haha$MSE)
rMSE = apply(MSE, 2, function(x) mean(sqrt(x)))
plot(rMSE)
absBest = min(rMSE)

MSEsd = apply(MSE, 2, function(x) sd(sqrt(x)))/sqrt(nfolds)
min(which((rMSE - MSEsd) < absBest)) - 1

ml.homo = continuum.ridge.fix(X, Y, G, lambda = lambda, gam = 1, om = ncomp.pls)
C = ml.homo$C
beta.C = C2beta(X, Y, C, lambda = lambda)$beta
Yhat.homo.list = lapply(1:G, function(g) X.list[[g]]%*%beta.C)
Y.heter.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.homo.list[[g]])
X.homo.list = lapply(1:G, function(g) X.list[[g]]%*%C%*%SOLVE(t(C)%*%C)%*%t(C))
X.heter.list = lapply(1:G, function(g) X.list[[g]] - X.homo.list[[g]])

temp = X.heter.list
for (g in 1:G){
  temp[[g]] <- temp[[g]]%*%(diag(p) - C%*%SOLVE(t(C)%*%C)%*%t(C))
}
  
ml.pls.list = lapply(1:G, function(g) plsr(Y.heter.list[[g]] ~ temp[[g]], validation = "CV", center = F, scale = F))
ncomp.pls.list = lapply(1:G, function(g) selectNcomp(ml.pls.list[[g]], method = "randomization", plot = T))
ncomp.pls.list
sapply(1:G, function(g) which.min(RMSEP(ml.pls.list[[g]])$val[1,,])-1)



