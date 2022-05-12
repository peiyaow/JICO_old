library(pls)
library(caret)
library(glmnet)
library(methods)

source("~/Documents/GitHub/continuum/function/JICO.R")


##==============Define the parameters for this toy example==============

G = 2                                         # Number of groups
n1 = 50                                       # sample size of the first group
n2 = 50                                       # sample size of the second group
n = n1 + n2                                   # sample size of all the data
p = 200                                       # dimension of the data
r = 1                                         # rank of the joint structure
r1 = 1                                        # rank of the 1st individual structure
r2 = 1                                        # rank of the 2nd individual structure
gam = 1e12                                    # infinity value for PCR model, set gam=0 for OLS model, gam=0.5 for PLS model, gam >= 1e10 for PCR model
cv.maxrank = 1                                # the maximum rank for joint and individual matrices when tuning
cv.gamma = 1e12                               # the gamma used when tuning
alpha = rep(1, r)                             # underlying true coefficient of the joint structure
alpha1 = rep(1, r1)                           # underlying true coefficient of the 1st individual structure
alpha2 = rep(1, r2)                           # underlying true coefficient of the 2nd individual structure

##==============Generate the data==============
X1 = mvrnorm(n1, rep(0, p), diag(p))          # Generate the observed covariate matrix of the first group
X2 = mvrnorm(n2, rep(0, p), diag(p))          # Generate the observed covariate matrix of the second group

X = rbind(X1, X2)
X.list = list(X1, X2)

V = matrix(svd(X)$v[,1:r], ncol = r)%*%rep(1/sqrt(r), r)                              # Use svd to generate the underlying true w in the joint model 
V1 = matrix(svd(X1%*%(diag(p) - V%*%t(V)))$v[,1:r1], ncol = r1)%*%rep(1/sqrt(r1), r1) # Use svd to generate the underlying true w in the 1st individual model 
V2 = matrix(svd(X2%*%(diag(p) - V%*%t(V)))$v[,1:r2], ncol = r2)%*%rep(1/sqrt(r2), r2) # Use svd to generate the underlying true w in the 2nd individual model 

e1 = rnorm(n1)*.2                             # Generate the error vector in the 1st group
Y1 = X1%*%V%*%alpha + X1%*%V1%*%alpha1 + e1   # Generate the response vector in the 1st group
e2 = rnorm(n2)*.2                             # Generate the error vector in the 2nd group
Y2 = X2%*%V%*%alpha + X2%*%V2%*%alpha2 + e2   # Generate the response vector in the 2nd group

Y = rbind(Y1, Y2)   
Y.list = list(Y1, Y2)

##==============Rank selection through cross validation==============

## Enumerate the combinations of hyperparameters to be tuned in a list: 
cv.parameter.set = parameter.set.G_2(maxrankA = cv.maxrank, maxrankJ = cv.maxrank, gamma = cv.gamma)

## fit the model and use CV to find the best parameters
cv.ml.JICO = cv.continnum.iter(X.list, Y.list, parameter.set = cv.parameter.set, 
                               criteria = "min", nfold = 5, maxiter = 300, 
                               center.X = F, scale.X = F, center.Y = F, scale.Y = F, 
                               orthIndiv = F)
print(paste0("The best parameters:"))
print(cv.ml.JICO$parameter)

##==============Model fitting example==============
# fitting oracle
ml.JICO = continuum.multigroup.iter(X.list, Y.list, maxiter = 300, gam = gam, rankJ = r, rankA = c(r1, r2),
                          center.X = F, scale.X = F, center.Y = F, scale.Y = F, 
                          orthIndiv = F)

# predicted response value
Yhat.list = lapply(1:G, function(g) as.numeric(ml.JICO$intercept[[g]])+ X.list[[g]]%*%ml.JICO$beta.C[[g]] + X.list[[g]]%*%ml.JICO$beta.Cind[[g]])






