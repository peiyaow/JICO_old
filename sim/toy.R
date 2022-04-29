library(pls)
library(caret)
library(glmnet)
library(methods)

source("~/Documents/GitHub/continuum-2/function/JICO.R")


##==============Define the parameters for this toy example==============

G = 2                                         # Number of groups
n1 = 50                                       # sample size of the first group
n2 = 50                                       # sample size of the second group
n = n1 + n2                                   # sample size of all the data
p = 200                                       # dimension of the data
r = 1                                         # rank of the joint structure
r1 = 1                                        # rank of the first individual structure
r2 = 1                                        # rank of the second individual structure
gam = 1e12                                    # infinity value for PCR model, set gam=0 for OLS model, gam=0.5 for PLS model, gam > 1e10 for PCR model
cv.maxrank = 2                                # the maximum rank for joint and individual matrices when tuning
cv.gamma = 1e12                               # the gamma used when tuning
alpha = rep(1, r)                             # underlying true coefficient of the joint structure
alpha1 = rep(1, r1)                           # underlying true coefficient of the 1st individual structure
alpha2 = rep(1, r2)                           # underlying true coefficient of the 2nd individual structure

##==============Generate the data==============
X1 = mvrnorm(n1, rep(0, p), diag(p))          # Generate the observed covariate matrix of the first group
X2 = mvrnorm(n2, rep(0, p), diag(p))          # Generate the observed covariate matrix of the second group
X = rbind(X1, X2)

q = r                                         # Define the underlying true K in our joint model
q1 = r1                                       # Define the underlying true K in our 1st individual model
q2 = r2                                       # Define the underlying true K in our 2nd individual model
V = matrix(svd(X)$v[,1:q], ncol = q)%*%rep(1/sqrt(q), q)                              #Use svd to generate the underlying true w in the joint model 
V1 = matrix(svd(X1%*%(diag(p) - V%*%t(V)))$v[,1:q1], ncol = q1)%*%rep(1/sqrt(q1), q1) #Use svd to generate the underlying true w in the 1st individual model 
V2 = matrix(svd(X2%*%(diag(p) - V%*%t(V)))$v[,1:q2], ncol = q2)%*%rep(1/sqrt(q2), q2) #Use svd to generate the underlying true w in the 2nd individual model 
e1 = rnorm(n1)*.2                             #Generate the error vector in the first group
Y1 = X1%*%V%*%alpha + X1%*%V1%*%alpha1 + e1   #Generate the response vector in the first group
e2 = rnorm(n2)*.2                             #Generate the error vector in the second group
Y2 = X2%*%V%*%alpha + X2%*%V2%*%alpha2 + e2   #Generate the response vector in the second group

Y = rbind(Y1, Y2)   
X.list = list(X1, X2)
Y.list = list(Y1, Y2)
r.list = list(r1, r2) 

##==============Rank selection through cross validation==============

## Define the tuning parameters list
cv.parameter.set <- list() 
for(rankA1 in 1:cv.maxrank)
  for(rankA2 in 1:cv.maxrank)
    for (rankJ in 1:cv.maxrank){
      cv.parameter.set = list.append(cv.parameter.set,list(rankA = c(rankA1,rankA2),rankJ = rankJ,gam = cv.gamma))
    }


## cross validation for finding the best parameters
cv.ml.JICO = cv.continnum.iter(X.list, Y.list, lambda = 0, parameter.set = cv.parameter.set,criteria =  "1se", nfold = 10, maxiter = 100,
                               center.X = F, scale.X = F, center.Y = F, scale.Y = F, orthIndiv = F)
print(paste0("The best parameters:"))
print(cv.ml.JICO$parameter)



##==============Model fitting example==============
# fitting oracle
ml.JICO = continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 300, gam = gam, rankJ = 1, rankA = c(1, 1),
                          center.X = F, scale.X = F, center.Y = F, scale.Y = F, orthIndiv = F)

#predicted response value
Yhat.list = lapply(1:G, function(g) as.numeric(ml.JICO$intercept[[g]])+ X.list[[g]]%*%ml.JICO$beta.C[[g]] + X.list[[g]]%*%ml.JICO$beta.Cind[[g]])






