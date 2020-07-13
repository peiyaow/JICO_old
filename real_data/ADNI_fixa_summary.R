n1 = 35
n2 = 34
n3 = 29
n = n1 + n2 + n3
name = sprintf(seq(0,1, by = 0.1), fmt = "%0.1f")
A = list()
B = list()
for (nn in name){
  path = paste0("~/Documents/GitHub/continuum/result/real_data/v7/ADNI_fixa/", nn)
  setwd(path)
  result <- read.csv("result.csv", header=FALSE)
  result = as.matrix(result)
  G = 3
  N = nrow(result)
  print(N)
  result.list = list()
  MSE.list = list()
  for (i in 1:N){
    result.list[[i]] = matrix(result[i,-1], ncol = G)
    MSE.list[[i]] = matrix(result[i,-1], ncol = G)%*%c(n1, n2, n3)/n
  }
  
  # overall
  MSE = apply(do.call(cbind, MSE.list), 1, mean)
  se = apply(do.call(cbind, MSE.list), 1, sd)/sqrt(N)
  A = list.append(A, round(cbind(MSE, se), digit = 3))
  
  # group
  MSE = matrix(apply(result, 2, mean)[-1], ncol = G)
  se = matrix(apply(result, 2, sd)[-1], ncol = G)/sqrt(N)
  ix.MSEse = do.call(c, lapply(1:G, function(g) c(1, 1+G) + g-1))
  B = list.append(B, round(cbind(MSE, se)[,ix.MSEse], digits = 3))
}
A = do.call(rbind, A)
B = do.call(rbind, B)

C = apply(cbind(B, A), 2, function(x) sprintf(x, fmt = "%0.3f"))
parenth = c(rep(c("\ (", ")\ &\ "), 3), c("\ (", ")\ \\\\"))
D = cbind(C, matrix(rep(parenth, nrow(C)), nrow = nrow(C), byrow = T))
printt = D[, do.call(c, lapply(1:8, function(i) c(i,(i+8))))]
for (i in 1: nrow(printt)){
  cat(c("& ", printt[i,]), sep = "", fill = T)
}


parenth = c(rep(" & ", 3), " \\\\")
D = cbind(RANK, matrix(rep(parenth, nrow(RANK)), nrow = nrow(RANK), byrow = T))
printt = D[, do.call(c, lapply(1:4, function(i) c(i,(i+4))))]
for (i in 1: nrow(printt)){
  cat(c(printt[i,]), sep = "", fill = T)
}


RANK = list()
for (nn in name){
  path = paste0("~/Documents/GitHub/continuum/result/real_data/v7/ADNI_fixa/", nn)
  setwd(path)
  rank <- unique(read.csv("rank_iter.csv", header=FALSE)[,-(1)])
  RANK = list.append(RANK, rank)
}

RANK = unique(do.call(rbind, RANK))
row.names(RANK) = NULL  
colnames(RANK) = NULL

save(RANK, file = "rank.RData")
# a 0.2 0.75 1 
# gamma 0.25, 3, 1e10
  
  
