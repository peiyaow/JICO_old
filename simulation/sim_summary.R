G = 2
result = as.matrix(result)
n = nrow(result)
result.list = list()
for (i in 1:n){
  result.list[[i]] = matrix(as.vector(result[i,-1]), ncol = G)
}

result.mean = matrix(apply(result, 2, mean)[-1], ncol = G)[c(L, (2*L):(2*L+6)), ]
result.sd = matrix(apply(result, 2, sd)[-1], ncol = G)[c(L, (2*L):(2*L+6)), ]
result.mean%*%c(n1, n2)/n

row.names(result.mean) = c("2step", "iter", "global.pls", "global.pcr", "global.ridge", "group.pls", "group.pcr", "group.ridge")

g = 1
plot(gam.list, matrix(apply(result, 2, mean)[-1], ncol = G)[1:(L-1),g], xlab = "gam", ylab = "MSE")
abline(h = matrix(apply(result, 2, mean)[-1], ncol = G)[L,g])

g = 2
plot(gam.list, matrix(apply(result, 2, mean)[-1], ncol = G)[1:(L-1),g], xlab = "gam", ylab = "MSE")
abline(h = matrix(apply(result, 2, mean)[-1], ncol = G)[L,g])

g = 1
plot(gam.list, matrix(apply(result, 2, mean)[-1], ncol = G)[L + (1:(L-1)),g], xlab = "gam", ylab = "MSE")
abline(h = matrix(apply(result, 2, mean)[-1], ncol = G)[L*2, g])

g = 2
plot(gam.list, matrix(apply(result, 2, mean)[-1], ncol = G)[L + (1:(L-1)),g], xlab = "gam", ylab = "MSE")
abline(h = matrix(apply(result, 2, mean)[-1], ncol = G)[L*2, g])

n = nrow(rank)
rank.list = list()
for (i in 1:n){
  rank.list[[i]] = matrix(rank[i,-1], ncol = G+1)
}
