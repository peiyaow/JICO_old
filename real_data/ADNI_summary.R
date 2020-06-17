result = as.matrix(result)
L = 50
n1 = 35
n2 = 34
n3 = 29
n = n1 + n2 + n3

G = 3
N = nrow(result)
result.list = list()
MSE.list = list()
for (i in 1:N){
  result.list[[i]] = matrix(result[i,-1], ncol = G)
  MSE.list[[i]] = matrix(result[i,-1], ncol = G)%*%c(n1, n2, n3)/n
}

# overall
MSE = apply(do.call(cbind, MSE.list), 1, mean)
se = apply(do.call(cbind, MSE.list), 1, sd)/sqrt(N)
A = round(cbind(MSE, se), digit = 3)[c(1,2,6:8,3:5),]

# per group
MSE = matrix(apply(result, 2, mean)[-1], ncol = G)
se = matrix(apply(result, 2, sd)[-1], ncol = G)/sqrt(N)
B = round(cbind(MSE, se)[,c(1,4,2,5,3,6)], digits = 3)[c(1,2,6:8,3:5),]

C = apply(cbind(B, A), 2, function(x) sprintf(x, fmt = "%0.3f"))
parenth = c(rep(c("\ (", ")\ &\ "), 3), c("\ (", ")\ \\\\"))
D = cbind(C, matrix(rep(parenth, 8), nrow = 8, byrow = T))
printt = D[, do.call(c, lapply(1:8, function(i) c(i,(i+8))))]
for (i in 1: nrow(printt)){
  cat(c("& ", printt[i,]), sep = "", fill = T)
}

apply(printt, 1, function(x) cat(x, sep = ""))

result.mean = matrix(apply(result, 2, mean)[-1], ncol = G)[c(L, (2*L):(2*L+6)), ]
result.sd = matrix(apply(result, 2, sd)[-1], ncol = G)[c(L, (2*L):(2*L+6)), ]
result.mean%*%c(n1, n2, n3)/n

row.names(result.mean) = c("2step", "iter", "global.pls", "global.pcr", "global.ridge", "group.pls", "group.pcr", "group.ridge")

g = 1
plot(gam.list, matrix(apply(result, 2, mean)[-1], ncol = G)[1:(L-1),g], xlab = "gam", ylab = "MSE")
abline(h = matrix(apply(result, 2, mean)[-1], ncol = G)[L,g])

g = 2
plot(gam.list, matrix(apply(result, 2, mean)[-1], ncol = G)[1:(L-1),g], xlab = "gam", ylab = "MSE")
abline(h = matrix(apply(result, 2, mean)[-1], ncol = G)[L,g])

g = 3
plot(gam.list, matrix(apply(result, 2, mean)[-1], ncol = G)[1:(L-1),g], xlab = "gam", ylab = "MSE")
abline(h = matrix(apply(result, 2, mean)[-1], ncol = G)[L,g])

g = 1
plot(gam.list, matrix(apply(result, 2, mean)[-1], ncol = G)[L + (1:(L-1)),g], xlab = "gam", ylab = "MSE")
abline(h = matrix(apply(result, 2, mean)[-1], ncol = G)[L*2, g])

g = 2
plot(gam.list, matrix(apply(result, 2, mean)[-1], ncol = G)[L + (1:(L-1)),g], xlab = "gam", ylab = "MSE")
abline(h = matrix(apply(result, 2, mean)[-1], ncol = G)[L*2, g])

g = 3
plot(gam.list, matrix(apply(result, 2, mean)[-1], ncol = G)[L + (1:(L-1)),g], xlab = "gam", ylab = "MSE")
abline(h = matrix(apply(result, 2, mean)[-1], ncol = G)[L*2, g])

save.image(file = "realdata_1se.RData")






