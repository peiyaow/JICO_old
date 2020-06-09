MSE = do.call(rbind, MSE)
L = 50
a = seq(0, 1, length.out = L+1)
# gam.list = a/(1-a)
# gam.list[L+1] = 1e10

MSE = as.matrix(result.list[[1]])

MSE = matrix(apply(result, 2, mean)[-1], ncol = G)
plot(a, MSE[1:(L+1),1], type = "l", ylim = c(min(MSE), max(MSE)))
lines(a, MSE[(L+1)+1:(L+1),1])
lines(a, MSE[2*(L+1)+1:(L+1),1])
lines(a, MSE[3*(L+1)+1:(L+1),1])
lines(a, MSE[4*(L+1)+1:(L+1),1])
abline(h = MSE[-(1:(5*(L+1))),1])

plot(a, MSE[1:(L+1),2], type = "l", ylim = c(min(MSE), max(MSE)))
lines(a, MSE[(L+1)+1:(L+1),2])
lines(a, MSE[2*(L+1)+1:(L+1),2])
lines(a, MSE[3*(L+1)+1:(L+1),2])
lines(a, MSE[4*(L+1)+1:(L+1),2])
abline(h = MSE[-(1:(5*(L+1))),2])

MSE = apply(MSE, 1, mean)
plot(a, MSE[1:(L+1)], type = "l", ylim = c(min(MSE), max(MSE)))
lines(a, MSE[(L+1)+1:(L+1)])
lines(a, MSE[2*(L+1)+1:(L+1)])
lines(a, MSE[3*(L+1)+1:(L+1)])
lines(a, MSE[4*(L+1)+1:(L+1)])
abline(h = MSE[-(1:(5*(L+1)))])

