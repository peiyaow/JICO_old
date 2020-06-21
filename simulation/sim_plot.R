# MSE = do.call(rbind, MSE)
# L = 50
# a = seq(0, 1, length.out = L+1)
# # gam.list = a/(1-a)
# # gam.list[L+1] = 1e10
# 
# MSE = as.matrix(result.list[[1]])
# 
# MSE = matrix(apply(result, 2, mean)[-1], ncol = G)
# plot(a, MSE[1:(L+1),1], type = "l", ylim = c(min(MSE), max(MSE)))
# lines(a, MSE[(L+1)+1:(L+1),1])
# lines(a, MSE[2*(L+1)+1:(L+1),1])
# lines(a, MSE[3*(L+1)+1:(L+1),1])
# lines(a, MSE[4*(L+1)+1:(L+1),1])
# abline(h = MSE[-(1:(5*(L+1))),1])
# 
# plot(a, MSE[1:(L+1),2], type = "l", ylim = c(min(MSE), max(MSE)))
# lines(a, MSE[(L+1)+1:(L+1),2])
# lines(a, MSE[2*(L+1)+1:(L+1),2])
# lines(a, MSE[3*(L+1)+1:(L+1),2])
# lines(a, MSE[4*(L+1)+1:(L+1),2])
# abline(h = MSE[-(1:(5*(L+1))),2])
# 
# MSE = apply(MSE, 1, mean)
# plot(a, MSE[1:(L+1)], type = "l", ylim = c(min(MSE), max(MSE)))
# lines(a, MSE[(L+1)+1:(L+1)])
# lines(a, MSE[2*(L+1)+1:(L+1)])
# lines(a, MSE[3*(L+1)+1:(L+1)])
# lines(a, MSE[4*(L+1)+1:(L+1)])
# abline(h = MSE[-(1:(5*(L+1)))])

###########################
result = as.matrix(result)
n = nrow(result)
L = 50
a = seq(0, 1, length.out = L+1)

MSE = matrix(apply(result, 2, mean)[-1], ncol = G)
se = matrix(apply(result, 2, sd)[-1], ncol = G)/sqrt(n)
y_max = max(MSE)
y_min = min(MSE)

################################ first group
i = 1
AA_mean = MSE[,i]
AA_sd = se[,i]

#100
plot(a, AA_mean[1:(L+1)], type = "l", 
     ylim = c(y_min, y_max), 
     lwd = 1, lty = 1, col = "green", xlab = "a", ylab = "MSE")
polygon(c(a, rev(a)), c(AA_mean[1:(L+1)]+AA_sd[1:(L+1)], rev(AA_mean[1:(L+1)]-AA_sd[1:(L+1)])), 
        col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)

#200
polygon(c(a, rev(a)), c(AA_mean[(L+1)+1:(L+1)]+AA_sd[(L+1)+1:(L+1)], rev(AA_mean[(L+1)+1:(L+1)]-AA_sd[(L+1)+1:(L+1)])), 
        col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(a, AA_mean[(L+1)+1:(L+1)],
      col = "red", lwd = 1, lty = 2)

#111
polygon(c(a, rev(a)), c(AA_mean[2*(L+1)+1:(L+1)]+AA_sd[2*(L+1)+1:(L+1)], rev(AA_mean[2*(L+1)+1:(L+1)]-AA_sd[2*(L+1)+1:(L+1)])), 
        col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(a, AA_mean[2*(L+1)+1:(L+1)],
      col = "blue", lwd = 1, lty = 3)

#011
polygon(c(a, rev(a)), c(AA_mean[3*(L+1)+1:(L+1)]+AA_sd[3*(L+1)+1:(L+1)], rev(AA_mean[3*(L+1)+1:(L+1)]-AA_sd[3*(L+1)+1:(L+1)])), 
        col = rgb(0, 255, 255, 20, maxColorValue=255), border = FALSE)
lines(a, AA_mean[3*(L+1)+1:(L+1)],
      col = rgb(0, 255, 255, 255, maxColorValue=255), lwd = 1, lty = 4)

#022
polygon(c(a, rev(a)), c(AA_mean[4*(L+1)+1:(L+1)]+AA_sd[4*(L+1)+1:(L+1)], rev(AA_mean[4*(L+1)+1:(L+1)]-AA_sd[4*(L+1)+1:(L+1)])), 
        col = rgb(255, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(a, AA_mean[4*(L+1)+1:(L+1)],
      col = rgb(255, 0, 255, 255, maxColorValue=255), lwd = 1, lty = 5)

legend("topright", inset=.1, legend=c("r = 1, r1 = 0, r2 = 0", 
                                      "r = 2, r1 = 0, r2 = 0", 
                                      "r = 1, r1 = 1, r2 = 1",
                                      "r = 0, r1 = 1, r2 = 1",
                                      "r = 0, r1 = 2, r2 = 2"),
       col=c("green", 
             "red", 
             "blue",
             rgb(0, 255, 255, 255, maxColorValue=255),
             rgb(255, 0, 255, 255, maxColorValue=255)), lty = seq(1,5), cex=1)

############################################### second group
i = 2
AA_mean = MSE[,i]
AA_sd = se[,i]

#100
plot(a, AA_mean[1:(L+1)], type = "l", 
     ylim = c(y_min, y_max), 
     lwd = 1, lty = 1, col = "green", xlab = "a", ylab = "MSE")
polygon(c(a, rev(a)), c(AA_mean[1:(L+1)]+AA_sd[1:(L+1)], rev(AA_mean[1:(L+1)]-AA_sd[1:(L+1)])), 
        col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)

#200
polygon(c(a, rev(a)), c(AA_mean[(L+1)+1:(L+1)]+AA_sd[(L+1)+1:(L+1)], rev(AA_mean[(L+1)+1:(L+1)]-AA_sd[(L+1)+1:(L+1)])), 
        col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(a, AA_mean[(L+1)+1:(L+1)],
      col = "red", lwd = 1, lty = 2)

#111
polygon(c(a, rev(a)), c(AA_mean[2*(L+1)+1:(L+1)]+AA_sd[2*(L+1)+1:(L+1)], rev(AA_mean[2*(L+1)+1:(L+1)]-AA_sd[2*(L+1)+1:(L+1)])), 
        col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(a, AA_mean[2*(L+1)+1:(L+1)],
      col = "blue", lwd = 1, lty = 3)

#011
polygon(c(a, rev(a)), c(AA_mean[3*(L+1)+1:(L+1)]+AA_sd[3*(L+1)+1:(L+1)], rev(AA_mean[3*(L+1)+1:(L+1)]-AA_sd[3*(L+1)+1:(L+1)])), 
        col = rgb(0, 255, 255, 20, maxColorValue=255), border = FALSE)
lines(a, AA_mean[3*(L+1)+1:(L+1)],
      col = rgb(0, 255, 255, 255, maxColorValue=255), lwd = 1, lty = 4)

#022
polygon(c(a, rev(a)), c(AA_mean[4*(L+1)+1:(L+1)]+AA_sd[4*(L+1)+1:(L+1)], rev(AA_mean[4*(L+1)+1:(L+1)]-AA_sd[4*(L+1)+1:(L+1)])), 
        col = rgb(255, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(a, AA_mean[4*(L+1)+1:(L+1)],
      col = rgb(255, 0, 255, 255, maxColorValue=255), lwd = 1, lty = 5)

legend("topright", inset=.1, legend=c("r = 1, r1 = 0, r2 = 0", 
                                      "r = 2, r1 = 0, r2 = 0", 
                                      "r = 1, r1 = 1, r2 = 1",
                                      "r = 0, r1 = 1, r2 = 1",
                                      "r = 0, r1 = 2, r2 = 2"),
       col=c("green", 
             "red", 
             "blue",
             rgb(0, 255, 255, 255, maxColorValue=255),
             rgb(255, 0, 255, 255, maxColorValue=255)), lty = seq(1,5), cex=1)


############## overall ###############
setwd("~/Documents/GitHub/continuum/result/sim/v4/")
result = as.matrix(result)
G = 2
n = nrow(result)
L = 50
a = seq(0, 1, length.out = L+1)
MSE.list = list()
for (i in 1:n){
  MSE.list[[i]] = apply(matrix(as.vector(result[i,-1]), ncol = G), 1, mean)
}
MSE = apply(do.call(rbind, MSE.list), 2, mean)
se = apply(do.call(rbind, MSE.list), 2, sd)/sqrt(n)

y_max = max(MSE[1:(5*(L+1))])
y_min = min(MSE[1:(5*(L+1))])

AA_mean = MSE
AA_sd = se

# ------------------------------- PCR/PLS -------------------------------
png("PLS.png", width=1000, height=800)

#100
plot(a, AA_mean[1:(L+1)], type = "l", 
     ylim = c(y_min, y_max), 
     lwd = 2, lty = 3, col = "gray", cex.lab = 1.3,
     xlab = expression(paste("a = ", gamma, "/(1+", gamma, ")")), ylab = "MSE")
polygon(c(a, rev(a)), c(AA_mean[1:(L+1)]+AA_sd[1:(L+1)], rev(AA_mean[1:(L+1)]-AA_sd[1:(L+1)])), 
        col = rgb(220, 220, 220, 80, maxColorValue=255), border = FALSE)

#200
polygon(c(a, rev(a)), c(AA_mean[(L+1)+1:(L+1)]+AA_sd[(L+1)+1:(L+1)], rev(AA_mean[(L+1)+1:(L+1)]-AA_sd[(L+1)+1:(L+1)])), 
        col = rgb(220, 220, 220, 80, maxColorValue=255), border = FALSE)
lines(a, AA_mean[(L+1)+1:(L+1)],
      col = "gray", lwd = 2, lty = 2)

#111
polygon(c(a, rev(a)), c(AA_mean[2*(L+1)+1:(L+1)]+AA_sd[2*(L+1)+1:(L+1)], rev(AA_mean[2*(L+1)+1:(L+1)]-AA_sd[2*(L+1)+1:(L+1)])), 
        col = rgb(255, 0, 0, 50, maxColorValue=255), border = FALSE)
lines(a, AA_mean[2*(L+1)+1:(L+1)],
      col = "red", lwd = 2, lty = 1)

#011
polygon(c(a, rev(a)), c(AA_mean[3*(L+1)+1:(L+1)]+AA_sd[3*(L+1)+1:(L+1)], rev(AA_mean[3*(L+1)+1:(L+1)]-AA_sd[3*(L+1)+1:(L+1)])), 
        col = rgb(220, 220, 220, 80, maxColorValue=255), border = FALSE)
lines(a, AA_mean[3*(L+1)+1:(L+1)],
      col = "gray", lwd = 2, lty = 4)

#022
polygon(c(a, rev(a)), c(AA_mean[4*(L+1)+1:(L+1)]+AA_sd[4*(L+1)+1:(L+1)], rev(AA_mean[4*(L+1)+1:(L+1)]-AA_sd[4*(L+1)+1:(L+1)])), 
        col = rgb(220, 220, 220, 80, maxColorValue=255), border = FALSE)
lines(a, AA_mean[4*(L+1)+1:(L+1)],
      col = "gray", lwd = 2, lty = 5)

legend("topleft", inset=.1, legend=c( "K = 1, K1 = 0, K2 = 0", 
                                      "K = 2, K1 = 0, K2 = 0", 
                                      "K = 1, K1 = 1, K2 = 1",
                                      "K = 0, K1 = 1, K2 = 1",
                                      "K = 0, K1 = 2, K2 = 2"),
       col=c("gray", 
             "gray", 
             "red",
             "gray",
             "gray"), lty = c(3,2,1,4,5), lwd = 2, cex=1.3)
dev.off()
# ------------------------------- OLS2 -------------------------------
png("OLS2.png", width=1000, height=800)

#100
plot(a, AA_mean[1:(L+1)], type = "l", 
     ylim = c(y_min, y_max), 
     lwd = 2, lty = 3, col = "gray", cex.lab = 1.3,
     xlab = expression(paste("a = ", gamma, "/(1+", gamma, ")")), ylab = "MSE")
polygon(c(a, rev(a)), c(AA_mean[1:(L+1)]+AA_sd[1:(L+1)], rev(AA_mean[1:(L+1)]-AA_sd[1:(L+1)])), 
        col = rgb(220, 220, 220, 80, maxColorValue=255), border = FALSE)

#200
polygon(c(a, rev(a)), c(AA_mean[(L+1)+1:(L+1)]+AA_sd[(L+1)+1:(L+1)], rev(AA_mean[(L+1)+1:(L+1)]-AA_sd[(L+1)+1:(L+1)])), 
        col = rgb(220, 220, 220, 80, maxColorValue=255), border = FALSE)
lines(a, AA_mean[(L+1)+1:(L+1)],
      col = "gray", lwd = 2, lty = 2)

#111
polygon(c(a, rev(a)), c(AA_mean[2*(L+1)+1:(L+1)]+AA_sd[2*(L+1)+1:(L+1)], rev(AA_mean[2*(L+1)+1:(L+1)]-AA_sd[2*(L+1)+1:(L+1)])), 
        col = rgb(220, 220, 220, 80, maxColorValue=255), border = FALSE)
lines(a, AA_mean[2*(L+1)+1:(L+1)],
      col = "gray", lwd = 2, lty = 4)

#022
polygon(c(a, rev(a)), c(AA_mean[4*(L+1)+1:(L+1)]+AA_sd[4*(L+1)+1:(L+1)], rev(AA_mean[4*(L+1)+1:(L+1)]-AA_sd[4*(L+1)+1:(L+1)])), 
        col = rgb(220, 220, 220, 80, maxColorValue=255), border = FALSE)
lines(a, AA_mean[4*(L+1)+1:(L+1)],
      col = "gray", lwd = 2, lty = 5)

#011
polygon(c(a, rev(a)), c(AA_mean[3*(L+1)+1:(L+1)]+AA_sd[3*(L+1)+1:(L+1)], rev(AA_mean[3*(L+1)+1:(L+1)]-AA_sd[3*(L+1)+1:(L+1)])), 
        col = rgb(255, 0, 0, 50, maxColorValue=255), border = FALSE)
lines(a, AA_mean[3*(L+1)+1:(L+1)],
      col = "red", lwd = 2, lty = 1)




legend("topright", inset=.1, legend=c("K = 1, K1 = 0, K2 = 0", 
                                      "K = 2, K1 = 0, K2 = 0", 
                                      "K = 1, K1 = 1, K2 = 1",
                                      "K = 0, K1 = 1, K2 = 1",
                                      "K = 0, K1 = 2, K2 = 2"),
       col=c("gray", 
             "gray", 
             "gray",
             "red",
             "gray"), lty = c(3,2,4,1,5), lwd = 2, cex=1.3)
dev.off()
# ------------------------------- OLS1 -------------------------------
png("OLS1.png", width=1000, height=800)

#200
plot(a, AA_mean[(L+1)+1:(L+1)], ylim = c(y_min, y_max), 
     type = "l", 
     col = "gray", lwd = 2, lty = 2, cex.lab = 1.3,
     xlab = expression(paste("a = ", gamma, "/(1+", gamma, ")")), ylab = "MSE")
polygon(c(a, rev(a)), c(AA_mean[(L+1)+1:(L+1)]+AA_sd[(L+1)+1:(L+1)], rev(AA_mean[(L+1)+1:(L+1)]-AA_sd[(L+1)+1:(L+1)])), 
        col = rgb(220, 220, 220, 80, maxColorValue=255), border = FALSE)

#111
polygon(c(a, rev(a)), c(AA_mean[2*(L+1)+1:(L+1)]+AA_sd[2*(L+1)+1:(L+1)], rev(AA_mean[2*(L+1)+1:(L+1)]-AA_sd[2*(L+1)+1:(L+1)])), 
        col = rgb(220, 220, 220, 80, maxColorValue=255), border = FALSE)
lines(a, AA_mean[2*(L+1)+1:(L+1)],
      lwd = 2, lty = 3, col = "gray")

#011
polygon(c(a, rev(a)), c(AA_mean[3*(L+1)+1:(L+1)]+AA_sd[3*(L+1)+1:(L+1)], rev(AA_mean[3*(L+1)+1:(L+1)]-AA_sd[3*(L+1)+1:(L+1)])), 
        col = rgb(220, 220, 220, 80, maxColorValue=255), border = FALSE)
lines(a, AA_mean[3*(L+1)+1:(L+1)],
      col = "gray", lwd = 2, lty = 4)

#022
polygon(c(a, rev(a)), c(AA_mean[4*(L+1)+1:(L+1)]+AA_sd[4*(L+1)+1:(L+1)], rev(AA_mean[4*(L+1)+1:(L+1)]-AA_sd[4*(L+1)+1:(L+1)])), 
        col = rgb(220, 220, 220, 80, maxColorValue=255), border = FALSE)
lines(a, AA_mean[4*(L+1)+1:(L+1)],
      col = "gray", lwd = 2, lty = 5)

#100
lines(a, AA_mean[1:(L+1)], ylim = c(y_min, y_max), 
      type = "l", 
      col = "red", lwd = 2, lty = 1)
polygon(c(a, rev(a)), c(AA_mean[1:(L+1)]+AA_sd[1:(L+1)], rev(AA_mean[1:(L+1)]-AA_sd[1:(L+1)])), 
        col = rgb(255, 0, 0, 50, maxColorValue=255), border = FALSE)

legend("topleft", inset=.1, legend=c("K = 1, K1 = 0, K2 = 0", 
                                     "K = 2, K1 = 0, K2 = 0", 
                                     "K = 1, K1 = 1, K2 = 1",
                                     "K = 0, K1 = 1, K2 = 1",
                                     "K = 0, K1 = 2, K2 = 2"),
       col=c("red", 
             "gray", 
             "gray",
             "gray",
             "gray"), lty = 1:5, lwd = 2, cex=1.3)
dev.off()

