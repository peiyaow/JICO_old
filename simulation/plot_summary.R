n = 10
x = c(0, 0.5, 1)

AA_mean = sapply(mse.global, function(x) apply(x, 2, mean))
AA_sd = sapply(mse.global, function(x) apply(x, 2, sd)/sqrt(n))

BB_mean = sapply(mse.group, function(x) apply(x, 2, mean))
BB_sd = sapply(mse.group, function(x) apply(x, 2, sd)/sqrt(n))

CC_mean = sapply(mse.int, function(x) apply(x, 2, mean))
CC_sd = sapply(mse.int, function(x) apply(x, 2, sd)/sqrt(n))

AA_mean = t(rbind(AA_mean, BB_mean, CC_mean))
AA_sd = t(rbind(AA_sd, BB_sd, CC_mean))

# global pcr
i = 1
plot(x, AA_mean[,i], type = "l", 
     ylim = c(min(AA_mean)-0.3, max(AA_mean)+0.3), 
     lwd = 2, lty = 5, col = "green", xlab = "h", ylab = "error")
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)

# global ridge
i = 2
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "red", AA_mean[,i], lwd = 2, lty = 5)

# global pls
i = 3
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "blue", AA_mean[,i], lwd = 2, lty = 5)

# group pcr
i = 4
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)
lines(x, AA_mean[,i], lwd = 2, lty = 4, col = "green")

# group ridge
i = 5
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(x, AA_mean[,i], lwd = 2, lty = 4, col = "red")

# global pls
i = 6
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "blue", AA_mean[,i], lwd = 2, lty = 4)

# ridge.maxmin
i = 7
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "red", AA_mean[,i], lwd = 2, lty = 3)

#mse.ridge.int
i = 8
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "red", AA_mean[,i], lwd = 2, lty = 1)

#mse.pls.maxmin
i = 9
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "blue", AA_mean[,i], lwd = 2, lty = 3)

#mse.pls.int
i = 10
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "blue", AA_mean[,i], lwd = 2, lty = 1)

legend("topleft", inset=.1, legend=c("Global", "Group", "Maxmin", "Int"), lwd = 2, lty=c(5, 4, 3, 1), cex=1.2)
legend("top", inset=.1, legend=c("PCR", "Ridge", "PLS"), lwd = 2, col=c("green", "red", "blue"), cex=1.2)



