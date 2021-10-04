PCR_mean = AA_mean[2*(L+1)+1:(L+1)]
PCR_sd = AA_sd[2*(L+1)+1:(L+1)]

PLS_mean = AA_mean[2*(L+1)+1:(L+1)]
PLS_sd = AA_sd[2*(L+1)+1:(L+1)]

OLS1_mean = AA_mean[1:(L+1)]
OLS1_sd = AA_sd[1:(L+1)]

OLS2_mean = AA_mean[3*(L+1)+1:(L+1)]
OLS2_sd = AA_sd[3*(L+1)+1:(L+1)]

MEAN = rbind(PCR_mean, PLS_mean, OLS1_mean, OLS2_mean)
SD = rbind(PCR_sd, PLS_sd, OLS1_sd, OLS2_sd)

y_min = min(MEAN)
y_max = max(MEAN)
setwd("~/Documents/GitHub/continuum/result/sim/v5/")
png("3curve.png", width=1000, height=800)
par(cex = 1.7, cex.lab = 1.2, mar = c(5, 4, 4, 2) + 0.1, lwd = 2)
plot(a, MEAN[1, ], type = "l", 
     ylim = c(y_min, y_max), 
     lwd = 2, lty = 2, col = 2, # cex.lab = 1.5,
     xlab = expression(paste("a = ", gamma, "/(1+", gamma, ")")), ylab = "MSE")
polygon(c(a, rev(a)), c(MEAN[1,]+SD[1,], rev(MEAN[1,]-SD[1,])), 
        col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)

polygon(c(a, rev(a)), c(MEAN[2,]+SD[2,], rev(MEAN[2,]-SD[2,])),
        col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)
lines(a, MEAN[2,],
      col = "green", lwd = 2, lty = 3)

polygon(c(a, rev(a)), c(MEAN[4,]+SD[4,], rev(MEAN[4,]-SD[4,])),
        col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(a, MEAN[4,],
      col = "blue", lwd = 2, lty = 4)

legend("topright", inset=.1, legend=c("PCR", 
                                      "PLS", 
                                      "OLS"),
       col=c("red", 
             "green", 
             "blue"), lty = c(2, 3, 5), lwd = 2)#, cex=1.3)
dev.off()
