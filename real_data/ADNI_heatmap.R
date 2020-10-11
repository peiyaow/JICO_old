library(pls)
library(caret)
library(glmnet)
library(methods)

source("~/Documents/GitHub/continuum/data/ADNI2/loaddata_mac.R")
source("~/Documents/GitHub/continuum/function/jive_continuum.R")
source("~/Documents/GitHub/continuum/function/cv_multigroup.R")

label.level = levels(label)
X.list = lapply(label.level, function(l) X[label == l,])
Y.list = lapply(label.level, function(l) matrix(Y[label == l]))
# rank_iter = as.matrix(rank_iter)
ml.iter = continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 200,
                                    gam = 1, 
                                    rankJ = 1, 
                                    rankA = c(1,1,1),
                                    orthIndiv = F)

J = do.call(rbind, ml.iter$J)
# dist(t(J))
hc = hclust(dist(t(J))^2,"ward.D")
hc3 = hclust(dist(ml.iter$J[[3]])^2,"ward.D")
hc2 = hclust(dist(ml.iter$J[[2]])^2,"ward.D")
hc1 = hclust(dist(ml.iter$J[[1]])^2,"ward.D")

library(r.jive)
png("JICOheatmap.png", width=400, height=600)
par(mfrow=c(1,2), mar = rep(0.5, 4))
matrix.J = rbind(ml.iter$J[[3]][hc3$order,hc$order], 
                 matrix(rep(0, 93*10), ncol = 93), 
                 ml.iter$J[[2]][rev(hc2$order),hc$order], 
                 matrix(rep(0, 93*10), ncol = 93), 
                 ml.iter$J[[1]][hc1$order,hc$order])
matrix.I = rbind(ml.iter$I[[3]][hc3$order,hc$order], 
                 matrix(rep(0, 93*10), ncol = 93), 
                 ml.iter$I[[2]][rev(hc2$order),hc$order], 
                 matrix(rep(0, 93*10), ncol = 93), 
                 ml.iter$I[[1]][hc1$order,hc$order])
show.image(matrix.J)
show.image(matrix.I)
dev.off()

# show.image(matrix(c(rep(1, 100), rep(-1, 300)), ncol = 10, byrow = T))
show.image= function (Image, ylab = ""){
  lower = mean(Image) - 3 * sd(Image)
  upper = mean(Image) + 3 * sd(Image)
  Image[Image < lower] = lower
  Image[Image > upper] = upper
  image(x = 1:dim(Image)[2], y = 1:dim(Image)[1], z = t(Image), 
        zlim = c(lower, upper), axes = FALSE, col = bluered(100), 
        xlab = "", ylab = ylab)
}

showpanel <- function(col)
{
  image(z=matrix(1:100, ncol=1), col=col, xaxt="n", yaxt="n" )
}
par(mfrow=c(1,1), mar = c(0, 4, 0, 4), cex = 4)
showpanel(bluered(100))
#title("Title text", adj = 1, line = -3)
#title("Title text", adj = 0, line = -3)
axis(2, at=0, labels= '-', tick = F, las= 1)
axis(4, at=0, labels= '+', tick = F, las= 1)
          