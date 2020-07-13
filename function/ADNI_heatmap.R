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
show.image(rbind(ml.iter$J[[3]][hc3$order,hc$order], matrix(rep(0, 93*10), ncol = 93), ml.iter$J[[2]][rev(hc2$order),hc$order], matrix(rep(0, 93*10), ncol = 93), ml.iter$J[[1]][hc1$order,hc$order]))
show.image(rbind(ml.iter$I[[3]][hc3$order,hc$order], matrix(rep(0, 93*10), ncol = 93), ml.iter$I[[2]][rev(hc2$order),hc$order], matrix(rep(0, 93*10), ncol = 93), ml.iter$I[[1]][hc1$order,hc$order]))
dev.off()

# show.image(matrix(c(rep(1, 100), rep(-1, 300)), ncol = 10, byrow = T))
