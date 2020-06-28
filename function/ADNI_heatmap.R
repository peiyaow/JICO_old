X.list = lapply(label.level, function(l) X[label == l,])
Y.list = lapply(label.level, function(l) matrix(Y[label == l]))
# rank_iter = as.matrix(rank_iter)
ml.iter = continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 200,
                                    gam = 1, 
                                    rankJ = 1, 
                                    rankA = c(1,1,1),
                                    orthIndiv = F)

# # library(r.jive)
# # ml.iter$J
# par(mfcol = c(3,2), mar = rep(0.5, 4))
# show.image((ml.iter$J[[1]]))
# show.image((ml.iter$J[[2]]))
# show.image((ml.iter$J[[3]]))
# 
# show.image((ml.iter$I[[1]]))
# show.image((ml.iter$I[[2]]))
# show.image((ml.iter$I[[3]]))
# 
# dev.off()
# 
# 
# par(mfrow=c(1,2))
# show.image((do.call(rbind, rev(ml.iter$J))))
# # abline(h = 145)
# # abline(h = 145+171)
# # abline(h = 145+171+178)
# show.image((do.call(rbind, rev(ml.iter$I))))
# #show.image((ml.iter$R))
# dev.off()

png("JICOheatmap.png", width=400, height=600)
par(mfrow=c(1,2), mar = rep(0.5, 4))
show.image(rbind(ml.iter$J[[3]], matrix(rep(0, 93*10), ncol = 93), ml.iter$J[[2]], matrix(rep(0, 93*10), ncol = 93), ml.iter$J[[1]]))
show.image(rbind(ml.iter$I[[3]], matrix(rep(0, 93*10), ncol = 93), ml.iter$I[[2]], matrix(rep(0, 93*10), ncol = 93), ml.iter$I[[1]]))
dev.off()

# show.image(matrix(c(rep(1, 100), rep(-1, 300)), ncol = 10, byrow = T))
