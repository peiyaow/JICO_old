library(readr)
plotPC = function(X){
  X.mean = apply(X, 2, mean)
  X = sweep(X, 2, X.mean)
  n = nrow(X)
  PCA.res = eigen(X%*%t(X)/n)
  par(mfrow = c(3,1))
  plot(PCA.res$vectors[,1]*sqrt(n), PCA.res$vectors[,2]*sqrt(n))
  hehe = boxplot(PCA.res$vectors[,1]*sqrt(n))$out
  haha = boxplot(PCA.res$vectors[,2]*sqrt(n))$out
  par(mfrow=c(1,1))
  return(list(!(PCA.res$vectors[,1]*sqrt(n))%in%hehe, !(PCA.res$vectors[,2]*sqrt(n))%in%haha))
}

setwd("~/Documents/GitHub/continuum/data/ADNI2")
X = as.matrix(read_table("X2.txt", col_names = F))
Y = as.matrix(read_table("Y2.txt", col_names = F))
label = as.ordered(read_table("label2.txt", col_names = F)$X1)

X = X[!is.na(Y),]
label = label[!is.na(Y)]
Y = Y[!is.na(Y)]

X = X[Y > 0,]
label = label[Y > 0]
Y = Y[Y > 0]

selection = (label!=1)&(label!=3)
X = X[selection,]
Y = Y[selection]
label = label[selection]
label = droplevels(label)

# remove first factor outliers
n = nrow(X)
selection = c(((1:n)[label==0])[plotPC(X[label==0,])[[1]]],
              ((1:n)[label==2])[plotPC(X[label==2,])[[1]]],
              ((1:n)[label==4])[plotPC(X[label==4,])[[1]]])
X = X[selection,]
Y = Y[selection]
label = label[selection]

# remove Y outliers
BOXPLOT = boxplot(Y~label)
selection = !((label == 0)&(Y%in%BOXPLOT$out[BOXPLOT$group == 1])
              |(label == 2)&(Y%in%BOXPLOT$out[BOXPLOT$group == 2])
              |(label == 4)&(Y%in%BOXPLOT$out[BOXPLOT$group == 3]))
  
X = X[selection,]
Y = Y[selection]
label = label[selection]

boxplot(Y~label)
