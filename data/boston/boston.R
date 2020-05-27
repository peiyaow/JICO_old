library(mlbench)

data(BostonHousing2)
num = c("lon", "lat", "crim", "zn", "indus", "chas", "nox", "rm", "age", "dis", "tax", "ptratio", "b", "lstat")  
cat = c("town", "chas", "rad")
X.num = BostonHousing2[num]
X.cat = apply(BostonHousing2[cat], 2, as.factor)

Y = BostonHousing2["cmedv"]
X = apply(as.matrix(X.num[, -(1:2)]), 2, as.numeric)
label = X.cat[, 3]
summary(X.cat)

X = X[label %in% c("24", " 5", " 4"),]
Y = as.matrix(Y[label %in% c("24", " 5", " 4"),])
label = as.ordered(label[label %in% c("24", " 5", " 4")])
