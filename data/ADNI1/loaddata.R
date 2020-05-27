library(readr)
setwd("~/Documents/GitHub/continuum/data/ADNI1")
t = 1
X0 = as.matrix(read_table("X1a.txt", col_names = F))
Y0 = as.matrix(read_table("Y5T.txt", col_names = F))[, 4+5*(t-1)]
label0 = as.ordered(read_table("label1.txt", col_names = F)$X1)

# id.cut1 = !(label0==4)
# X0 = X0[id.cut1, ]
# Y0 = Y0[id.cut1]
# label0 = label0[id.cut1]
# label0 = droplevels(label0)

# # remove NAs in Y
# id.cut1 = !is.na(Y0)
# X0 = X0[id.cut1, ]
# Y0 = Y0[id.cut1]
# label0 = label0[id.cut1]

# remove negatives in Y
id.cut2 = (Y0 > 0)
X0 = X0[id.cut2,]
Y0 = Y0[id.cut2]
label0 = label0[id.cut2]

# only select MRI
Xs = X0[, 1:93]

X = Xs
Y = Y0
label = label0
