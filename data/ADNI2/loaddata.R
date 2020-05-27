library(readr)
# setwd("~/Documents/GitHub/continuum/data/ADNI2")
# setwd("/nas/longleaf/home/peiyao/continuum")

X = as.matrix(read_table("/nas/longleaf/home/peiyao/continuum/data/ADNI2/X2.txt", col_names = F))
Y = as.matrix(read_table("/nas/longleaf/home/peiyao/continuum/data/ADNI2/Y2.txt", col_names = F))
label = as.ordered(read_table("/nas/longleaf/home/peiyao/continuum/data/ADNI2/label2.txt", col_names = F)$X1)

X = X[!is.na(Y),]
label = label[!is.na(Y)]
Y = Y[!is.na(Y)]

X = X[Y > 0,]
label = label[Y > 0]
Y = Y[Y > 0]


#boxplot(Y~label)

selection = !((label==0)&(Y > 11))
X = X[selection,]
Y = Y[selection]
label = label[selection]
#boxplot(Y~label)

selection = !((label==2)&(Y > 16))
X = X[selection,]
Y = Y[selection]
label = label[selection]
#boxplot(Y~label)

selection = !((label==3)&(Y > 24))
X = X[selection,]
Y = Y[selection]
label = label[selection]
#boxplot(Y~label)

# selection = !((label==4)&(Y > 32))
# X = X[selection,]
# Y = Y[selection]
# label = label[selection]
# boxplot(Y~label)

selection = (label!=1)&(label!=3)
X = X[selection,]
Y = Y[selection]
label = label[selection]

label = droplevels(label)
# Y = log2(Y+1)
# boxplot(log2(Y+1)~label)



