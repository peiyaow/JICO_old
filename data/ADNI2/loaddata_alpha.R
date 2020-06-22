load("/nas/longleaf/home/peiyao/alpha/data/ADNI2_clean3.RData")
# load("~/Documents/GitHub/alpha/data/ADNI2_clean3.RData")

selection = (label!=1)&(label!=3)
X = X[selection,]
Y = Y[selection]
label = label[selection]

label = droplevels(label)

selection = !((label == 2) & (Y >=16))
X = X[selection,]
Y = Y[selection]
label = label[selection]

selection = !((label == 4) & (Y >= 36))
X = X[selection,]
Y = Y[selection]
label = label[selection]
# boxplot(Y~label)




