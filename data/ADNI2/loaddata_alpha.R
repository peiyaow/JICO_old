load("/nas/longleaf/home/peiyao/alpha/data/ADNI2_clean3.RData")
selection = (label!=1)&(label!=3)
X = X[selection,]
Y = Y[selection]
label = label[selection]

label = droplevels(label)
