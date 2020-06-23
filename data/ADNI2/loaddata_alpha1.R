load("/nas/longleaf/home/peiyao/alpha/data/ADNI2_clean3.RData")
# load("~/Documents/GitHub/alpha/data/ADNI2_clean3.RData")

selection = (label!=1)&(label!=3)
X = X[selection,]
Y = Y[selection]
label = label[selection]
label = droplevels(label)

BOXPLOT = boxplot(Y~label)

selection = !((label == 0)&(Y%in%BOXPLOT$out[BOXPLOT$group == 1])
              |(label == 2)&(Y%in%BOXPLOT$out[BOXPLOT$group == 2])
              |(label == 4)&(Y%in%BOXPLOT$out[BOXPLOT$group == 3]))

X = X[selection,]
Y = Y[selection]
label = label[selection]



