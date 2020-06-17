L = 50
n1 = 35
n2 = 34
n3 = 29
n = n1 + n2 + n3
result = result_iter
result.mean = matrix(apply(result, 2, mean)[-1], ncol = G)%*%c(n1, n2, n3)/n

RANK[3,]
a
# gam.list
1:5

nrow(result.mean)
plot(a, result.mean[1:51], type = "l")
lines(a, result.mean[51+1:51], type = "l")
lines(a, result.mean[51*2+1:51], type = "l")
lines(a, result.mean[51*3+1:51], type = "l")
lines(a, result.mean[51*4+1:51], type = "l")

#result.mean[6:15]
plot(a, result.mean[51*5 + 1:51], type = "l", ylim = c(min(result.mean[(51*5+1):(51*15)]), max(result.mean[(51*5+1):(51*15)])))
lines(a, result.mean[51*6+1:51], type = "l")
lines(a, result.mean[51*7+1:51], type = "l")
lines(a, result.mean[51*8+1:51], type = "l")
lines(a, result.mean[51*9+1:51], type = "l")
lines(a, result.mean[51*10+1:51], type = "l")
lines(a, result.mean[51*11+1:51], type = "l")
lines(a, result.mean[51*12+1:51], type = "l")
lines(a, result.mean[51*13+1:51], type = "l")
lines(a, result.mean[51*14+1:51], type = "l")


16:35
plot(a, result.mean[51*15 + 1:51], type = "l", ylim = c(min(result.mean[(51*15+1):(51*35)]), max(result.mean[(51*15+1):(51*35)])))
lines(a, result.mean[51*16+1:51], type = "l")
lines(a, result.mean[51*17+1:51], type = "l")
lines(a, result.mean[51*18+1:51], type = "l")
lines(a, result.mean[51*19+1:51], type = "l")
lines(a, result.mean[51*20+1:51], type = "l")
lines(a, result.mean[51*21+1:51], type = "l")
lines(a, result.mean[51*22+1:51], type = "l")
lines(a, result.mean[51*23+1:51], type = "l")
lines(a, result.mean[51*24+1:51], type = "l")
lines(a, result.mean[51*25+1:51], type = "l")
lines(a, result.mean[51*26+1:51], type = "l")
lines(a, result.mean[51*27+1:51], type = "l")
lines(a, result.mean[51*28+1:51], type = "l")
lines(a, result.mean[51*29+1:51], type = "l")
lines(a, result.mean[51*30+1:51], type = "l")
lines(a, result.mean[51*31+1:51], type = "l")
lines(a, result.mean[51*32+1:51], type = "l")
lines(a, result.mean[51*33+1:51], type = "l")
lines(a, result.mean[51*34+1:51], type = "l")
lines(a, result.mean[51*35+1:51], type = "l")

