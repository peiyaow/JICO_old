load("myseeds.RData")
nsim = length(myseeds)
for (i in 1:nsim){
	system(paste0('sbatch -o iter_', myseeds[i], '.out -t 05:00:00 --wrap="Rscript ADNI_iter.R myseed=', myseeds[i], '"'))
#	system(paste0('sbatch -o main.out -t 02:00:00 --wrap="Rscript ADNI.R myseed=', myseeds[i], '"'))
}
