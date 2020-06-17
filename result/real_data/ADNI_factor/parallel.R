load("myseeds.RData")
nsim = length(myseeds)
for (i in 1:nsim){
#	system(paste0('sbatch -o main_', myseeds[i], '.out -t 05:00:00 --wrap="Rscript ADNI.R myseed=', myseeds[i], '"'))
	system(paste0('sbatch -o main.out -t 00:30:00 --wrap="Rscript ADNI_factor.R myseed=', myseeds[i], '"'))
}
