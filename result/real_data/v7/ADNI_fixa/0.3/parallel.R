load("myseeds.RData")
nsim = length(myseeds)
for (i in 1:nsim){
	system(paste0('sbatch -o main_', myseeds[i], '.out -t 12:00:00 --mem=16g --wrap="Rscript ADNI_fixa.R myseed=', myseeds[i], '"'))
#	system(paste0('sbatch -o main.out -t 02:00:00 --wrap="Rscript ADNI.R myseed=', myseeds[i], '"'))
}
