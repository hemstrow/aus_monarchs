# define paths, library snpR
.libPaths(c(.libPaths(), "/home/hemstrow/R/x86_64-pc-linux-gnu-library/3.6", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library", "/share/apps/rmodules"))
library(snpR); source("scripts/rf_refinement.R")

run <- as.character(args[1])
outfile <- as.character(args[2])

# prep data
sample.meta <- read.table("data/sample_metadata.txt", header = T, stringsAsFactors = F)
genos <- read.table("data/genotypes.geno", stringsAsFactors = F)
snp.meta <- genos[,1:2]
colnames(snp.meta) <- c("scaffold", "position")
dat <- import.snpR.data(genos[,-c(1:2)], snp.meta, sample.meta)
miss.samples <- which(is.na(dat@sample.meta$likely_pheno))
dat <- subset_snpR_data(dat, samps = (1:ncol(dat))[-miss.samples])
dat <- filter_snps(dat, maf = 0.05, hf_hets = 0.55, min_ind = 0.5, min_loci = 0.4)


# run parameters
trim_cuttoffs <- c(1000, 200) # above the first level, will trim at the first percentage, below the last, will trim a single at a time
trim <- c(.9, .1) # trim at .9 before 1000 SNPs, at .1 before 200, and one per run below this.
response <- "likely_pheno"
search_cuttoff <- 1 # keep searching at least until this many SNPs are left
num.trees <- 50000 # number of trees
init.mtry <- nrow(dat) # initial mtry
init.num.trees <- 10000 # number of trees in initial run. Probably doesn't need to be huge!
par <- FALSE

# run the first random forest
rf <- run_random_forest(dat, response = "likely_pheno", mtry = init.mtry, num.trees = init.num.trees)

# run the refinement
refined_model <- refine_rf(rf = rf, 
                           response = response, 
                           trim_cuttoffs = trim_cuttoffs, 
                           num.trees = num.trees,
                           trim = trim,
                           seach_cuttoff = search_cuttoff,
                           par = par)

# save outputs
error.delta <- refined_model$error_delta
error.delta <- cbind(run = run, error.delta)
write.table(error.delta, paste0(outfile, "_delta_", run, ".txt"), quote = F, row.names = F, col.names = F)
saveRDS(refined_model, paste0(outfile, "_model_", run, ".RDS"))