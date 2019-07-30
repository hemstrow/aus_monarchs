# define paths, library snpR
.libPaths(c(.libPaths(), "/home/hemstrow/R/x86_64-pc-linux-gnu-library/3.6", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library", "/share/apps/rmodules"))
library(snpR); library(ranger);

args <- commandArgs(TRUE)
run <- as.character(args[1])
outfile <- as.character(args[2])
cross_sample <- as.numeric(as.character(args[3]))

# prep data
sample.meta <- read.table("~/monarch/aus_monarchs/aus_monarchs/data/sample_metadata.txt", header = T, stringsAsFactors = F)
genos <- read.table("~/monarch/aus_monarchs/aus_monarchs/data/genotypes.geno", stringsAsFactors = F)
snp.meta <- genos[,1:2]
colnames(snp.meta) <- c("scaffold", "position")
dat <- import.snpR.data(genos[,-c(1:2)], snp.meta, sample.meta)
miss.samples <- which(is.na(dat@sample.meta$likely_pheno))
dat <- subset_snpR_data(dat, samps = (1:ncol(dat))[-miss.samples])
dat <- filter_snps(dat, maf = 0.05, hf_hets = 0.55, min_ind = 0.5, min_loci = 0.4)

## remove the cross sample
rdat <- subset_snpR_data(dat, samps = (1:ncol(dat))[-cross_sample])

# run parameters
trim_cuttoffs <- c(1000, 200) # above the first level, will trim at the first percentage, below the last, will trim a single at a time
trim <- c(.9, .1) # trim at .9 before 1000 SNPs, at .1 before 200, and one per run below this.
response <- "likely_pheno"
num.trees <- 50000 # number of trees
init.mtry <- nrow(dat) # initial mtry
init.num.trees <- 50000 # number of trees in initial run. Probably doesn't need to be huge!
par <- 11

# run the first random forest
rf <- run_random_forest(rdat, response = "likely_pheno", mtry = init.mtry, num.trees = init.num.trees)

# run the refinement
refined_model <- refine_rf(rf = rf, 
                           response = response, 
                           trim_cuttoffs = trim_cuttoffs, 
                           num.trees = num.trees,
                           trim = trim,
                           par = par)

# save outputs
error.delta <- refined_model$error_delta
error.delta <- cbind(run = run, error.delta)
write.table(error.delta, paste0(outfile, "_delta_", run, ".txt"), quote = F, row.names = F, col.names = F)
saveRDS(refined_model, paste0(outfile, "_model_", run, ".RDS"))

# prediction, need to re-run the model since the impurity_corrected importance can apparently cause issues. Run with permutation.
kept.snps <- refined_model$best_model$data@snp.meta$.snp.id
kept.snps <- match(kept.snps, dat@snp.meta$.snp.id)

predict.rf <- run_random_forest(refined_model$best_model$data, response = response, mtry = nrow(refined_model$best_model$data), num.trees = num.trees, importance = "permutation", pvals = F)
sn <- format_snps(dat, "sn")
sn <- sn[,-c(1:3)]
sn <- sn[kept.snps, cross_sample]
sn <- as.data.frame(t(sn))
colnames(sn) <- refined_model$best_model$models$.base_.base$model$forest$independent.variable.names
out <- predict(predict.rf$models$.base_.base$model, sn)

write.table(data.frame(sample = cross_sample, phenotype = dat@sample.meta[cross_sample,response], prediction = out$predictions),
            paste0(outfile, "_prediction_", run, ".txt"), quote = F, col.names = F, row.names = F)
