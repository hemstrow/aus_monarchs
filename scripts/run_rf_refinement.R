# run parameters
trim <- .9
response <- "likely_pheno"
mtry.scale <- 1 # initial mtry
num.trees <- 50000 # number of trees
par <- F



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


# run the first random forest
rf <- run_random_forest(rdat, response = "likely_pheno", mtry = nrow(rdat)*mtry.scale, num.trees = num.trees)




# run the refinement
best.imp <- quantile(abs(rf$models$.base_.base$model$variable.importance), trim)
best.imp <- which(rf$models$.base_.base$model$variable.importance >= best.imp[1])
rdat <- subset_snpR_data(rdat, best.imp)
rf <- run_random_forest(rdat, response = response, mtry = nrow(rdat)*mtry.scale, num.trees = num.trees, 
                        par = par, importance = "permutation", pvals = F)



saveRDS(rf, paste0(outfile, "_model_", run, ".RDS"))

# prediction, need to re-run the model since the impurity_corrected importance can apparently cause issues. Run with permutation.
kept.snps <- rf$data@snp.meta$.snp.id
kept.snps <- match(kept.snps, dat@snp.meta$.snp.id)

sn <- format_snps(dat, "sn")
sn <- sn[,-c(1:3)]
sn <- sn[kept.snps, cross_sample]
sn <- as.data.frame(t(sn))
colnames(sn) <- rf$models$.base_.base$model$forest$independent.variable.names
out <- predict(predict.rf$models$.base_.base$model, sn)

write.table(data.frame(sample = cross_sample, phenotype = dat@sample.meta[cross_sample,response], prediction = out$predictions),
            paste0(outfile, "_prediction_", run, ".txt"), quote = F, col.names = F, row.names = F)
