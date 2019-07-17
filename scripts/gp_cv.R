# runs cross validation on filtered snps. Arguments are: 1) which sample(s) should be used for cross validation? 2) what is the outfile name?

# get arguments, define paths, library snpR
args <- commandArgs(TRUE)

.libPaths(c(.libPaths(), "/home/hemstrow/R/x86_64-pc-linux-gnu-library/3.6", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library", "/share/apps/rmodules"))
library(snpR)

cross_samples <- as.numeric(as.character(args[1]))
outfile <- as.character(args[2])


# read in data
sample.meta <- read.table("~/monarch/aus_monarchs/aus_monarchs/data/sample_metadata.txt", header = T, stringsAsFactors = F)
genos <- read.table("~/monarch/aus_monarchs/aus_monarchs/data/sample_metadata.txt", stringsAsFactors = F)
snp.meta <- genos[,1:2]
colnames(snp.meta) <- c("scaffold", "position")
dat <- import.snpR.data(genos[,-c(1:2)], snp.meta, sample.meta)
miss.samples <- which(is.na(dat@sample.meta$likely_pheno))
dat <- subset_snpR_data(dat, samps = (1:ncol(dat))[-miss.samples])
dat <- filter_snps(dat, maf = 0.05, hf_hets = 0.55, min_ind = 0.5, min_loci = 0.4)


# run the cross_validation
cv_out <- cross_validate_genomic_prediction(dat, "likely_pheno",
                                            iterations = 1000000,
                                            burn_in = 3000,
                                            thin = 300, cross_samples = cross_samples, plot = F)


# save output
out <- cbind(sampleID = rownames(cv_out$comparison), samplenum = paste0(cross_samples, collapse = "_"), cv_out$comparison, h2 = cv_out$model$h2)

write.table(out, outfile, quote = F, col.names = F, row.names = F, append = T)

outRDS <- paste0(gsub(".txt", "_", outfile), cross_samples, ".RDS")

saveRDS(cv_out, outRDS)
