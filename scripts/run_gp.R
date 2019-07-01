library(snpR)

# import data
sample.meta <- read.table("data/sample_metadata.txt", header = T, stringsAsFactors = F)
genos <- read.table("data/genotypes.geno", stringsAsFactors = F)
snp.meta <- genos[,1:2]
colnames(snp.meta) <- c("scaffold", "position")

dat <- import.snpR.data(genos[,-c(1:2)], snp.meta, sample.meta)

# filter
fdat <- filter_snps(dat, maf = 0.05, hf_hets = 0.55, min_ind = 0.5, min_loci = 0.4)

miss.samples <- which(is.na(fdat@sample.meta$likely_pheno))

#p1 <- plot_clusters(fdat, "likely_pheno")

# run
fdat2 <- subset_snpR_data(fdat, samps = (1:ncol(fdat))[-miss.samples])
gp <- run_genomic_prediction(fdat2, "quant_pheno", iterations = 100000, burn_in = 1000, thin = 200)

cv <- cross_validate_genomic_prediction(fdat2, "quant_pheno", iterations = 100000, burn_in = 1000, thin = 200, cross_percentage = .9)

saveRDS(list(gp = gp, cv = cv), "trial_gp_cv_100k.RDS")
