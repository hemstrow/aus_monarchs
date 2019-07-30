library(snpR)

# import data
sample.meta <- read.table("data/sample_metadata.txt", header = T, stringsAsFactors = F)
genos <- read.table("data/genotypes.geno", stringsAsFactors = F)
snp.meta <- genos[,1:2]
colnames(snp.meta) <- c("scaffold", "position")

dat <- import.snpR.data(genos[,-c(1:2)], snp.meta, sample.meta)

# filter, subset, and set seed.
miss.samples <- which(is.na(dat@sample.meta$likely_pheno))
dat <- subset_snpR_data(dat, samps = (1:ncol(dat))[-miss.samples])
dat <- filter_snps(dat, maf = 0.10, hf_hets = 0.55, min_ind = 0.5, min_loci = 0.4)
keep.snps <- sample(nrow(dat), 30000)
dat <- subset_snpR_data(dat, snps = keep.snps)
write.table(matrix(keep.snps, ncol = 1), file = "colony/colony_snps.txt", quote = F, col.names = F, row.names = F)
seed <- sample(10000, 1)
write(seed, "colony/colony_seed.txt")

# run colony
col_res <- run_colony(dat, colony_path = "~/bin/colony2s.ifort.out", seed = seed, precision = 3, run_length = 3)

# save results
saveRDS(list(col_res = col_res, data = dat), "colony_R_outputs.RDS")
