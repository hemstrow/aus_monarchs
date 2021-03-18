library(snpR)

# import data
sample.meta <- readr::read_delim("data/sample_metadata.txt", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
genos <- read.table("data/v4/genotypes.geno", stringsAsFactors = F)
snp.meta <- genos[,1:2]
colnames(snp.meta) <- c("scaffold", "position")

dat <- import.snpR.data(genos[,-c(1:2)], snp.meta, sample.meta)

# filter, subset, and set seed.
miss.samples <- which(is.na(dat@sample.meta$likely_pheno))
dat <- subset_snpR_data(dat, samps = (1:ncol(dat))[-miss.samples])
dat <- filter_snps(dat, maf = 0.10, HWE = 0.000001, min_ind = .5)
keep.snps <- sample(nrow(dat), 30000)
dat <- subset_snpR_data(dat, snps = keep.snps)
write.table(matrix(keep.snps, ncol = 1), file = "colony/colony_snps.txt", quote = F, col.names = F, row.names = F)
seed <- sample(10000, 1)
write(seed, "colony/colony_seed.txt")

# run colony
col_res <- run_colony(dat, colony_path = "D:/usr/bin/Colony/colony2s.exe", seed = seed, precision = 3, run_length = 3)

# save results
saveRDS(list(col_res = col_res, data = dat), "colony_R_outputs.RDS")
