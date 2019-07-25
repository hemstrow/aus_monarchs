library(snpR); library(qqman)
# prep data
sample.meta <- read.table("data/sample_metadata.txt", header = T, stringsAsFactors = F)
genos <- read.table("data/genotypes.geno", stringsAsFactors = F)
snp.meta <- genos[,1:2]
colnames(snp.meta) <- c("scaffold", "position")
dat <- import.snpR.data(genos[,-c(1:2)], snp.meta, sample.meta)
miss.samples <- which(is.na(dat@sample.meta$likely_pheno))
dat <- subset_snpR_data(dat, samps = (1:ncol(dat))[-miss.samples])
dat <- filter_snps(dat, maf = 0.05, hf_hets = 0.55, min_ind = 0.5, min_loci = 0.4)


# run assocition test
#dat <- calc_association(dat, response = "relaxed_pheno", maxiter = 1000)
dat <- calc_association(dat, response = "likely_pheno", maxiter = 1000)
#dat <- calc_association(dat, response = "strict_pheno", maxiter = 1000)
#dat <- calc_association(dat, response = "quant_pheno", maxiter = 1000)
dat <- calc_association(dat, response = "likely_pheno", method = "armitage")

# plot
p <- plot_manhattan(dat, "gmmat_pval_likely_pheno", chr = "scaffold", significant = 0.00001, suggestive = 0.0001, log.p = T)
p2 <- plot_manhattan(dat, "p_armitage_likely_pheno", chr = "scaffold", significant = 0.00001, suggestive = 0.0001, log.p = T)

#qq plots
qq(get.snpR.stats(x)$gmmat_pval_likely_pheno)
qq(get.snpR.stats(x)$p_armitage_likely_pheno)


# from angsd
adat <- read.table("results/association_out_likely.lrt0", header = T)
adat$pval <- pchisq(adat$LRT, 1, lower.tail = FALSE)
adat$lpval <- -log(adat$pval)

qqman::qq(adat$pval)
plot_manhattan(adat, "pval", chr = "Chromosome", bp = "Position", log.p = T, significant = 0.00000005, suggestive = 0.0000005)


# gp
gp <- run_genomic_prediction(dat, "likely_pheno", iterations = 10000, burn_in = 1000, thin = 100)
cvgp <- cross_validate_genomic_prediction(dat, "likely_pheno")
ggplot(gp$predictions, aes(x = phenotype, y = predicted_BV)) + geom_point() + theme_bw()

gp.effects <- cbind(dat@snp.meta, effect = gp$model$ETA[[1]]$b)
gp.effects$effect <- abs(gp.effects$effect)
stats <- get.snpR.stats(dat)
stats <- dplyr::arrange(stats, gmmat_pval_likely_pheno)
hsnps <- head(stats$.snp.id, 10)
hsnps <- which(gp.effects$.snp.id %in% hsnps)
plot_manhattan(gp.effects, "effect", chr = "scaffold", highlight = hsnps)


# rf
rf <- run_random_forest(dat, response = "likely_pheno", num.trees = 500000, par = 10) # probably want to try with higher mtry, and could purge out very unimportant markers as well.
p <- plot_manhattan(rf$data, "likely_pheno_RF_importance", chr = "scaffold", abs = T, significant = 0.02)
# note, the gene on scaff 64 (DPOGS208560) has a hit to a silkworm gene that seems to be related to molt-intermolt transition.