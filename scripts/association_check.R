library(snpR); library(qqman)
# prep data
sample.meta <- read.table("data/sample_metadata.txt", header = T, stringsAsFactors = F)
genos <- read.table("data/v4/genotypes.geno", stringsAsFactors = F)
snp.meta <- genos[,1:2]
colnames(snp.meta) <- c("scaffold", "position")
genos <- genos[,-c(1:2)]
genos <- genos[,-which(is.na(sample.meta$relaxed_pheno))]
sample.meta <- sample.meta[-which(is.na(sample.meta$relaxed_pheno)),]
dat <- import.snpR.data(genos, snp.meta, sample.meta)
# dat <- dat[,-which(sample.meta(dat)$year != 2018)] # note: same locus pops out if you only use my samples, but is stronger with both. Good sign!
# dat <- filter_snps(dat, maf = 0.05, hf_hets = .55)
dat <- filter_snps(dat, maf = 0.05, hwe = 0.000001)

# dat <- calc_pairwise_fst(dat, "year", boot = 100, boot_par = 3)

# calc_association(dat, response = "likely_pheno", method = "chisq")

# run assocition test
#dat <- calc_association(dat, response = "relaxed_pheno", maxiter = 1000)
dat <- calc_association(dat, response = "likely_pheno", maxiter = 1000)
#dat <- calc_association(dat, response = "strict_pheno", maxiter = 1000)
#dat <- calc_association(dat, response = "quant_pheno", maxiter = 1000)
#dat <- calc_association(dat, response = "likely_pheno", method = "armitage")

# plot
p <- plot_manhattan(dat, "gmmat_pval_likely_pheno", chr = "scaffold", significant = 0.00001, suggestive = 0.0001, log.p = T)
p$plot + ggplot2::ylab("-log10(p-value)") + ggplot2::xlab("Position") + ggplot2::theme(strip.text = ggplot2::element_blank())
#p2 <- plot_manhattan(dat, "p_armitage_likely_pheno", chr = "scaffold", significant = 0.00001, suggestive = 0.0001, log.p = T)

#qq plots
qq(get.snpR.stats(dat)$gmmat_pval_likely_pheno)
qq(get.snpR.stats(dat)$p_armitage_likely_pheno)

dat <- calc_pairwise_fst(dat, c("year", "likely_pheno", "year.likely_pheno"))
dat <- calc_pairwise_fst(dat, "year", boot = 1000)

get.snpR.stats(dat, c("year", "likely_pheno", "year.likely_pheno"), "fst")$weighted.means










# from angsd
adat <- read.table("results/association_out_likely.lrt0", header = T)
adat$pval <- pchisq(adat$LRT, 1, lower.tail = FALSE)
adat$lpval <- -log(adat$pval)

q <- qqman::qq(adat$pval)
p3 <- plot_manhattan(adat, "pval", chr = "Chromosome", bp = "Position", log.p = T, significant = 0.000000005, suggestive = 0.00000005)

# save PDFs
pdf("association_test_gmmat.pdf")
p$plot
pdf("association_test_gmmat_qq.pdf")
qq(get.snpR.stats(dat)$gmmat_pval_likely_pheno)
pdf("association_test_armitage.pdf")
p2$plot
pdf("association_test_armitage_qq.pdf")
qq(get.snpR.stats(dat)$p_armitage_likely_pheno)
pdf("association_test_LRT.pdf")
p3$plot
pdf("association_test_LRT_qq.pdf")
qqman::qq(adat$pval)
dev.off(); dev.off();dev.off();dev.off();dev.off();dev.off();dev.off();


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
