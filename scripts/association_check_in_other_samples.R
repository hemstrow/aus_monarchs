library(snpR);library(data.table);library(ggplot2)
# #===========prep data================
# # 2018 data
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

s <- get.snpR.stats(dat)
s <- s[s$scaffold == "DPSCF300064",]
ggplot(s, aes(x = position, y = -log10(gmmat_pval_likely_pheno))) + geom_point() +
  scale_x_continuous(limits = c(1000000, 1250000))

# pull in pacific data, no minor allele frequency filter
pdat <- readRDS("../F-H_2018/Raw_data/nomaf_snps_snpR.RDS")

# 
# #=============merge===================
# # genotypes
# dat.geno <- as.data.table(cbind(dat@snp.meta, dat))
# pdat.geno <- as.data.table(cbind(pdat@snp.meta, pdat))
# colnames(pdat.geno)[1] <- "scaffold"
# pdat.geno <- pdat.geno[,-".snp.id"]
# dat.geno <- dat.geno[,-".snp.id"]
# mdat.geno <- merge(dat.geno, pdat.geno, by = c("scaffold", "position"), all = T)
# rm(dat.geno, pdat.geno); gc(); gc()
# 
# # sample meta
# dat.sa.meta <- as.data.table(dat@sample.meta)
# pdat.sa.meta <- as.data.table(pdat@sample.meta)
# colnames(dat.sa.meta)[1] <- "samp"
# dat.sa.meta$pop <- "QLD"
# dat.sa.meta$source <- "migration.test"
# pdat.sa.meta$source <- "pacific.samples"
# pdat.sa.meta <- pdat.sa.meta[,-".sample.id"]
# dat.sa.meta <- dat.sa.meta[,-".sample.id"]
# mdat.sa.meta <- merge(dat.sa.meta, pdat.sa.meta, all = T)
# rm(dat.sa.meta, pdat.sa.meta); gc(); gc()
# 
# # snp meta
# mdat.sn.meta <- mdat.geno[,1:2]
# rm(pdat); gc(); gc();
# 
# 
# #==============get the snp of interest===========
# # consult test results, grab top ten
# i.snps <- dplyr::arrange(get.snpR.stats(dat), gmmat_pval_likely_pheno)[1:10,c(3:4, 10)]
# i.snps <- as.data.table(i.snps)
# 
# # grab these sites from main dataset
# i.snps <- merge(i.snps, mdat.geno, by = c("scaffold", "position"), all.x = T, all.y = F)
# 
# # put into snpRdata object, just these snps!
# i.snps <- na.omit(i.snps)
# idat <- import.snpR.data(i.snps[,-c(1:3)], i.snps[,1:3], mdat.sa.meta)
# saveRDS(idat, "results/mon_isnps_snpRdat.RDS") # save
idat <- readRDS("results/mon_isnps_snpRdat.RDS") # read in, skip the prior parts of script!

# check that global majors and minors are the same
idat <- calc_maf(idat)
min.ident.check <- merge(get.snpR.stats(dat), get.snpR.stats(idat), by = c("facet", "subfacet", "scaffold", "position"))
min.ident.check <- min.ident.check[,c("facet", "subfacet", "scaffold", "position", "major.y", "major.x")]

#==============get minor allele frequencies per population==========
idat <- calc_maf(idat, "pop.source")
plot_dat <- get.snpR.stats(idat, "pop.source")

#=============plot==================================================
plot_dat$source <- gsub("^.{4}", "", plot_dat$subfacet)
plot_dat$pop <- substr(plot_dat$subfacet, 1, 3)
plot_dat$snp <- paste0(plot_dat$scaffold, "_", plot_dat$position)
plot_dat$snp <- factor(plot_dat$snp, levels = paste0(idat@snp.meta$scaffold, "_", idat@snp.meta$position)[order(idat@snp.meta$gmmat_pval_likely_pheno)])

mpdat <- plot_dat[,c(2,11,14)]
mpdat$minor <- mpdat$maf
mpdat$major<- 1 - mpdat$maf
mpdat$maf <- NULL
mpdat$n <- plot_dat$maj.count + plot_dat$min.count
mpdat <- reshape2::melt(mpdat, id.vars = c("subfacet", "snp", "n"))
mpdat$subfacet <- gsub(".pacific.samples", "w", mpdat$subfacet)
mpdat$subfacet <- gsub(".migration.test", "t", mpdat$subfacet)
mpdat$snp <- factor(mpdat$snp, levels = paste0(idat@snp.meta$scaffold, "_", idat@snp.meta$position)[order(idat@snp.meta$gmmat_pval_likely_pheno)])

# all pops
ggplot(mpdat, aes(x = variable, y = value)) + 
  geom_bar(stat = "identity") + theme_bw() +
  facet_grid(snp~subfacet)
ggplot(plot_dat, aes(x = pop, y = maf, fill = source)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  facet_grid(snp~., ) + theme_bw()

# only decent sample sizes
ggplot(mpdat[mpdat$n >= 20,], aes(x = variable, y = value)) + 
  geom_bar(stat = "identity") + theme_bw() +
  facet_grid(snp~subfacet)
ggplot(plot_dat[plot_dat$maj.count + plot_dat$min.count >= 20,], aes(x = pop, y = maf, fill = source)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  facet_grid(snp~., ) + theme_bw()

#================tables of genotype counts===============
cbind(idat@facet.meta, idat@geno.tables$gs)[idat@facet.meta$position == 1180485,]
tdat <- subset_snpR_data(idat, facets = "source", subfacets = "migration.test")
tdat <- calc_maf(tdat, "year")
cbind(tdat@facet.meta, tdat@geno.tables$gs)[tdat@facet.meta$position == 1180485,]


id2 <- subset_snpR_data(idat, facets = "source", subfacets = "migration.test")

gp <- run_genomic_prediction(id2, "likely_pheno", 100000, 1000, thin = 200)
rf <- run_random_forest(id2, response = "likely_pheno", mtry = nrow(id2), num.trees = 1000000)

sn <- format_snps(id2, output = "sn", interpolate = F)

sn <- t(sn[,-c(1:3)])
sn <- cbind(id2@sample.meta, sn)
colnames(sn)[13:ncol(sn)] <- paste0("snp_", 1:9)

mod <- glm(as.factor(likely_pheno)~, sn, family = binomial(link = "logit"))
