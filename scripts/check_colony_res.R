library(snpR)
# import data
sample.meta <- read.table("data/sample_metadata.txt", header = T, stringsAsFactors = F)
genos <- read.table("data/genotypes.geno", stringsAsFactors = F)
snp.meta <- genos[,1:2]
colnames(snp.meta) <- c("scaffold", "position")

dat <- import.snpR.data(genos[,-c(1:2)], snp.meta, sample.meta)


# filter
miss.samples <- which(is.na(dat@sample.meta$likely_pheno))
dat <- subset_snpR_data(dat, samps = (1:ncol(dat))[-miss.samples])
dat <- filter_snps(dat, maf = 0.10, hf_hets = 0.55, min_ind = 0.5, min_loci = 0.4)

csnps <- read.table("colony/colony_snps.txt")
csnps <- as.numeric(csnps[,1])
csnps <- sort(csnps)
dat <- subset_snpR_data(dat, csnps)


# pull in colony results
pc <- parse_colony("colony_input", x = dat, path = "colony/")
pc$all_pairs$sampleID <- dat@sample.meta$SampleID

pc$all_pairs$Sample1 <- dat@sample.meta$SampleID[match(pc$all_pairs$Sample1, colnames(dat))]
pc$all_pairs$Sample2 <- dat@sample.meta$SampleID[match(pc$all_pairs$Sample2, colnames(dat))]


pc$dyads$Sample1 <- dat@sample.meta$SampleID[match(pc$dyads$Sample1, colnames(dat))]
pc$dyads$Sample2 <- dat@sample.meta$SampleID[match(pc$dyads$Sample2, colnames(dat))]


library(ggplot2)
ggplot(pc$dyads, aes(x = as.factor(Sample1), y = as.factor(Sample2), fill = type, color = type)) + theme_bw() + geom_tile()
