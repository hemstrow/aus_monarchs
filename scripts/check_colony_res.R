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
#pc$all_pairs$sampleID <- dat@sample.meta$SampleID

pc$all_pairs$Sample1 <- dat@sample.meta$SampleID[match(pc$all_pairs$Sample1, colnames(dat))]
pc$all_pairs$Sample2 <- dat@sample.meta$SampleID[match(pc$all_pairs$Sample2, colnames(dat))]


pc$dyads$Sample1 <- dat@sample.meta$SampleID[match(pc$dyads$Sample1, colnames(dat))]
pc$dyads$Sample2 <- dat@sample.meta$SampleID[match(pc$dyads$Sample2, colnames(dat))]

pc$dyads$m1 <- dat@sample.meta$likely_pheno[match(pc$dyads$Sample1,  dat@sample.meta$SampleID)]
pc$dyads$m2 <- dat@sample.meta$likely_pheno[match(pc$dyads$Sample2,  dat@sample.meta$SampleID)]


View(pc$dyads[pc$dyads$type == "FullSib",])

clusters <- readr::read_table2("colony/colony_input.BestConfig_Ordered", col_names = T)
clusters$OffspringID <- dat@sample.meta$SampleID[match(clusters$OffspringID, colnames(dat))]
clusters$m <- dat@sample.meta$likely_pheno[match(clusters$OffspringID,  dat@sample.meta$SampleID)]


keeps <- clusters$OffspringID[clusters$m == "M"]
res <- clusters[clusters$m == "R",]
cap <- 96
while(length(keeps) < cap){
  uclusts <- unique(res$ClusterIndex)
  tkeep <- character()
  for(i in 1:length(uclusts)){
    tres <- res[res$ClusterIndex == uclusts[i],]
    tkeep <- c(tkeep, tres[sample(nrow(tres), 1), 1])
  }
  if(length(keeps) + length(uclusts) > cap){
    tkeep <- tkeep[sample(length(tkeep), cap - length(keeps), replace = F)]
  }
  
  keeps <- c(keeps, unlist(tkeep))
  
  res <- res[-which(res$OffspringID %in% tkeep),] 
}

check <- data.frame(keeps = keeps, fam = clusters$ClusterIndex[match(keeps, clusters$OffspringID)])
check$m <- clusters$m[match(keeps, clusters$OffspringID)]
check_tab <- table(check[,2:3])
margin.table(check_tab, 1)

library(ggplot2)
ggplot(pc$dyads, aes(x = as.factor(Sample1), y = as.factor(Sample2), fill = type, color = type)) + theme_bw() + geom_tile() +
  theme(panel.grid = element_blank())




