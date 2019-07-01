#ybins
sample.meta <- read.table("data/sample_metadata.txt", header = T, stringsAsFactors = F)
osm <- sample.meta
bad <- which(is.na(sample.meta$SampleID))
sample.meta <- sample.meta[-bad,]
ybin <- matrix(NA, nrow(sample.meta), 4)
colnames(ybin) <- c("relaxed", "strict", "likely", "quant")
ybin[,"relaxed"] <- ifelse(sample.meta$relaxed_pheno == "M", 0, 1)
ybin[,"strict"] <- ifelse(sample.meta$strict_pheno == "M", 0, 1)
ybin[,"likely"] <- ifelse(sample.meta$likely_pheno == "M", 0, 1)
ybin[,"quant"] <- sample.meta$quant_pheno


write.table(ybin[,1], "data/relaxed_ybin.txt", quote = F, col.names = F, row.names = F)
write.table(ybin[,2], "data/strict_ybin.txt", quote = F, col.names = F, row.names = F)
write.table(ybin[,3], "data/likely_ybin.txt", quote = F, col.names = F, row.names = F)
write.table(ybin[,4], "data/quant_ybin.txt", quote = F, col.names = F, row.names = F)

# bamlist
bl <- read.table("data/gbamlist.txt", stringsAsFactors = F)
write.table(bl[-bad,], "data/gbamlist_minus_bads.txt", quote = F, col.names = F, row.names = F)
