library(readr)

# import data
barcodes <- read_delim("data/sample_IDs_barcodes.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
meta18 <- read_delim("data/2018_mon_meta.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
meta16 <- read_delim("data/2016_mon_meta.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
bamlist <- read.table("data/gbamlist.txt", stringsAsFactors = F)

# condense phenotype calls
# strict: yolked at all means resident
# relaxed: only chorionated are resident
# likely: best guess
# quantitative: 0 if no yolk, 1 if yolk, 2 if chorionated
meta16$strict_pheno <- ifelse(meta16$Yolked_Oocytes == 0, "M", "R")
meta16$likely_pheno <- meta16$strict_pheno
meta16$relaxed_pheno <- ifelse(meta16$Mature_Oocytes == 0, "M", "R")
meta16$quant_pheno <- ifelse(meta16$Mature_Oocytes != 0, 2, ifelse(meta16$Yolked_Oocytes != 0, 1, 0))

strict_migrant_phenos <- c("I", "Y(barely)", "Y(b)", "N")
relaxed_migrant_phenos <- c("I", "Y(barely)", "Y(b)", "Y", "Y(s)", "Y(m)", "Y(slightly)", "N")
meta18$strict_pheno <- ifelse(meta18$Phenotype %in% strict_migrant_phenos, "M", "R")
meta18$likely_pheno <- meta18$Pheno.Call
meta18$relaxed_pheno <- ifelse(meta18$Phenotype %in% relaxed_migrant_phenos, "M", "R")
yolked_cases <- c("Y(b)", "Y", "Y(mostly)", "Y(barely)", "Y(s)", "Y(m)", "Y(slightly)")
chor_cases <- c("C", "?")
meta18$quant_pheno <- ifelse(meta18$Phenotype %in% chor_cases, 2, ifelse(meta18$Phenotype %in% yolked_cases, 1, 0))

# make a note of any wierd, possibly reject samples in 18
meta18$reject <- 0
meta18$reject[meta18$Phenotype == "N"] <- 1
meta18$reject[meta18$ID == "FIT?"] <- 1

# add year notes
meta16$year <- "2016"
meta18$year <- "2018"

# pull out barcode and plate from bamlist
bamtable <- data.frame(Plate = substr(bamlist[,1], 6, 12),
                       Index = substr(bamlist[,1], 32, 39),
                       ord = 1:nrow(bamlist))
bamtable$Plate <- gsub("_", "", bamtable$Plate)
barcodes$Plate <- tolower(barcodes$Plate)

# merge
sample.meta <- merge(bamtable, barcodes, by = c("Plate", "Index"), sort = F, all.x = T)
sample.meta <- merge(sample.meta, meta18, by.x = "SampleID", by.y = "ID", sort = F, all.x = T)
sample.meta <- sample.meta[-which(duplicated(sample.meta$ord)),] # due to the fact that I had two vials with the same ID. They have different dates and unique metadata, so it's not a huge deal. Just screws up the merge.
sample.meta <- merge(sample.meta, meta16, by.x = "SampleID", by.y = "ID", sort = F, all.x = T)
sample.meta <- sample.meta[order(sample.meta$ord),]

# fix pheno columns
sample.meta$strict_pheno <- sample.meta$strict_pheno.x
sample.meta$strict_pheno[which(is.na(sample.meta$strict_pheno))] <- 
  sample.meta$strict_pheno.y[which(is.na(sample.meta$strict_pheno))]

sample.meta$relaxed_pheno <- sample.meta$relaxed_pheno.x
sample.meta$relaxed_pheno[which(is.na(sample.meta$relaxed_pheno))] <- 
  sample.meta$relaxed_pheno.y[which(is.na(sample.meta$relaxed_pheno))]

sample.meta$likely_pheno <- sample.meta$likely_pheno.x
sample.meta$likely_pheno[which(is.na(sample.meta$likely_pheno))] <- 
  sample.meta$likely_pheno.y[which(is.na(sample.meta$likely_pheno))]

sample.meta$quant_pheno <- sample.meta$quant_pheno.x
sample.meta$quant_pheno[which(is.na(sample.meta$quant_pheno))] <- 
  sample.meta$quant_pheno.y[which(is.na(sample.meta$quant_pheno))]

sample.meta$year <- sample.meta$year.x
sample.meta$year[which(is.na(sample.meta$year))] <- 
  sample.meta$year.y[which(is.na(sample.meta$year))]

# clean
sample.meta <- sample.meta[,c(1, 2, 3, 5, 29:33)]

# save
write.table(sample.meta, "data/sample_metadata.txt", sep = "\t", quote = F, col.names = T)


