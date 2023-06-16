library(readr)
#===============pull in metadata and associate with barcodes
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
                       bamname = bamlist[,1],
                       ord = 1:nrow(bamlist))
bamtable$Plate <- gsub("_", "", bamtable$Plate)
barcodes$Plate <- tolower(barcodes$Plate)

# merge
meta <- merge(bamtable, barcodes, by = c("Plate", "Index"), sort = F, all.x = T)
meta <- merge(meta, meta18, by.x = "SampleID", by.y = "ID", sort = F, all.x = T)
meta <- meta[-which(duplicated(meta$ord)),] # due to the fact that I had two vials with the same ID. They have different dates and unique metadata, so it's not a huge deal. Just screws up the merge.
meta <- merge(meta, meta16, by.x = "SampleID", by.y = "ID", sort = F, all.x = T)
meta <- meta[order(meta$ord),]

# fix pheno columns
meta$strict_pheno <- meta$strict_pheno.x
meta$strict_pheno[which(is.na(meta$strict_pheno))] <- 
  meta$strict_pheno.y[which(is.na(meta$strict_pheno))]

meta$relaxed_pheno <- meta$relaxed_pheno.x
meta$relaxed_pheno[which(is.na(meta$relaxed_pheno))] <- 
  meta$relaxed_pheno.y[which(is.na(meta$relaxed_pheno))]

meta$likely_pheno <- meta$likely_pheno.x
meta$likely_pheno[which(is.na(meta$likely_pheno))] <- 
  meta$likely_pheno.y[which(is.na(meta$likely_pheno))]

meta$quant_pheno <- meta$quant_pheno.x
meta$quant_pheno[which(is.na(meta$quant_pheno))] <- 
  meta$quant_pheno.y[which(is.na(meta$quant_pheno))]

meta$year <- meta$year.x
meta$year[which(is.na(meta$year))] <- 
  meta$year.y[which(is.na(meta$year))]


# clean
colnames(meta)[7] <- "emm_date"
meta <- meta[-which(is.na(meta$year)),]

#=================biosample============
meta$sample_name <- paste0(meta$SampleID, "_", meta$Plate, "_", meta$Index)
organism <- "Danaus plexippus"
isolation_source <- "Raised from eggs produced by females collected on milkweed patches at Pinjara Hills, Queensland, Australia, -27.540750, 152.906306"
geo_loc_name <- "Australia: Queensland: -27.540750, 152.906306"
tissue <- ifelse(meta$year == 2018, "leg", "wing_clip")
collection_date <- ifelse(meta$year == 2018, "July-2018", "June-2019")
isolate <- meta$SampleID
isolate[which(meta$sample_name == "9_2_3_plate2_CGACTGGA")] <- "9_2_3.2"
sex <- "Female"
collected_by <- ifelse(meta$year == 2018, "William Hemstrom", "Micah Freedman")
description <- paste0("Development: ", ifelse(meta$likely_pheno == "M", "immature", "mature"),
                      ". Raised from eggs from female collected collected on milkweed patches at Pinjara Hills, Queensland. Individuals raised under declining photoperiod and phenotyped for reproductive development. Photo ID: ",
                      meta$Photo, ", see Dryad.\n")


biosample <- data.frame(sample_name = meta$sample_name,
                        organism = organism,
                        isolate = isolate,
                        isolation_source = isolation_source,
                        geo_loc_name = geo_loc_name,
                        tissue = tissue,
                        collection_date = collection_date,
                        collected_by = collected_by,
                        description = description)
write.table(biosample, "data/BioSample_meta.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

#================SRA===============
files <- read.table("data/fastqlist.txt")
patterns <- gsub("_R.+_", "_", files$V1)
patterns <- gsub(".fastq", "", patterns)
SRA <- data.frame(library_ID = gsub("rad_split_out", "", patterns))
SRA$library_ID <- gsub("_GG", "_", SRA$library_ID)
SRA$library_ID <- gsub("TGCAG$", "", SRA$library_ID)
SRA$Plate <- stringi::stri_extract_first_regex(SRA$library_ID, "plate_[1-3]")
SRA$Plate <- gsub("_", "", SRA$Plate)
SRA$Index <- gsub(".+_", "", SRA$library_ID)
SRA <- merge(SRA, meta[,c("Plate", "Index", "sample_name")], by = c("Plate", "Index"), all = TRUE)
SRA$Index <- NULL
SRA$Plate <- NULL
SRA$filename <- files$V1
SRA$filename2 <- files$V2
SRA <- na.omit(SRA)
SRA$title <- "RAD-Seq of Danaus plexippus: Legs and wings, samples from Queensland, Australia"
SRA$library_strategy <- "RAD-Seq"
SRA$library_source <- "GENOMIC"
SRA$library_selection <- "Reduced Representation"
SRA$library_layout <- "paired"
SRA$platform <- "ILLUMINA"
SRA$instrument_model <- "Illumina HiSeq 4000"
SRA$design_description <- "RAD-Seq according to Ali et al. (2016)"
SRA$filetype <- "fastq"


SRA <- SRA[,c(1:2,  5:ncol(SRA), 3:4)]
write.table(SRA, "data/SRA_metadata.txt", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
writeLines(c(SRA$filename, SRA$filename2), "data/files_for_ncbi.txt")
