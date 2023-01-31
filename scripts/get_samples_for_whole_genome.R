library(data.table)
# colony results
dat <- readRDS("colony/colony_R_outputs.RDS")
colony <- dat$col_res
# data
dat <- dat$data
sm <- snpR::sample.meta(dat)
sm$colID <- colnames(dat)
# merge
sm <- merge(sm, colony$clusters[, c(1,2,3)], by.x = "colID", by.y = "OffspringID")
sm <- sm[order(sm$ClusterIndex),]


# pick samples to minimize family representation
n <- 92

migs <- as.data.table(sm[sm$likely_pheno == "M",])
res <- as.data.table(sm[sm$likely_pheno == "R",])

## function to pick evenly
pick_evenly <- function(x, n, by){
  x <- x[sample(nrow(x), nrow(x), F),]
  idstring <- stringi::stri_rand_strings(1, 5)
  x[,eval(idstring) := 1:nrow(x)]
  
  # grab
  keep <- as.data.table(matrix(NA, 0, ncol(x)))
  colnames(keep) <- colnames(x)
  while(nrow(keep) <= n){
    keep <- rbind(keep, x[,head(.SD, 1), by = by])
    
    x <- x[-which(x[,..idstring][[1]] %in% keep[,..idstring][[1]]),]
  }
  
  # trim
  keep <- head(keep, n)
  keep <- keep[,-..idstring]
  return(keep)
}

mig_keep <- pick_evenly(migs, n/2, "ClusterIndex")
res_keep <- pick_evenly(res, n/2, "ClusterIndex")
keeps <- rbind(mig_keep, res_keep)
keeps <- dplyr::arrange(keeps, Plate, Well)
keeps$PlateID <- NA
keeps$PlateID[which(keeps$Plate == "plate1")] <- "DNAA_530"
keeps$PlateID[which(keeps$Plate == "plate2")] <- "DNAA_535"
keeps$PlateID[which(keeps$Plate == "plate3")] <- "DNAA_536"

keep.cols <- c("PlateID", "Well", "SampleID", "ClusterIndex")
ret <- keeps[,..keep.cols]
write.table(ret, "data/samples_for_whole_genome.txt", quote = F, col.names = T, row.names = F, sep = "\t")
