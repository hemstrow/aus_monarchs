likely <- read.table("results/association_out_likely.lrt0", header = T, stringsAsFactors = F)
relaxed <- read.table("results/association_out_relaxed.lrt0", header = T, stringsAsFactors = F)
strict <- read.table("results/association_out_strict.lrt0", header = T, stringsAsFactors = F)
quant <- read.table("results/association_out_quant.lrt0", header = T, stringsAsFactors = F)

#colnames(combined)[6] <- "LRT_likely"
combined <- merge(likely, relaxed, by = colnames(likely)[1:5], sort = F)
combined <- merge(combined, strict, by = colnames(combined)[1:5], sort = F)
combined <- merge(combined, quant, by = colnames(combined)[1:5], sort = F)
colnames(combined)[6:9] <- paste0("LRT_", c("likely", "relaxed", "strict", "quant"))


scomb <- combined[order(combined$LRT_likely, decreasing = T),]

library(qqman)

mcomb <- combined
mcomb$Chromosome <- gsub("DPSCF3", "", mcomb$Chromosome)
mcomb$Chromosome <- as.numeric(mcomb$Chromosome)
mcomb$Position <- as.numeric(mcomb$Position)
mcomb$P_likely <- pchisq(mcomb$LRT_likely, 1, lower.tail = F)
mcomb$P_relaxed <- pchisq(mcomb$LRT_relaxed, 1, lower.tail = F)
mcomb$P_strict <- pchisq(mcomb$LRT_strict, 1, lower.tail = F)
mcomb$P_quant <- pchisq(mcomb$LRT_quant, 1, lower.tail = F)
colnames(mcomb)[1:2] <- c("CHR", "BP")
plikely <- qqman::manhattan(mcomb, p = "P_likely")
pstrict <- qqman::manhattan(mcomb, p = "P_strict")
prelaxed <- qqman::manhattan(mcomb, p = "P_relaxed")
pquant <- qqman::manhattan(mcomb, p = "P_quant")
