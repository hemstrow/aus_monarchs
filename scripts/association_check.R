library(snpR); library(qqman)
# from angsd--quant pheno
chr_key <- fread("data/chr_key.txt")
adat_quant <- data.table::fread("scripts/association/association_out_quant_cov.lrt0.gz", header = T)
adat_quant$chr <- chr_key$V1[match(adat_quant$Chromosome, chr_key$V3)]
adat_quant$pval <- pchisq(adat_quant$LRT, 1, lower.tail = FALSE)
adat_quant$lpval <- -log(adat_quant$pval)
adat_quant$fdr <- p.adjust(adat_quant$pval, "fdr")
adat_quant <- dplyr::arrange(adat_quant, pval)

q_quant <- plot_qq(adat_quant, "pval")
p_quant <- plot_manhattan(adat_quant, "fdr", chr = "Chromosome", bp = "Position", log.p = T, significant = 0.05, suggestive = 0.01)