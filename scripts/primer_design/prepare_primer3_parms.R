args <- commandArgs(T)
args <- as.character(args)
parms <- readLines(args[1])
regions <- readLines(args[2])
fasta <- readLines(args[3])

parms <- c(paste0("SEQUENCE_TEMPLATE=", fasta), parms)
out <- vector("list", length(regions))

for(i in 1:length(regions)){
  out[[i]] <- c(paste0("SEQUENCE_ID=", args[4], "_", i), paste0("SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=", regions[i]), parms)
}

writeLines(unlist(out), paste0("parms_", args[4]), sep = "\n")
