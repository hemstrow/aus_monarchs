find.match <- function(x, y, sample.meta){
  out <- vector("list", length(x))
  out.best <- data.frame(sample = x, matches = character(length(x)),
                         well = character(length(x)),
                         plate = character(length(x)),
                         percentage = numeric(length(x)), stringsAsFactors = F)
  
  for(i in 1:length(x)){
    t.samp.id <- which(sample.meta$Plate == "plate3" & sample.meta$Well== x[i])
    if(length(t.samp.id) == 0){
      out.best$matches[i] <- "no_x"
      next()
    }
    t.samp <- y[,t.samp.id]
    miss <- t.samp != "NN"
    out[[i]]$hits <- numeric(ncol(y))
    names(out[[i]]$hits) <- sample.meta$SampleID
    
    for(j in (1:ncol(y))){
      if(j == t.samp.id){
        out[[i]]$hits[j] <- 0
        next()
      }
      c.samp <- y[,j]
      c.miss <- c.samp != "NN"
      u.miss <- which(c.miss & miss)
      
      out[[i]]$hits[j] <- sum(t.samp[u.miss] == c.samp[u.miss])/length(u.miss)
    }
    best <- which.max(out[[i]]$hits)
    out[[i]]$best <- out[[i]]$hits[best]
    names(out[[i]]$best) <- names(out[[i]]$hits)[best]
    
    out.best$matches[i] <- paste0(names(out[[i]]$best), collapse = ", ")
    out.best$percentage[i] <- out[[i]]$best
    out.best$well[i] <- sample.meta$Well[best]
    out.best$plate[i] <- sample.meta$Plate[best]
  }
  return(list(best = out.best, data = out))
}