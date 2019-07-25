# function to run:
refine_rf <- function(rf, response, num.trees, trim_cuttoffs = NULL, trim = 0.5, seach_cuttoff = NULL, par = FALSE, ...){
  #=============subfunctions:=========
  # trim snps below a given importance quantile
  remove_trim <- function(imps, trim){
    best.imp <- quantile(imps, trim)
    best.imp <- which(imps >= best.imp[1])
    return(best.imp)
  }
  # trim to a specific number of snps
  remove_set_to_n_snps <- function(imps, n_snps){
    imps <- data.frame(imps, ord = 1:length(imps))
    imps <- dplyr::arrange(imps, desc(imps))
    imps <- imps[1:n_snps,]
    imps <- dplyr::arrange(imps, ord)
    best.imp <- imps$ord
    return(best.imp)
  }
  # remove a single snp
  remove_single <- function(imps){
    best.imp <- which.min(imps)
    best.imp <- (1:nrow(rf$data@stats))[-best.imp]
    return(best.imp)
  }
  # set to the next level, then trim a given percentage
  remove_set_and_trim <- function(imps, n_snps, trim){
    # set
    imps <- data.frame(imps, ord = 1:length(imps))
    imps <- dplyr::arrange(imps, desc(imps))
    imps <- imps[1:n_snps,]
    imps <- dplyr::arrange(imps, ord)
    
    # trim
    best.imp <- quantile(imps[,1], trim)
    best.imp <- imps[which(imps[,1] >= best.imp[1]),]
    
    # return
    best.imp$ord
  }
  
  # remove according to cuttoffs and trim levels
  remove_snps <- function(imps, trim_cuttoffs, trim){

    # figure out which "bin" we are in
    n_snps <- length(imps)
    bin.logi <- which(n_snps > trim_cuttoffs)
    if(length(bin.logi) > 0){
      
      # grab the current bin and trim
      bin <- min(bin.logi)
      best.imp <- remove_trim(imps, trim[bin])
      
      # check if we changed bins!
      bin.t.trim <- which(length(best.imp) > trim_cuttoffs)
      if(!identical(bin.t.trim, bin.logi)){
        
        # if we haven't hit the last bin
        if(length(bin.t.trim) != 0){
          # we either want to set this to - the current trim or to the cuttoff - the next trim, whichever has more snps.
          
          best.imp.current.trim <- remove_trim(imps, trim[bin])
          best.imp.set.plus.next.trim <- remove_set_and_trim(imps, trim_cuttoffs[min(bin.t.trim) - 1], trim[min(bin.t.trim)])
          if(length(best.imp.current.trim) > length(best.imp.set.plus.next.trim)){
            best.imp <- best.imp.current.trim
          }
          else{
            best.imp <- best.imp.set.plus.next.trim
          }
        }
        
        
        ## if we've hit the last bin and are using a single cutoff, set to single cuttoff
        else if(length(bin.t.trim) == 0 & trim_single){
          best.imp <- remove_set_to_n_snps(imps, trim_cuttoffs[length(trim_cuttoffs)])
          if(length(best.imp) == n_snps){
            best.imp <- remove_single(imps)
          }
        }

        ## if we've hit the last bin but aren't using a single cuttoff, keep trimming as usual
        else if(length(bin.t.trim) == 0){
          best.imp.current.trim <- remove_trim(imps, trim[bin])
          best.imp.set.plus.next.trim <- remove_set_and_trim(imps, trim_cuttoffs[bin],  trim[length(trim)])
          if(length(best.imp.current.trim) > length(best.imp.set.plus.next.trim)){
            best.imp <- best.imp.current.trim
          }
          else{
            best.imp <- best.imp.set.plus.next.trim
          }
          
        }
      }
    }
    # if we've hit the last bin...
    else{
      if(trim_single){
        best.imp <- remove_single(imps)
      }
      else{
        best.imp <- remove_trim(imps, trim[length(trim)])
      }
    }
    
    
    
    return(best.imp)
  }
  
  #==========initialize============
  
  
  
  # intialize output
  out <- matrix(NA, 1000, 2)
  colnames(out) <- c("n_snps", "prediction_error")
  out[1,] <- c(nrow(rf$data), rf$models[[1]]$model$prediction.error)
  
  # initialize output confusion array
  conf.out <- array(NA, dim = c(dim(rf$models[[1]]$model$confusion.matrix), 1000))
  conf.out[,,1] <- rf$models[[1]]$model$confusion.matrix
  
  
  # intialize difference
  diff <- 1
  
  # find the target column name
  tar.col <- paste0(response, "_RF_importance")
  
  # fix the single run_cuttoff if NULL
  if(is.null(trim_cuttoffs)){
    trim_cuttoffs <- 1
  }
  
  # if we are provided with less trim_cuttoffs than trim levels, we assume no single trimming
  if(length(trim_cuttoffs) < length(trim)){
    trim_single <- FALSE
  }
  # if the same, we trim single SNPs beneath the lowest cutoff
  else if (length(trim_cuttoffs) == length(trim)){
    trim_single <- TRUE
  }
  else{
    stop("Not enough trim levels provided for given trim cuttoffs.\n")
  }
  
  #==========while loop================
  i <- 2
  while(diff > 0 & nrow(rf$data) > 1){
    
    # if we somehow reach the end of the output storage, make it bigger!
    if(i == nrow(out)){
      out <- rbind(out, matrix(NA, 100000000, 2))
      conf.out <- c(conf.out, array(NA, dim = c(dim(rf$models[[1]]$model$confusion.matrix), 1000)))
      conf.out <- array(conf.out, c(dim(rf$models[[1]]$model$confusion.matrix), 2000))
    }
    
    imps <- abs(rf$data@stats[[tar.col]])
    
    # get the snps to keep
    best.imp <- remove_snps(imps, trim_cuttoffs, trim)
    
    # subset
    suppressWarnings(input <- subset_snpR_data(rf$data, snps = best.imp))
    
    # report some things
    out[i,1] <- nrow(input)
    cat("Refinement: ", i - 1, "\n\tStarting prediction error: ", out[i-1, 2],
        "\n\tNumber of remaining SNPs: ", out[i,1], "\n\tBeginning rf...\n")
    
    # run rf
    suppressWarnings(rf <- run_random_forest(x = input, response = response, num.trees = num.trees, mtry = nrow(input), par = par, pvals = FALSE, ...))
    
    # save outputs
    out[i,2] <- rf$models[[1]]$model$prediction.error
    conf.out[,,i] <- rf$models[[1]]$model$confusion.matrix
    
    # set diff
    diff <- out[i-1,2] - out[i,2]
    
    if(!is.null(search_cuttoff) & diff <= 0){
      if(out[i,1] > search_cuttoff){
        diff<- 1
        warning("Increase in prediction error, but search cuttoff number of SNPs not passed, continuing refinement.\n")
      }
    }
    
    i <- i + 1
  }
  
  #==============return==============
  # return final model and outputs
  empties <- which(is.na(out[,1]))
  out <- out[-empties,]
  conf.out <- conf.out[,,-empties]
  
  return(list(model = rf, error_delta = out, confusion_matrices = conf.out))
}
