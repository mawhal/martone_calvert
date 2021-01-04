mcmc.list2array <-
  function(x) {
    arr <- array(NA_real_,
                 dim = c(nrow(x[[1]]), length(x), ncol(x[[1]])))
    dimnames(arr)[[3]] <- colnames(x[[1]])
    for(i in seq_along(x)) arr[,i,] <- as.matrix(x[[i]])
    return(arr)
  }