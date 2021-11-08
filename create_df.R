# create raw dataframes
# 
create_df <- function(n_count, n_sample, n_x, X_dist=NA, mu=rep(0,n_x),
                      sigma=diag(n_x), important, intercept=0, variance=1){
  # e.g. X[,1]+X[,2]
  important_expr <- NULL
  temp1 <- NULL
  if(length(important) > 1){
    for(i in 1:length(important)-1){
      temp1[i] <- paste0("X[,", important[i], "]+")
    }
    temp1[length(important)] <- paste0("X[,", important[length(important)], "]")
    for(i in 1:length(important)){
      important_expr <- paste0(important_expr, temp1[i])
    }
  }else
    important_expr <- paste0("X[,", important, "]")
  
  # e.g. X1+X2
  lm_x <- 1:n_x
  lm_expr <- NULL
  temp2 <- NULL
  if(length(lm_x) > 1){
    for(i in 1:length(lm_x)-1){
      temp2[i] <- paste0("X", lm_x[i], "+")
    }
    temp2[length(lm_x)] <- paste0("X", lm_x[length(lm_x)])
    for(i in 1:length(lm_x)){
      lm_expr <- paste0(lm_expr, temp2[i])
    }
  }else
    lm_expr <- paste0("X", lm_x)
  
  # if n_count = 5, then create 5 raw dataframes
  Xs <- list()
  dfs <- list()
  errors <- list()
  for(i in 1:n_count){
    if (is.na(X_dist[1] == T)){
      X <- MASS::mvrnorm(n = n_sample, mu = mu, Sigma = sigma)
      colnames(X) <- paste0("X", 1:n_x)
      error <- rnorm(n=n_sample,intercept,variance)
      Y <- as.data.frame(eval(parse(text = important_expr))+error)
      colnames(Y) <- "Y"
      df <- cbind(X, Y)
    }else{
      X <- data.frame()
      for(j in 1:length(X_dist)){
        X <- rbind(X, eval(parse(text = X_dist[j])))
      }
      X <- t(X)
      colnames(X) <- paste0("X", 1:n_x)
      error <- rnorm(n=n_sample,intercept,variance)
      Y <- as.data.frame(eval(parse(text = important_expr))+error)
      colnames(Y) <- "Y"
      df <- cbind(X, Y)
    }
    Xs[[i]] <- X
    dfs[[i]] <- df 
    errors[[i]] <- error
  }
  result <- list("Xs"=Xs, "dfs"=dfs, "errors"=errors,
                 "model"=paste0("Y=",intercept,"+",important_expr,"+error"),
                 "lm"=lm_expr)
  return(result)
}