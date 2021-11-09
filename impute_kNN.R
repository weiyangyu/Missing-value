library(VIM)
# impute with kNN and analyze the result
impute_kNN <- function(n_count, n_x, ampdf, k,
                       dist_var=colnames(ampdf[[1]]),
                       weights=NULL, numFun=median, catFun=maxCat,
                       weightDist=FALSE,
                       alpha=0.05, lm_expr){
  # mean
  estimate <- rep(0,n_x+1)
  lambda <- rep(0,n_x+1)
  ubar <- rep(0,n_x+1)
  b <- rep(0,n_x+1)
  t <- rep(0,n_x+1)
  
  # df
  estimatedf <- data.frame()
  lambdadf <- data.frame()
  ubardf <- data.frame()
  bdf <- data.frame()
  tdf <- data.frame()
  p.valuedf <- data.frame()
  
  # list
  imp_n_count <- list()
  imp_k <- list()
  result_n_count <- list()
  result_k <- list()
  pool_list <- list()
  
  expr <- paste0("lm(Y~",lm_expr, ")")
  
  start <- Sys.time()
  for(i in 1:n_count){
    for(j in 1:length(k)){
      imp <- kNN(data=ampdf[[i]], k=k[j], dist_var=dist_var,
                 weights=weights, numFun=numFun, catFun=maxCat,
                 weightDist=weightDist,
                 imp_var=F)
      imp_k[[j]] <- imp
      result_k[[j]] <- with(imp_k[[j]], eval(parse(text = expr)))
    }
    imp_n_count[[i]] <- imp_k
    result_n_count[[i]] <- result_k
    pool_list[[i]] <- pool(result_n_count[[i]])
    
    
    estimatedf <- rbind(estimatedf, pool_list[[i]]$pooled$estimate)
    lambdadf <- rbind(lambdadf, pool_list[[i]]$pooled$lambda)
    ubardf <- rbind(ubardf, pool_list[[i]]$pooled$ubar)
    bdf <- rbind(bdf, pool_list[[i]]$pooled$b)
    tdf <- rbind(tdf, pool_list[[i]]$pooled$t)

    estimate <- pool_list[[i]]$pooled$estimate + estimate
    lambda <- pool_list[[i]]$pooled$lambda +lambda
    ubar <- pool_list[[i]]$pooled$ubar + ubar
    b <- pool_list[[i]]$pooled$b + b
    t <- pool_list[[i]]$pooled$t + t
    
    p.valuedf <- rbind(p.valuedf, summary(pool_list[[i]])$p.value)
  }
  end <- Sys.time()
  runtime <- end - start
  
  # add name 
  var_name <- c("intercept", paste0("X", 1:n_x))
  colnames(estimatedf) <- var_name
  colnames(lambdadf) <- var_name
  colnames(ubardf) <- var_name
  colnames(bdf) <- var_name
  colnames(tdf) <- var_name
  colnames(p.valuedf) <- var_name
  names(estimate) <- var_name
  names(ubar) <- var_name
  names(b) <- var_name
  names(t) <- var_name
  names(lambda) <- var_name
  
  # calsulate the mean of each value
  m_est <- estimate/n_count
  m_ubar <- ubar/n_count
  m_b <- b/n_count
  m_t <- t/n_count
  m_lambda <- lambda/n_count
  
  # calculate frequency where p.value is lower that alpha
  p.value <- apply(p.valuedf, 2, function(x){x<alpha})
  p.value <- apply(p.value, 2, sum)
  
  imp_list <- list("imputation" = imp_n_count, 
                   "analysis" = result_n_count,
                   "pool" = pool_list)
  
  mean_list <- list("mean_estimate"=m_est,
                    "mean_ubar" = m_ubar,
                    "mean_b" = m_b,
                    "mean_t" = m_t,
                    "mean_lambda" = m_lambda)
  df_list <- list("estimatedf"=estimatedf, 
                  "lambdadf" = lambdadf, 
                  "ubardf" = ubardf, "bdf" = bdf, "tdf" = tdf)
  result <- list("imp_output"=imp_list, "mean_output"=mean_list,
                 "df_output"=df_list,"p.value_freq"=p.value,
                 "runtime" = runtime)
  return(result)
}



