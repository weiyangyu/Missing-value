library(mice)
# impute with MICE and analyze the result
impute_mice <- function(n_count, n_x, ampdf, predictor_matrix, m=5,
                        meth=NULL, lm_expr, alpha=0.05){
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
  
  # imputed list
  imp_list <- list()
  pool_list <- list()
  
  # expression
  expr <- paste0("lm(Y~",lm_expr, ")")
  
  start <- Sys.time()
  for(i in 1:n_count){
    imp <- mice(ampdf[[i]], print=FALSE,
                predictorMatrix=predictor_matrix, m=m,
                method=meth)
    imp_list[[i]] <- imp
    result <- with(imp, eval(parse(text = expr)))
    pool_list[[i]] <- pool(result)
    
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
  
  impute_list <- list("imputation" = imp_list,
                      "pool" = pool_list)
  
  mean_list <- list("mean_estimate" = m_est,
                    "mean_ubar" = m_ubar,
                    "mean_b" = m_b,
                    "mean_t" = m_t,
                    "mean_lambda" = m_lambda)
  df_list <- list("estimatedf" = estimatedf, "lambdadf" = lambdadf, 
                  "ubardf" = ubardf, "bdf" = bdf, "tdf" = tdf)
  result <- list("imp_output"=impute_list, "mean_output"=mean_list,
                 "df_output"=df_list, "p.value_freq"=p.value, 
                 "runtime"=runtime)
  return(result)
}