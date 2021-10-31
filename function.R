library(mice)
library(miselect)

# all normal
Simu_lambda_norm <- function(n_count, n_sample, mu, sigma, n_x, important,
                             predictor_M, est_x, alpha, pattern, variance,
                             pf, adWeight, m, meth){
  # value
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
  galassodf <- list()
  for (i in 1:100){
    galassodf[[i]] <- data.frame()
  }
  # important
  z <- NULL
  d <- NULL
  if(length(important) > 1){
    for(i in 1:length(important)-1){
      d[i] <- paste0("X[,", important[i], "]+")
    }
    d[length(important)] <- paste0("X[,", important[length(important)], "]")
    for(i in 1:length(important)){
      z <- paste0(z, d[i])
    }
  }else
    z <- paste0("X[,", important, "]")
  # est_x
  a <- NULL
  d <- NULL
  if(length(est_x) > 1){
    for(i in 1:length(est_x)-1){
      d[i] <- paste0("X", est_x[i], "+")
    }
    d[length(est_x)] <- paste0("X", est_x[length(est_x)])
    for(i in 1:length(est_x)){
      a <- paste0(a, d[i])
    }
  }else
    a <- paste0("X", est_x)
  
  start <- Sys.time()
  for(i in 1:n_count){
    # create original data
    X <- MASS::mvrnorm(n = n_sample, mu = mu, Sigma = sigma)
    colnames(X) <- paste0("X", 1:n_x)
    error <- rnorm(n=n_sample,0,variance)
    Y <- as.data.frame(eval(parse(text = z))+error)
    colnames(Y) <- "Y"
    df <- cbind(X, Y)

    # create amputed data
    amp <- ampute(X, patterns = pattern)
    ampdf <- cbind(amp$amp, Y)
    
    # execute imputation
    imp <- mice(ampdf, print = FALSE,
                predictorMatrix = predictor_M, m = m,
                method = meth)
    
    # create the df of the pooled result
    expr <- paste0("lm(Y~",a, ")")
    result <- with(imp, eval(parse(text = expr)))
    pool <- pool(result)
    estimatedf <- rbind(estimatedf, pool$pooled$estimate)
    ubardf <- rbind(ubardf, pool$pooled$ubar)
    bdf <- rbind(bdf, pool$pooled$b)
    tdf <- rbind(tdf, pool$pooled$t)
    
    # create galasso
    galasso_dfs <- lapply(1:5, function(i) complete(imp, action = i))
    galasso_X <- list()
    galasso_Y <- list()
    for (i in 1:5) {
      galasso_X[[i]] <- as.matrix(galasso_dfs[[i]][, paste0("X", 1:n_x)])
      galasso_Y[[i]] <- galasso_dfs[[i]]$Y
    }
    galasso_fit <- galasso(galasso_X, galasso_Y, pf, adWeight)
    for (i in 1:length(galasso_fit$lambda)){
      galasso_coef <- coef(galasso_fit, lambda = galasso_fit$lambda[i])
      galassodf[[i]] <- rbind(galassodf[[i]], galasso_coef)
    }
    
    # create the pooled estimate
    estimate <- pool$pooled$estimate + estimate
    ubar <- pool$pooled$ubar + ubar
    b <- pool$pooled$b + b
    t <- pool$pooled$t + t
    
    # test regression coefficient
    p.valuedf <- rbind(p.valuedf, summary(pool)$p.value)
    
    # lambda
    lambdadf <- rbind(lambdadf, pool$pooled$lambda)
    lambda <- pool$pooled$lambda +lambda
  }
  end <- Sys.time()
  
  # example
  e_reg <- lm(Y~., data = df)
  e_ampreg <- lm(Y~., data = ampdf)
  
  var_name <- c("intercept", paste0("X", 1:n_x))
  colnames(estimatedf) <- var_name
  colnames(lambdadf) <- var_name
  colnames(ubardf) <- var_name
  colnames(bdf) <- var_name
  colnames(tdf) <- var_name
  colnames(p.valuedf) <- var_name
  for(i in 1:100){
    colnames(galassodf[[i]]) <- var_name
  }
  names(estimate) <- var_name
  names(ubar) <- var_name
  names(b) <- var_name
  names(t) <- var_name
  names(lambda) <- var_name
  m_est <- estimate/n_count
  m_ubar <- ubar/n_count
  m_b <- b/n_count
  m_t <- t/n_count
  m_lambda <- lambda/n_count
  p.value <- apply(p.valuedf, 2, function(x){x<alpha})
  p.value <- apply(p.value, 2, sum)
  runtime <- end-start           
  true_model <- paste0("Y = ", z, "+error")
  est_model <- paste0("lm(Y~",a, ")")
  
  mean_list <- list("mean_estimate" = m_est,
                      "mean_ubar" = m_ubar,
                      "mean_b" = m_b,
                      "mean_t" = m_t,
                      "mean_lambda" = m_lambda,
                      "p_value" = p.value)
  df_list <- list("estimatedf" = estimatedf, "lambdadf" = lambdadf, 
                  "ubardf" = ubardf, "bdf" = bdf, "tdf" = tdf,
                  "p.valuedf" = p.valuedf, "galassodf" = galassodf)
  all_equal_list <- list("est" = all.equal(as.numeric(estimate),as.numeric(apply(estimatedf, 2, sum))), 
                         "lambda" = all.equal(as.numeric(lambda), as.numeric(apply(lambdadf, 2, sum))),
                         "ubar" = all.equal(as.numeric(ubar), as.numeric(apply(ubardf, 2, sum))),
                         "b" = all.equal(as.numeric(b), as.numeric(apply(bdf, 2, sum))),
                         "t" = all.equal(as.numeric(t), as.numeric(apply(tdf, 2, sum))))
  example_list <- list("e_df" = df, "e_error" = error, "e_ampdf" = ampdf, "e_reg_summary" = summary(e_reg),
                       "e_ampreg_summary" = summary(e_ampreg), "e_pool_summary" = summary(pool))
  result_list <- list("model"=true_model, "mean"=mean_list, "df"=df_list,
                      "all_equal"=all_equal_list, "runtime"=runtime, "example" = example_list,
                      "missing_pattern" = md.pattern(ampdf, plot = F))
  
  return(result_list)
}

# not all normal
Simu_lambda1_other <- function(n_count, n_sample, n_x, X_str,
                         important, predictor_M, est_x, alpha, pattern, variance,
                         pf, adWeight, m, meth){
  # value
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
  galassodf <- list()
  for (i in 1:100){
    galassodf[[i]] <- data.frame()
  }
  # important
  z <- NULL
  d <- NULL
  if(length(important) > 1){
    for(i in 1:length(important)-1){
      d[i] <- paste0("X[,", important[i], "]+")
    }
    d[length(important)] <- paste0("X[,", important[length(important)], "]")
    for(i in 1:length(important)){
      z <- paste0(z, d[i])
    }
  }else
    z <- paste0("X[,", important, "]")
  # impute_x
  a <- NULL
  d <- NULL
  if(length(est_x) > 1){
    for(i in 1:length(est_x)-1){
      d[i] <- paste0("X", est_x[i], "+")
    }
    d[length(est_x)] <- paste0("X", est_x[length(est_x)])
    for(i in 1:length(est_x)){
      a <- paste0(a, d[i])
    }
  }else
    a <- paste0("X", est_x)
  
  start <- Sys.time()
  for(i in 1:n_count){
    # create original data
    X <- data.frame()
    for(j in 1:length(X_str)){
      X <- rbind(X, eval(parse(text = X_str[j])))
    }
    X <- t(X)
    colnames(X) <- paste0("X", 1:n_x)
    error <- rnorm(n=n_sample,0,variance)
    Y <- as.data.frame(eval(parse(text = z))+error)
    colnames(Y) <- "Y"
    df <- cbind(X, Y)

    # create amputed data
    amp <- ampute(X, patterns = pattern)
    ampdf <- cbind(amp$amp, Y)
    
    # execute imputation
    imp <- mice(ampdf, print = FALSE,
                predictorMatrix = predictor_M, m = m,
                method = meth)
    
    # estimate regression coefficient
    expr <- paste0("lm(Y~",a, ")")
    result <- with(imp, eval(parse(text = expr)))
    pool <- pool(result)
    estimatedf <- rbind(estimatedf, pool$pooled$estimate)
    ubardf <- rbind(ubardf, pool$pooled$ubar)
    bdf <- rbind(bdf, pool$pooled$b)
    tdf <- rbind(tdf, pool$pooled$t)
    
    # create galasso
    galasso_dfs <- lapply(1:5, function(i) complete(imp, action = i))
    galasso_X <- list()
    galasso_Y <- list()
    for (i in 1:5) {
      galasso_X[[i]] <- as.matrix(galasso_dfs[[i]][, paste0("X", 1:n_x)])
      galasso_Y[[i]] <- galasso_dfs[[i]]$Y
    }
    galasso_fit <- galasso(galasso_X, galasso_Y, pf, adWeight)
    for (i in 1:length(galasso_fit$lambda)){
      galasso_coef <- coef(galasso_fit, lambda = galasso_fit$lambda[i])
      galassodf[[i]] <- rbind(galassodf[[i]], galasso_coef)
    }
    
    # create the pooled estimate
    estimate <- pool$pooled$estimate + estimate
    ubar <- pool$pooled$ubar + ubar
    b <- pool$pooled$b + b
    t <- pool$pooled$t + t
    
    # test regression coefficient
    p.valuedf <- rbind(p.valuedf, summary(pool)$p.value)
    
    # lambda
    lambdadf <- rbind(lambdadf, pool$pooled$lambda)
    lambda <- pool$pooled$lambda +lambda
  }
  end <- Sys.time()
  
  # example
  e_reg <- lm(Y~., data = df)
  e_ampreg <- lm(Y~., data = ampdf)
  
  var_name <- c("intercept", paste0("X", 1:n_x))
  colnames(estimatedf) <- var_name
  colnames(lambdadf) <- var_name
  colnames(ubardf) <- var_name
  colnames(bdf) <- var_name
  colnames(tdf) <- var_name
  colnames(p.valuedf) <- var_name
  for(i in 1:100){
    colnames(galassodf[[i]]) <- var_name
  }
  names(estimate) <- var_name
  names(ubar) <- var_name
  names(b) <- var_name
  names(t) <- var_name
  names(lambda) <- var_name
  m_est <- estimate/n_count
  m_ubar <- ubar/n_count
  m_b <- b/n_count
  m_t <- t/n_count
  m_lambda <- lambda/n_count
  p.value <- apply(p.valuedf, 2, function(x){x<alpha})
  p.value <- apply(p.value, 2, sum)
  runtime <- end-start           
  true_model <- paste0("Y = ", z, "+error")
  est_model <- paste0("lm(Y~",a, ")")
  
  mean_list <- list("mean_estimate" = m_est,
                    "mean_ubar" = m_ubar,
                    "mean_b" = m_b,
                    "mean_t" = m_t,
                    "mean_lambda" = m_lambda,
                    "p_value" = p.value)
  df_list <- list("estimatedf" = estimatedf, "lambdadf" = lambdadf, 
                  "ubardf" = ubardf, "bdf" = bdf, "tdf" = tdf,
                  "p.valuedf" = p.valuedf, "galassodf" = galassodf)
  all_equal_list <- list("est" = all.equal(as.numeric(estimate),as.numeric(apply(estimatedf, 2, sum))), 
                         "lambda" = all.equal(as.numeric(lambda), as.numeric(apply(lambdadf, 2, sum))),
                         "ubar" = all.equal(as.numeric(ubar), as.numeric(apply(ubardf, 2, sum))),
                         "b" = all.equal(as.numeric(b), as.numeric(apply(bdf, 2, sum))),
                         "t" = all.equal(as.numeric(t), as.numeric(apply(tdf, 2, sum))))
  example_list <- list("e_df" = df, "e_error" = error, "e_ampdf" = ampdf, "e_imp" = imp,
                       "e_reg_summary" = summary(e_reg),
                       "e_ampreg_summary" = summary(e_ampreg), "e_pool_summary" = summary(pool))
  result_list <- list("model"=true_model, "mean"=mean_list, "df"=df_list,
                      "all_equal"=all_equal_list, "runtime"=runtime, "example" = example_list,
                      "missing_pattern" = md.pattern(ampdf, plot = F))
  
  return(result_list)
}

Simu_lambda_final <- function(all_normal, n_count, n_sample, mu, sigma, n_x,
                              X_str, important, predictor_M, est_x, alpha,
                              pattern, variance, pf, adWeight, m, meth){
  if(all_normal == TRUE)
    result <- Simu_lambda_norm(n_count, n_sample, mu, sigma, n_x, important,
                               predictor_M, est_x, alpha, pattern, variance,
                               pf, adWeight, m, meth)
  else if(all_normal == FALSE)
    result <- Simu_lambda1_other(n_count, n_sample, n_x, X_str,
                                 important, predictor_M, est_x, alpha, pattern,
                                 variance, pf, adWeight, m, meth)
  else
    print("Input legitimate value of all_normal")
  return(result)
}

galasso_mean <- function(galassodf, lamb){
  apply(galassodf[[lamb]], 2, mean)
}