source("C:/Users/user/Desktop/研究助理/simulation/function_10_29.R", echo = F)

# default
mypattern_default <- matrix(c(0,1,1,1,1,
                              1,0,1,1,1,
                              1,1,0,1,1,
                              1,1,1,0,1,
                              1,1,1,1,0), byrow = T, ncol = 5)
mypred <- matrix(c(0,1,1,1,1,1,
                   1,0,1,1,1,1,
                   1,1,0,1,1,1,
                   1,1,1,0,1,1,
                   1,1,1,1,0,1,
                   1,1,1,1,1,0), byrow = T, ncol = 6)
colnames(mypred) <- c(paste0("X", 1:5), "Y")
rownames(mypred) <- c(paste0("X", 1:5), "Y")

# P1
my_X_P1 <- c("rnorm(n=n_sample, 10,10)", "rnorm(n=n_sample, 0,100)",
          "rnorm(n=n_sample, 10,10)", "rnorm(n=n_sample, 0,100)",
          "rexp(n=n_sample, 1)")

P1 <- Simu_lambda_final(all_norm = F, n_count = 1000,
                        n_sample = 100, n_x = 5,
                        X_str = my_X_P1, important = c(1,5),
                        predictor_M = mypred, est_x = 1:5,
                        alpha = 0.05, pattern = mypattern_default,
                        variance = 1, pf = rep(1,5), adWeight = rep(1,5),
                        m = 5, meth = "pmm")

# P1_ans
P1_ans <- Simu_lambda_final(all_norm = F, n_count = 1000,
                        n_sample = 100, n_x = 5,
                        X_str = my_X_P1, important = c(1,5),
                        predictor_M = mypred, est_x = 1:5,
                        alpha = 0.05, pattern = mypattern_default,
                        variance = 1, pf = rep(1,5), adWeight = rep(1,5),
                        m = 5, meth = "norm")

# P1_ans_1
mypred_P1_ans_1 <- matrix(c(0,0,0,0,1,1,
                            1,0,0,0,1,1,
                            1,0,0,0,1,1,
                            1,0,0,0,1,1,
                            1,0,0,0,0,1,
                            1,1,1,1,1,0), byrow = T, ncol = 6)
colnames(mypred_P1_ans_1) <- c(paste0("X", 1:5), "Y")
rownames(mypred_P1_ans_1) <- c(paste0("X", 1:5), "Y")

P1_ans_1 <- Simu_lambda_final(all_norm = F, n_count = 1000,
                            n_sample = 100, n_x = 5,
                            X_str = my_X_P1, important = c(1,5),
                            predictor_M = mypred_P1_ans_1, est_x = 1:5,
                            alpha = 0.05, pattern = mypattern_default,
                            variance = 1, pf = rep(1,5),
                            adWeight = rep(1,5), m = 5, meth = "pmm")


# P2
my_X_P2 <- c("rnorm(n=n_sample, 10,10)", "rnorm(n=n_sample, 0,100)",
             "rt(n=n_sample, 1)", "rgamma(n=n_sample, shape = 2, scale = 3)",
             "rexp(n=n_sample, 1)")
P2 <- Simu_lambda_final(all_norm = F, n_count = 1000,
                        n_sample = 100, n_x = 5,
                        X_str = my_X_P2, important = c(1,5),
                        predictor_M = mypred, est_x = 1:5,
                        alpha = 0.05, pattern = mypattern_default,
                        variance = 1, pf = rep(1,5), adWeight = rep(1,5),
                        m = 5, meth = "pmm")

