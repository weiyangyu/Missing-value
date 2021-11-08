library(mice)
# ampute created dataframes
ampute_df <- function(n_count, Xs, dfs, prop=0.5, patterns=NULL,
                      freq=NULL, bycases=TRUE){
  ampdf <- list()
  missing_pattern <- list()
  for(i in 1:n_count){
    amp <- ampute(data=Xs[[i]], prop=prop, patterns=patterns,
                  freq=freq, bycases=bycases)
    Y <- dfs[[i]]$Y
    ampdf[[i]] <- cbind(amp$amp, Y)
    missing_pattern[[i]] <- md.pattern(ampdf[[i]], plot = F)
  }
  result <- list("ampdf"=ampdf, "missing_pattern"=missing_pattern)
  return(result)
}
