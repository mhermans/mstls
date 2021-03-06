#' Two-region observed population characteristics for Yugoslavian females, 1961.
#' 
#' A dataset containing observed 1961 female population, births, deaths and migration between two regions, Slovenia
#' and the rest of Yugoslavia. Original study is, analyzed as demonstration dataset in Rogers (1975) 
#' and Willekens & Rogers (1978). 
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 36 rows and 7 variables
#' @name YU61female
NULL


collapse_interval <- function(df, interval=5) {
  # collapse a DF with 1 row per age to intervals
  
  max_age <- nrow(df)
  ages <- seq(0, max_age, interval)
  
  l <- list()
  for (i in ages) {
    l <- cbind(l, colSums(df[i:(i+(interval-1)),]))
  }
  
  df_i <- data.frame(matrix(unlist(l), nrow=length(ages), byrow=T))
  rownames(df_i) <- paste(ages, ages + (interval-1), sep='-')
  colnames(df_i) <- colnames(df)
  
  df_i
}


state_table <- function(data, nstate, rlbl, clbl) { 
  # Extract single state from multistate nested matrix
  
  tab <- do.call(rbind, lapply(data, function(x) x[,nstate])) 
  
  if ( missing(rlbl) ) {
    rlbl <- 1:length(data)
  }
  
  if ( missing(clbl) ) {
    clbl <- colnames(data[[1]])
  }
  
  # TODO: add col, rownames parameter
  rownames(tab) <- rlbl
  colnames(tab) <- clbl
  
  tab
  
}

stablepop_pct <- function(init_pop, n_states, pmat, n_steps=500) {
  # calculate pct. distr. for approx. stable equvalent population
  
  tn <- project(init_pop=init_pop, n_states=n_states, pmat=pmat, n_steps=n_steps)
  tn_1 <- t(tn[,seq(1,ncol(tn),2)])
  tn_2 <- t(tn[,seq(2,ncol(tn),2)])
  n1 <- tn_1[,ncol(tn_1)]
  n2 <- tn_2[,ncol(tn_2)]
  
  spop_pct <- cbind(n1/sum(n1),n2/sum(n2))
  
  spop_pct
}
