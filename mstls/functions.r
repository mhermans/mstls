# Projection functions
# --------------------


library(ggplot2)
library(reshape)
library(stringr)
#library(expm) # needed for %^% operator

project <- function(init, pmat, nsteps, lbls=NULL) {
  # linear, time-invariant projection
  
  nclasses <- length(init)
  pops <- matrix(nrow=nsteps + 1, ncol=nclasses)

  colnames(pops) <- lbls
  rownames(pops) <- paste('t', seq(0, nsteps), sep='')
  
  pops[1,] <- init

  
  # iterate for intermediate values, else: A %^% nsteps %*% n0
  i <- nsteps
  n <- init
  while (i > 0) {
    n <- pmat %*% n
    pops[nsteps + 2 - i,] <- n
    i <- i - 1
  }
  
  pops
}


plot_proj <- function(proj_result){
  # projection results: plot absolute numbers and proportions
  
  p_num <- melt(proj_result)
  p_num$type <- rep('num', nrow(p_num))
  
  p_prop <- proj_result/rowSums(proj_result)
  p_prop <- melt(p_prop)
  p_prop$type <- rep('prop', nrow(p_prop))
  
  p <- rbind(p_num, p_prop)
  names(p) <- c('tlabel', 'class', 'value', 'type')
  p$time <- as.integer(str_replace(p$tlabel, 't', ''))
  
  q <- ggplot(p, aes(x=time, y=value, group=class, colour=class))
  q <- q + geom_line() + facet_grid(.~type)
  
  q
}


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


leslie_old <- function(F, P) {
  
  #construct (non-generalized) leslie matrix
  n_ageclass <- length(F)
  
  # treat last P-cat as open-ended? -> not included in Leslie?
  P <- head(P, -1) # drop last P
  
  L <- diag(P)
  L <- cbind(L, rep(0, n_ageclass-1))
  L <- rbind(F, L)
  
  L
}



leslie <- function(fertility, survivorship) {
  # Construction a generalized Leslie matrix
  
  B <- fertility
  S <- survivorship
  
  # S & B matrices do not contain last interval => n-1
  n_ages <- length(B) 
  
  G <- diag(rep(0, n_ages*2))
  
  j <- 1
  for (i in seq(1, n_ages*2, 2)) {
    G[i:(i+1), i:(i+1)] <- S[[j]]
    j <- j+1
  }
  
  G <- rbind(do.call(cbind, B), G) # add B-matrices on top
  
  # add double column of zero's at the end
  G <- cbind(G, matrix(rep(0, 36*2), ncol=2, nrow=36))
  
  G
  
}


# ########################################## #
# Replicatiefuncties voor Willekens & Rogers #
# ########################################## #


transfer_rates <- function(death_rates, migration_rates) {
  # construct M_x, the matrix of observed transfer rates, 
  # based on the observed death and outmigration rates.
  
  n_age <- nrow(death_rates)
  M <- list()
  
  OR <- death_rates + migration_rates
  
  for (x in 1:n_age) {
    # TODO: generalize from 2 regions/states to n states
    Mx <- diag(OR[x,])
    Mx[2,1] <- -MR[x,1]
    Mx[1,2] <- -MR[x,2]
    
    colnames(Mx) <- colnames(death_rates)
    
    M[[x]] <- Mx
    
  }
  
  M
  
}

transfer_prob <- function(transfer_rates) {
  # Calculate transfer probabilities Px,
  # based on observed transfer rates matrices Mx
  
  n_age <- length(transfer_rates)
  n_state <- ncol(transfer_rates[[1]])
  I <- diag(rep(1, n_state))
  
  P <- lapply(transfer_rates, function(Mx) 
  { solve( I + 2.5 * Mx) %*% ( I - 2.5 * Mx ) } )
  
  P[[n_age]] <- P[[n_age]] * 0 # all prob. are zero in final age category
  
  P
}


expected_survivors <- function(transition_probs, radix=100000) {
  # Calculate the expected number of survivors at exact age x,
  # based on the transition probabilities give radix/unit
  
  n_state <- ncol(transition_probs[[1]])
  
  l0 <- diag(rep(radix, n_state)) # default radix: 100k  
  L <- list(l0)
  
  for (x in 1:(length(transition_probs)-1)) {
    L[[x+1]] <- transition_probs[[x]] %*% L[[x]]
  }
  
  L
  
}


years_lived <- function(expected_survivors) {
  # Calculate number of years lived in each state by unit birth cohort
  # based on the matrices giving the expected number of survivors  
  
  n_ages <- length(L.surv)
  
  # linear method, cf. pg 70 (Schoen, 1988)
  
  L.surv <- expected_survivors
  L.dur <- list()
  
  for (x in 1:(n_ages-1)) {
    L.dur[[x]] <- 2.5 * (L.surv[[x]] + L.surv[[x+1]]) %*% solve(L.surv[[1]])
  }

  
  # Last half-open interval: L_z = M_z^{-1} l_z l*0^{-1}
  L.dur[[n_ages]] <- solve(M[[n_ages]]) %*% L.surv[[n_ages]] %*% solve(L.surv[[1]])
    
  L.dur
    
}

survivor_prop <- function(years_lived) {
  # Calculate the proportion of survivors based on the number
  # of years lived.
  
  # TODO: formule 3.6 pg 58: van M naar S
  
  S <- list()
  
  for (x in 1:(length(years_lived)-1)) {
    S[[x]] <- years_lived[[x+1]] %*% solve(years_lived[[x]])
  }
  
  S
  
}


birth_prop <- function(birth_rates, transition_probs, survivor_prop) {
  # Calculate proportion of births based on the observed fertility rate,
  # the transition probabilities and the proportion of survivors
  
  # based on: B_x = \fraq{5}{4} ( P_0 + I ) ( F{x+5} S_x )
  # TODO: equivalent voor multiple transitions?
  
  n_age <- nrow(birth_rates)
  n_state <- ncol(birth_rates)
  I <- diag(rep(1, n_state))
  
  F <- birth_rates
  P <- transition_probs
  S <- survivor_prop
  B <- list()
  
  for (x in 1:(n_age - 1)) {
    Fx <- diag(unlist(F[x,1:n_state]))
    Fx5 <- diag(unlist(F[x+1,1:n_state]))
    
    B[[x]] <- 5/4 * (P[[1]] + I) %*% ( Fx + Fx5 %*% S[[x]] )
    
  }
  
  B
  
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