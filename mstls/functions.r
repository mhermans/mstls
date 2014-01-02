# Projection functions
# --------------------

library(ggplot2)
library(reshape)
library(stringr)
#library(expm) # needed for %^% operator

# Based on:
# - Willekens, F. & Rogers, A. (1978), Spatial Population Analyssi: Methods and Computer Programs. IASSA: Laxenburg, Austria.
# - Rogers, A. (1975), Introduction to Multiregional Mathematical Demography. John Wiley & Sons: New York.
# - Schoen, R. (1988), Modeling Multigroup Populations. Plenus Press: New York.

# ============================== #
# MULTISTATE LIFETABLE FUNCTIONS #
# ============================== #


transfer_matrix <- function(observed_rates, multiple=TRUE, absorbing_state="last") {
  # Construct Mx, the matrix of observed transfer rates based on the 
  # observed death and outmigration rates. 
  
  # number non-absorbing states
  n_states <- ncol(observed_rates) - 3 - 1
  
  n_ages <- nrow(observed_rates) / n_states
  ages <- unique(observed_rates[,1])
  
  # TODO: make abs.state column selectable (currently must have one)
  if (absorbing_state == 'last') { absorbing <- n_states + 1 }
  
  
  # Convert "long columns" input format to lists of matrices #
  # -------------------------------------------------------- #
  
  # absorbing/death state mdx
  Md.cols <- matrix(
    observed_rates[, 3 + n_states + 1],
    ncol=n_states) 
  
  # non-absorbing transition rates
  Mx.cols <- observed_rates[, 4:(3+n_states)]
  
  # re-arrange Md cols to matrices
  Md <- list()
  for (i in 1:n_ages) {
    Md[[i]] <- diag(Md.cols[i,])
  }
  
  # re-arrange Mx cols to matrices
  Mx <- list()
  i <- 1  
  for (x in ages) {
    Mx[[i]] <- t(Mx.cols[observed_rates[,1] == x, ])
    i <- i+1
  }
  
  
  # Combine Mdx, Mx into transition matrix M #
  # ---------------------------------------- #
  
  # Diagonal matrices with total outmigration rate per state
  Mx_it <- lapply(lapply(Mx, colSums), diag)
  
  M <- list()
  for (i in 1:n_ages) {
    
    # diagonal elements are sum of deathr + total outmigr. per state
    M[[i]] <- Md[[i]] + Mx_it[[i]]
    
    # off-diagonal elements are negative mx-elements
    M[[i]] <- M[[i]] + (Mx[[i]] * -1)
    
    #rownames(M[[i]]) <- c('M_i1', 'M_i2', 'M_i3')
    #colnames(M[[i]]) <- c('M_1j', 'M_3j', 'M_3j')
  }
  
  
  # Correct M for single transition scenario #
  # ---------------------------------------- #
  
  # if no multiple transitions per interval: transition matrix for last, open interval
  # is a diagonal matrix with the regional death rates (can't die & migrate).
  # (W&R, pg. 52-53; Rogers, pg. 65)
  
  if (!multiple) { M[[n_ages]] <- Md[[n_ages]] }
  
  
  M
  
}



transfer_prob <- function(transfer_rates, multiple=TRUE) {
  # ================================================================== #
  # Calculate age-specific probabilities of first transfer Px,         #
  # based on observed observed transfer (death and outmigration) rates #
  # structured in the Mx matrix (allows multiple transitions)          #
  # ================================================================== #
  
  n_ages <- length(transfer_rates)
  n_states <- ncol(transfer_rates[[1]])
  I <- diag(rep(1, n_states))
  
  
  # ---------------------------------------------- #
  # Calculate probability, only single transitions #
  # uses Schoen (1988), formula 4.57, p. 82        #
  # ---------------------------------------------- #
  if (!multiple) {
    
    # TODO: raise error if Mz is not diagonal?
    
    P <- lapply(transfer_rates, function(Mx) {
      Mx_t <- t(Mx)
      Mx_t.diag <- diag(diag(Mx_t)) # need diagonal matrix
      t(I - 5 * solve(I + 2.5 * Mx_t.diag) %*% Mx_t)
    })
    
  }
  
  # ----------------------------------------------------------- #
  # Calculate probability, allowing for multiple transitions    #
  # Uses "option 3", cf. Willekens & Rogers (1978), p. 49-51.   #
  # ----------------------------------------------------------- #
  if (multiple) {  
    
    P <- lapply(transfer_rates, function(Mx) {
      solve( I + 2.5 * Mx) %*% ( I - 2.5 * Mx ) 
    })
    
  }
  
  # all prob. are zero in final age category
  P[[n_ages]] <- P[[n_ages]] * 0 
  
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


years_lived <- function(expected_survivors, transfer_rates) {
  # Calculate number of years lived in each state by unit birth cohort
  # based on the matrices giving the expected number of survivors,
  # and the transfer rates for the open interval
  
  L.surv <- expected_survivors
  Mx <- transfer_rates
  n_ages <- length(L.surv)
  
  # linear method, cf. pg 70 (Schoen, 1988)
  
  L.surv <- expected_survivors
  L.dur <- list()
  
  for (x in 1:(n_ages-1)) {
    L.dur[[x]] <- 2.5 * (L.surv[[x]] + L.surv[[x+1]]) %*% solve(L.surv[[1]])
  }
  
  
  # Last half-open interval: L_z = M_z^{-1} l_z l*0^{-1}
  # (in the case of single transitions, the last M-matrix is a diagonal death rate matrix)
  L.dur[[n_ages]] <- solve(Mx[[n_ages]]) %*% L.surv[[n_ages]] %*% solve(L.surv[[1]])
  
  
  L.dur
  
}


survivor_prop <- function(years_lived) {
  # Calculate the proportion of survivors based on the number
  # of years lived.
  
  # TODO: formule 3.6 pg 58 W&R: rechstreeks van M naar S?
  # TODO: vgl formule 3.5 proportions vs. probabilities, pg. 57 W&R
  
  S <- list()
  
  for (x in 1:(length(years_lived)-1)) {
    S[[x]] <- years_lived[[x+1]] %*% solve(years_lived[[x]])
  }
  
  S
  
}


birth_prop <- function(birth_rates, transition_probs, survivor_prop) {
  # Calculate proportion of births based on the observed fertility rate,
  # the transition probabilities and the proportion of survivors 
  
  n_ages <- length(transition_probs)
  n_states <- ncol(transition_probs[[1]])
  I <- diag(rep(1, n_states))
  
  P <- transition_probs
  S <- survivor_prop
  
  # rearrange Fx-vector to columns 
  F <- matrix(birth_rates, ncol=n_states)
  
  B <- list()
  for (x in 1:(n_ages - 1)) {
    Fx <- diag(unlist(F[x,1:n_states]))
    Fx5 <- diag(unlist(F[x+1,1:n_states]))
    
    # W&R: B_x = \fraq{5}{4} ( P_0 + I ) ( F{x+5} S_x )
    B[[x]] <- 5/4 * (P[[1]] + I) %*% ( Fx + Fx5 %*% S[[x]] )
    
  }
  
  B
  
}


# ================= #
# UTILITY FUNCTIONS #
# ================= #


projection_matrix <- function(fertility, survivorship) {
  # Construction a generalized Leslie matrix/projection matrix G
  
  B <- fertility
  S <- survivorship
  
  # S & B matrices do not contain last interval => n-1
  n_ages <- length(B) + 1 
  n_states <- ncol(B[[1]])
  
  G <- diag(rep(0, n_ages*n_states))
  
  j <- 1
  i_start <- head(seq(1, (n_ages*n_states), n_states), -1)
  mwidth <- n_states - 1
  for (i in i_start) {
    G[i:(i+mwidth), i:(i+mwidth)] <- S[[j]]
    j <- j+1
  }
  
  # add 0-matrix for Bz at the end
  B[[(n_ages)]] <- diag(rep(0, n_states)) 
  
  # add B-matrices on top of G, and drop superfluous last rows
  G <- rbind(do.call(cbind, B), G)
  G <- G[1:(n_ages*n_states),]
  
  G
  
}


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

stablepop_pct <- function(init, pmat, nsteps=500) {
  # calculate pct. distr. for approx. stable equvalent population
  
  tn <- project(init=init, pmat=pmat, nsteps=nsteps)
  tn_1 <- t(tn[,seq(1,ncol(tn),2)])
  tn_2 <- t(tn[,seq(2,ncol(tn),2)])
  n1 <- tn_1[,ncol(tn_1)]
  n2 <- tn_2[,ncol(tn_2)]
  
  spop_pct <- cbind(n1/sum(n1),n2/sum(n2))
  
  spop_pct
}
