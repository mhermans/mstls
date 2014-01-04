
#' @title Construct Mx, the matrix of observed transfer rates
#' 
#' @description \code{transfer_matrix} returns a list with age-specific transfer-rate matrices Mx,
#' based on the observed death and outmigration rates per state.
#' 
#' @param observed_rates dataframe or matrix with the observed rates, stacked per state (see Details) 
#' @param multiple allow for multiple transitions between states during an interval (default TRUE)
#' @param absorbing_state Indicate which state, if any, is the absorbing state (currently not used)
#' @return list of Mx matrices
#' @keywords manip
#' @export
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


#' @title Calculate age-specific state transfer probabilities Px.
#' 
#' @description \code{transfer_prob} returns a list of age-specific transfer probabilities matrices Px,
#' based on observed transfer rates (death and outmigration) structured in the transition matrices Mx. 
#' It allows for single or multiple transitions per interval.
#' @param transfer_rates list of transfer matrices Mx, orderd by age x
#' @param multiple allow for multiple transitions between states during an interval (default TRUE)
#' @return list of Px matrices
#' @keywords manip
#' @export 
transfer_prob <- function(transfer_rates, multiple=TRUE) {
  
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


#' @title Calculate lx, the expected number of survivors at exact age x
#' 
#' @description \code{expected_survivors} returns a list with \code{lx} matrices, containing the expected number of survivors 
#' at exact age x, based on the list of transition probability matrices \code{Px}.
#' 
#' @param transition_probs list of transition probabilities Px, ordered by age x 
#' @param radix radix (default 100000)
#' @return list of lx matrices
#' @keywords manip
#' @export
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


#' @title Calculate Lx, number of years lived in each state by unit birth cohort
#' 
#' @description \code{years_lived} returns a list with \code{Lx} matrices, containing the number of years lived 
#' in each state by unit birth cohort, based on the number of survivors \code{lx} and the transfer matrices \code{Mx}.
#' 
#' @param expected_survivors list of number of survivor matrices \code{lx}, ordered by age x 
#' @param transfer_rates list of transfer-rate matrices Mx, orderd by age x
#' @return list of Lx matrices
#' @keywords manip
#' @export
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


#' @title Calculate Sx, the proportion of survivors for each age x and state
#' 
#' @description \code{survivor_prop} returns a list with \code{Sx} matrices, containing the proportion of survivors,
#' based on the years lived in each state \code{Lx}.
#' 
#' @param years_lived list of years lived matrices \code{Lx}, ordered by age x 
#' @return list of Sx matrices
#' @keywords manip
#' @export
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


#' @title Calculate Bx, the proportion of births for each age x and state
#' 
#' @description \code{birth_prop} returns a list with \code{Bx} matrices containing the proportion of births,
#' based on the birth rates \code{Fx}, the transition probabilities \code{Px} and the survivorship proportion \code{Sx}.
#' 
#' @param birth_rates a vector of age-specific birth rates, stacked per state 
#' @param transition_probs list of transtion probability matrices Mx, orderd by age x
#' @param survivor_prop list of survivorship proportions matrices Sx, orderd by age x
#' @return list of Bx matrices
#' @keywords manip
#' @export
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


#' @title Construct the generalized Leslie matrix G
#' 
#' @description \code{projection_matrix} combines the birth (\code{Bx}) and survivorship (\code{Sx}) proportions into
#' a generalized Leslie matrix \code{G}.
#' 
#' @param birth_prop list of birth proportions matrices Bx, orderd by age x
#' @param survivor_prop list of survivorship proportions matrices Sx, orderd by age x
#' @return matrix G
#' @keywords manip
#' @export
projection_matrix <- function(birth_prop, survivor_prop) {
  # Construction a generalized Leslie matrix/projection matrix G
  
  B <- birth_prop
  S <- survivor_prop
  
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
