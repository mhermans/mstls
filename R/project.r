
#' @title Project an intial multistate population forward using G
#' 
#' @description \code{project} takes a intitial multistate population vector and projects it forward using the generalized
#' Leslie matrix \code{G}.
#' 
#' @param init_pop initial population vector, stacked per state
#' @param n_states number of states
#' @param pmat projection matrix/generalized Leslie matrix G
#' @param n_steps number of steps to project forward
#' @param lbls optional labels for the age groups
#' @return matrix with the stacked population vectors in the columns
#' @keywords manip
#' @export
project <- function(init_pop, n_states, pmat, n_steps, lbls=NULL) {
  # linear, time-invariant projection
  
  n_ages <- length(init_pop) / n_states
  
  #split and interleave population vector
  N0 <- split(init_pop, rep(seq(n_states), each=length(init_pop) / n_states))
  N0 <- ggplot2:::interleave(N0)
  N.projected <- matrix(nrow= n_steps + 1, ncol= n_ages * n_states )
  
  colnames(N.projected) <- rep(lbls, n_states)
  rownames(N.projected) <- paste('t', seq(0, n_steps), sep='')
  
  N.projected[1,] <- N0 # assign pop at t0 to matrix
  
  
  # iterate for intermediate values, else: A %^% nsteps %*% n0
  i <- n_steps
  n <- N0
  while (i > 0) {
    n <- pmat %*% n
    N.projected[n_steps + 2 - i,] <- n
    i <- i - 1
  }
  
  # reorder with intervals as columns
  N.list <- list()
  N.projected <- t(N.projected)
  for (i in 1:n_states) {
    N.list[[i]] <- N.projected[seq(i, nrow(N.projected), n_states),]
  }
  
  # stack states
  N.projected <- do.call(rbind, N.list)
  
  N.projected
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


