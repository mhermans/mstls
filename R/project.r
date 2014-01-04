
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


