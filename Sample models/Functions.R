
#---------------------------------------------------------------------------#
#### R function to sample states for multiple individuals simultaneously ####
#---------------------------------------------------------------------------#

# Krijkamp EM, Alarid-Escudero F, Enns EA, Jalal HJ, Hunink MGM, Pechlivanoglou P. 
# Microsimulation modeling for health decision sciences using R: A tutorial. 
# Med Decis Making. 2018;38(3):400–22. https://www.ncbi.nlm.nih.gov/pubmed/29587047

################################################################################
# Copyright 2017, THE HOSPITAL FOR SICK CHILDREN AND THE COLLABORATING INSTITUTIONS. 
# All rights reserved in Canada, the United States and worldwide.  
# Copyright, trademarks, trade names and any and all associated intellectual property 
# are exclusively owned by THE HOSPITAL FOR SICK CHILDREN and the collaborating 
# institutions and may not be used, reproduced, modified, distributed or adapted 
# in any way without appropriate citation.

################################################################################
# Developed by Petros Pechlivanoglou

samplev <- function(m.Probs, m) {
# Arguments
 # m.Probs: matrix with probabilities (n.i * n.s)
 # m:       number of states than need to be sampled per individual  
# Return
  # ran:    n.i x m matrix filled with sampled health state(s) per individual
  
  d <- dim(m.Probs)  # dimensions of the matrix filled with the multinomical probabilities for the health states 
  n <- d[1]          # first dimension - number of rows (number of individuals to sample for)
  k <- d[2]          # second dimension - number of columns (number of health states considered)
  lev <- dimnames(m.Probs)[[2]]  # extract the names of the health states considered for sampling
  if (!length(lev))  # in case names for the health states are missing, use numbers to specify the health states
    lev <- 1:k       # create a sequence from 1:k (number of health states considered)
  # create a matrix 
  ran <- matrix(lev[1], ncol = m, nrow = n) # create the matrix ran, filled with the first health state of the levels 
  U <- t(m.Probs)    # transposed m.Probs matrix n.i x n.s --> n.s x n.i 
  
  for(i in 2:k) {    # start loop, from the 2nd health states
    U[i, ] <- U[i, ] + U[i - 1, ] # start summing the probabilities of the different health states per individual 
  }
  if (any((U[k, ] - 1) > 1e-05))  # sum of all probs per individual - 1 should be 0 (use 1e-05 for rounding issues), else print the error statement
    stop("error in multinom: probabilities do not sum to 1")
  
  for (j in 1:m) {   # start loop of the state that needs to be sampled (m)
    un <- rep(runif(n), rep(k, n))       # sample from a uniform distribution of length n*k
    ran[, j] <- lev[1 + colSums(un > U)] # store the health state at the jth column of the U matrix
  }
  ran # return the new health state per individual n.i x m
} # close the function 

# plot density of total cost
plot_tc <- function(tc) {
  # Histogram showing variability in individual total costs
  plot(density(tc), main = paste("Total cost per person"), xlab = "Cost ($)")
}

# plot density of total QALYs
plot_te <- function(te) {
  # Histogram showing variability in individual total QALYs
  plot(density(te), main = paste("Total QALYs per person"), xlab = "QALYs")
}

# plot health state trace
plot_m_TR <- function(m_M) {
  # plot the distribution of the population across health states over time (trace)
  # count the number of individuals in each health state at each cycle
  m_TR <- t(apply(m_M, 2, function(x) table(factor(x, levels = v_n, ordered = TRUE)))) 
  m_TR <- m_TR / n_i                                       # calculate the proportion of individuals 
  colnames(m_TR) <- v_n                                    # name the rows of the matrix
  rownames(m_TR) <- paste("Cycle", 0:n_t, sep = " ")       # name the columns of the matrix
  # Plot trace of first health state
  plot(0:n_t, m_TR[, 1], type = "l", main = "Health state trace", 
       ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
  # add a line for each additional state
  for (n_states in 2:length(v_n)) {
    lines(0:n_t, m_TR[, n_states], col = n_states)  # adds a line to current plot
  }
  legend("topright", v_n, col = 1:4,    # add a legend to current plot
         lty = rep(1, 3), bty = "n", cex = 0.65)
  
}

#-----------------------------------------------------------------------------------------------#
#### R function to extract the parameters of a beta distribution from mean and st. deviation ####
#-----------------------------------------------------------------------------------------------#
#' @param m mean 
#' @param s standard deviation
#' 
betaPar <- function(m, s) 
{
  a <- m * ((m * (1 - m) / s ^ 2) - 1)
  b <- (1 - m) * ((m * (1 - m) / s ^ 2) - 1)
  list(a = a, b = b)
}

#-------------------------------------------------------------------------------------------------#
#### R function to extract the parameters of a gamma distribution from mean and st. deviation  ####
#-------------------------------------------------------------------------------------------------#
#' @param m mean 
#' @param s standard deviation
#' 
gammaPar <- function(m, s) {   
  # m: mean  
  # s: standard deviation 
  shape <- m ^ 2 / s ^ 2
  scale <- s ^ 2 / m
  list(shape = shape, scale = scale)
}

