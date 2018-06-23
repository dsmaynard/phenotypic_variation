
# this function constructs the matrix of competitive abilities H
# each row/column of the matrix represents a phenotype.
# The phenotypes are ordered by species, such that all phenotypes
# belonging to the same species are adjacent.
# The input parameters are:
# m  ->  Membership vector, e.g., m = c(1, 3, 4). 
#        Each coefficient specified the number of phenotypes in 
#        the corresponding species. In the example, the first species
#        has a single phenotype, the second has three phenotypes
# tau -> Parameter controlling the similarity among the phenotypes
#        of a species. tau = 0 means phenotypes are uncorrelated;
#        tau = 1 that they are identical.
# tau_within -> parameter controlling the similarity within the 
#        diagonal block of each species; this controls how 
#        competition among the phenotypes of the same species
#        is regulated.
construct_matrix_H <- function(m, 
                               tau = 0, 
                               tau_within = NULL
                               ){
  n <- length(m) # number of species
  nm <- sum(m) # number of phenotypes (in total)
  if (tau < 0.0 | tau > 1.0) {
    stop("tau must be between zero and one")
  }
  if (is.null(tau_within)) {
    tau_within <- tau # if not specified use tau for the whole matrix
  }
  H <- NULL 
  for(i in 1:n){
    rowmat <- NULL
    zrow <- NULL
    for (j in 1:n) {
      # select the average value for the interactions between species i and j
      tau_use <- tau
      # if this is the intra-specific coefficient
      if (i == j) {
        tau_use <- tau_within
        s1 <- 0.5
      }
      else{
        s1 <- runif(1)
      }
      # get the distance from the boundary
      mv <- max(1 - s1, s1)
      # sample uniformly, with the bounds scaled by tau_use
      # this takes the block of matrix H corresponding to the
      # interactions between the m[i] phenotypes of species i
      # and the m[j] phenotypes of species j
      subblock <- matrix(runif(m[i] * m[j], 
                               max(0, s1 - mv * (1 - tau_use)),
                               min(1, s1 + mv * (1-tau_use))),
                         m[i], m[j])
      # bind the new columns to the previous ones
      rowmat <- cbind(rowmat, subblock)
    }
    # bind the new rows to the previous ones
    H <- rbind(H, rowmat)
  }
  # now that we have set the values, normalize them such that
  # H + H^t is a matrix of all ones
  H[lower.tri(H)] <- 1 - t(H)[lower.tri(t(H))]
  diag(H) <- 0.5 # each phenotype has probability 1/2 of winning with themselves
  return(H)
}

# construct a Q matrix with variable entries if desired
construct_matrix_Q<-function(m,
                             p,
                             range_p = 0.0,
                             range_q = 0.0,
                             min_above = 0.05,
                             max_p = 0.95, 
                             ignore_min=FALSE){
  n <- length(m) # number of species
  nm<-sum(m) # total number of phenotypes
  # create an expanded vector of memberships
  membership <- numeric(0)
  for (i in 1:n) {
    membership <- c(membership, rep(i, m[i]))
  }
  # set maximum phenotypic memory
  max_p <- max(max_p, p)
  # set the lower bound for the diagonal elements
  min_vals <- rep(0, nm)
  if (!ignore_min) {
    min_vals <- rep(1 / m, m)
  }
  lower <- rep(p, nm) - range_p / 2
  low_lim <- min_vals + min_above
  lower[lower < low_lim] <- low_lim[lower < low_lim]
  # set the upper bounds
  upper <- lower + range_p
  upper[upper > max_p] <- max_p
  lower[lower > max_p] <- max_p
  # check to see if there is only one species
  upper[rep(m, m) == 1] <- 1
  lower[rep(m, m) == 1] <- 1
  # now sample the diagonal elements from a uniform distribution
  Q <- matrix(0, nm, nm)
  for (i in 1:nm) {
    Q[i, i]<- ifelse(upper[i] - lower[i] == 0, upper[i], runif(1, lower[i], upper[i]))
  }
  # warming message in case a condition was missed
  if (any(diag(Q) < min_vals) & !ignore_min) {
    stop("diagonal is less than 1/m")
  }
  # now add the off diagonals
  for (i in 1:nm) {
      if (diag(Q)[i] < 1){
      # relatives of the diagonal species
      relatives <- membership == membership[i] & (1:nm)!=i
      # equal weighted q
      qvals <- (1 - diag(Q)[i]) / (rep(m, m)[i] - 1)
      if(range_q > 0){
        qvals <- 2
        count <- 0
        # sample, then normalize
        while (any(qvals > diag(Q)[i]) & !ignore_min) {
          if (count > 20) {
            qvals <- (1 - diag(Q)[i]) / (rep(m, m)[i] - 1)
            break()
          }	
          qvals <- runif(sum(relatives), max(0, qvals - range_q), qvals + range_q)
          qvals <- qvals / sum(qvals) * (1 - diag(Q)[i])
          count <- count+1
        }
      }
      if(any(qvals > diag(Q)[i]) & !ignore_min){
        stop("One of the off diagonals is larger than the diagonal")
      }
      Q[relatives, i] <- qvals
    }
  }
  if(!all(round(colSums(Q), 12) == 1)){
    stop("Q is not stochastic, something went wrong...")
  }
  return(Q)
}

