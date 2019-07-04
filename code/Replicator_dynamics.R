integrate_dynamics <- function(x0 = NULL, # starting value
                               pars, # model paramters
                               int_function = zero_sum_phenotypes_variation, # integration function to use
                               int_time = 500, # length of intergration
                               int_steps = 1000, # number of integration steps
                               int_method = "lsoda", # can also be ode45, which is more precise but takes much longer
                               thresh = .Machine$double.eps*10) # threshold for extinction, can use zero
							   								   # but this induces computional noise near 
							   								   # machine epsilon
{

  # initial starting condition	
  if(is.null(x0)){
    x0<-runif(nrow(pars$H))
  }
  # add uniform death and fecundity rates if not supplied	
  if(is.null(pars$f)){
  	pars$f <- rep(1,nrow(pars$H))
  }
  if(is.null(pars$d)){
  	pars$d <- rep(1,nrow(pars$H))
  }
  # scale to sum to 1	
  x0<-x0/sum(x0)
  # get the integration times
  times <- seq(0, int_time, length=int_steps)
  # set the threshold for pruning
  pars$THRESH<-thresh
  # integrate!
  out <- as.matrix(ode(x0, times, int_function, pars, method = int_method))
  # clean up the output
  out<-out[!is.na(rowSums(as.matrix(out[,2:ncol(out)]))),] 
  # threshold once again to remove noise near machine epsilon
  out[out<thresh]<-0
  return(out)
}


## generalized model allowing for variable fitness. set d=f=(1,,,1)' for traditional model
zero_sum_phenotypes_variation<- function(time, x, params){
  with(as.list(params), {
    x[x<THRESH]<-0 # prune those below machine epsilon
    x<-x/sum(x) # rescale to avoid numerical errors
    dxdt<-sum(d*x)*(2 * Q%*%diag((f*x/sum(f*x)))%*%H%*%(f*x/sum(f*x))) - d*x
    return(list(dxdt))
  })
}

