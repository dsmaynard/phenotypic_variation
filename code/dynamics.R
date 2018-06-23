integrate_dynamics <- function(x0 = NULL,
                               pars,
                               func_use = zero_sum_phenotypes,
                               maxtime = 500,
                               maxsteps = 10000,
                               lengthtime = 1000,
                               thresh = 1e-7, 
                               method = "ode45", 
                               rtol = NULL, 
                               atol = NULL){
  if(is.null(x0)){
    x0<-runif(nrow(pars$H))
  }
  x0<-x0/sum(x0)
  times <- seq(0, maxtime, length=lengthtime)
  pars$THRESH<-thresh
  if(!is.null(rtol) & !is.null(atol)){
    out <- as.matrix(ode(x0, times, func_use, pars, method = method, maxsteps = maxsteps,
                         rtol = rtol, atol = atol))  
  }
  else{
    out <- as.matrix(ode(x0, times, func_use, pars,method = method, maxsteps = maxsteps))
  }
  out<-out[!is.na(rowSums(as.matrix(out[,2:ncol(out)]))),] 
  out[out<thresh]<-0
  return(out)
}

zero_sum_phenotypes<- function(time, x, params){
  with(as.list(params), {
    x[x < THRESH] <- 0 # extinct phenotypes
    x <- x / sum(x) # keep on the simplex
    dxdt<- 2 * diag(as.numeric(Q %*% x)) %*% H %*% (Q %*% x) - x
    return(list(dxdt))
  })
}