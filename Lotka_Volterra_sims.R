rm(list=ls())
library(tidyverse)
library(deSolve)
library(rootSolve)

# sequence of rho values 
rhoseq <- sort(c(seq(0,1,length=25),0.999))
# sequence of p values
pseq <- sort(c(seq(0.2,1,length=25),0.999))
# fixed n and m 
nspp <- 30
mper <- 5
# how much to perturb parameters, as proportion
pertval <- 0.20
# threshold for extinction  
THRESH <- 1e-7
# max death rate, as prop of growth
max_d <- 0.1

# create the list of simulation combinations
sim_opts <- expand.grid(rho = rhoseq, p = pseq, pert = pertval, m = mper) %>% mutate(no_ph = ifelse(rho==1 & p==1,1,0))#, data.frame(rho = 0, p = 0, pert  = rep(pertval,5), m = 1, no_ph = 1))


#### AUX functions
make.comp<-function(X){
  X[lower.tri(X)]<-1-t(X)[lower.tri(t(X))]
  return(X)  
}


# calculate a membership vector from the m=(m1,m2,...mn) vector
get_membership<-function(m){
  
  ## create an expanded vector of memberships
  membership<-NULL
  n<-length(m)
  for(i in 1:n){
    membership<-c(membership,rep(i,m[i]))
  }
  
  return(membership)
}


# give it a membership vector, a rho between zero and one, and should it have 0.5 between species?
construct_A_matrix<-function(m,rho=0){
  n<-length(m)
  nm<-sum(m)
  if(rho<0 | rho>1){
    stop("rho must be between zero and one")
  }
  # create expanded A matrix
  avec <- data.frame(expand.grid(rep(1:n,m), rep(1:n, m)), val=NA) %>% setNames(c("row","col","coef"))
  
  for(i in 1:n){
    for(j in 1:n){
      # select the average value for this block
      rho_use<-rho
      # first see if we are in the diagonal, if so use rho_within
      s1<-runif(1)
      #get the maximum dist from the boundary
      mv<-max(1-s1,s1)
      if(i!=j){
        avec$coef[avec$row==i & avec$col==j] <- runif(nrow(avec[avec$row==i & avec$col==j,]),max(0,s1-mv*(1-rho_use)),min(1,s1+mv*(1-rho_use)))
      }
      else{
        # else take a weighted average between 1 and unif(0.95,1)
        avec$coef[avec$row==i & avec$col==j] <- runif(length(avec$coef[avec$row==i & avec$col==j]), 0.9,1)
      }
    }
    
  }
  # wrap into matrix
  A <- avec %>% bind_cols(expand.grid(1:nm, 1:nm)) %>% select(-row, -col) %>% spread(Var2, coef) %>% select(-Var1) %>% as.matrix()
  # diag(A) <- 1
  return(-A)
}



# construct a Q matrix with variable entries if desired
construct_Q_matrix<-function(m,p,range_p=0.1,range_q=0.1,min_above=0.05,max_p=0.95,ignore_min=FALSE){
  membership<-get_membership(m)
  max_p<-max(max_p,p)
  n<-length(m)
  nm<-sum(m)
  # set the lower bound for the diagonals
  min_vals<-rep(0,nm)
  if(!ignore_min){
    min_vals<-rep(1/m,m)
  }
  lower<-rep(p,nm)-range_p/2
  low_lim<-min_vals+min_above
  lower[lower<low_lim]<-low_lim[lower<low_lim]
  # set the upper bounds
  upper<-lower+range_p
  upper[upper>max_p]<-max_p
  lower[lower>max_p]<-max_p
  # check to see if there is only one species
  upper[rep(m,m)==1]<-lower[rep(m,m)==1]<-1
  # now sample the diagonals from a uniform dist
  Q<-matrix(0,nm,nm)
  for(i in 1:nm){
    Q[i,i]<-ifelse(upper[i]-lower[i]==0,upper[i],runif(1,lower[i],upper[i]))
  }
  # warning message in case a condition was missed
  if(any(diag(Q)<min_vals) & !ignore_min){
    stop("diagonal is less than 1/m")
  }
  # now add the off diagonals
  for(i in 1:nm){
    if(diag(Q)[i]<1){
      # relatives of the diagonal species
      relatives<-membership==membership[i] & (1:nm)!=i
      # equal weighted q
      qvals<-(1-diag(Q)[i])/(rep(m,m)[i]-1)
      if(range_q>0){
        qvals<-2
        count<-0
        # sample then normalize
        while(any(qvals>diag(Q)[i]) & !ignore_min){
          if(count>20){
            qvals<-(1-diag(Q)[i])/(rep(m,m)[i]-1)
            break()
          }  
          qvals<-runif(sum(relatives),max(0,qvals-range_q),qvals+range_q)
          qvals<-qvals/sum(qvals)*(1-diag(Q)[i])
          count<-count+1
        }
      }
      if(any(qvals>diag(Q)[i]) & !ignore_min){
        stop("One of the off diagonals is larger than the diagonal")
      }
      Q[relatives,i]<-qvals
    }
  }
  return(Q)
}

# LV dynamics with mutation
LV_Q<- function(time, x, params){
  with(as.list(params), {
    x[x<THRESH]<-0
    dxdt<- Q%*%(r*x)+x*(A%*%x)-d*x
    return(list(dxdt))
  })
}

# output matrix
res_mat <- matrix(NA, nrow=nrow(sim_opts), ncol=5)

# loop through all the simulation configurations
for(ii in 1:nrow(sim_opts)){
  print(ii)
  n<-nspp
  m<-rep(sim_opts$m[ii], n)
  n0 <- length(m)
  m0 <- sum(m)
  A <- construct_A_matrix(m,rho=sim_opts$rho[ii])
  
  # get get Q matrix
  if(sim_opts$p[ii]==1){
    Q <- diag(m0)
  }
  else{
    Q <- construct_Q_matrix(m,p=sim_opts$p[ii],range_p = 0.15,range_q = 0.15,ignore_min = F)      
  }
  # get growth rate, perturbed from species average
  r <- rep(runif(n),m)*sim_opts$rho[ii]+runif(sum(m))*(1-sim_opts$rho[ii])
  d <- r*(rep(runif(n, 0, max_d),m)*sim_opts$rho[ii]+runif(sum(m),0,max_d)*(1-sim_opts$rho[ii]))
  # }
  # integrate dynamics
  x0 <- rep(runif(n),m)
  x0 <- x0/sum(x0)*n
  dyn <- ode(y = x0,times = seq(0,1e7, length=10), parms = list(r=r, A=A, Q=Q, d=d, THRESH = 1e-9), func = LV_Q, method = "lsoda")
  # endpoint
  xs <- as.numeric(dyn[nrow(dyn),-1])
  
  # get phenotypes and species alive
  alive <- xs>THRESH
  mlabs <- rep(1:n,m)
  mlabs <- mlabs[alive]
  sp_alive <- unique(mlabs)
  sp_labs <- data.frame(ph = 1:ncol(A), sp = rep(1:n,m))
  ph_alive <- sp_labs %>% filter(sp%in%sp_alive)
  # subset, get new community size
  A <- matrix(A[ph_alive$ph,ph_alive$ph],nrow(ph_alive), nrow(ph_alive))
  r <- r[ph_alive$ph]
  d <- d[ph_alive$ph]      
  Q <- matrix(Q[ph_alive$ph, ph_alive$ph],nrow(ph_alive), nrow(ph_alive))
  xs0 <- xs[ph_alive$ph]
  m <- as.numeric(table(ph_alive$sp))
  n <- length(m)
  # initialize summary stats
  lam1 <- n_dead <- NA
  
  if(any(as.numeric(LV_Q(time = 0, x = xs0, params = list(r=r, A=A, Q=Q, d=d, THRESH = 1e-9))[[1]])>1e-8)){
    
    dyn <- ode(y = xs0,times = seq(0,1e7, length=10), parms = list(r=r, A=A, Q=Q, d=d, THRESH = 1e-9), func = LV_Q, method = "lsoda")
    # endpoint
    xs <- as.numeric(dyn[nrow(dyn),-1])
    
    # get phenotypes and species alive
    alive <- xs>THRESH
    mlabs <- rep(1:n,m)
    mlabs <- mlabs[alive]
    sp_alive <- unique(mlabs)
    sp_labs <- data.frame(ph = 1:ncol(A), sp = rep(1:n,m))
    ph_alive <- sp_labs %>% filter(sp%in%sp_alive)
    # subset, get new community size
    A <- matrix(A[ph_alive$ph,ph_alive$ph],nrow(ph_alive), nrow(ph_alive))
    r <- r[ph_alive$ph]
    d <- d[ph_alive$ph]      
    Q <- matrix(Q[ph_alive$ph, ph_alive$ph],nrow(ph_alive), nrow(ph_alive))
    xs0 <- xs[ph_alive$ph]
    m <- as.numeric(table(ph_alive$sp))
    n <- length(m)
    # initialize summary stats
    lam1 <- n_dead <- NA
  }
  
  # if more than 2 species, perturb the parameters and integrate
  if(n>=2){
    ############# Jacobian stability
    J <- jacobian.full(xs0, func = LV_Q, time=0, parms = list(r=r, A=A, Q=Q, THRESH = 1e-9))
    # get max real eigen
    lam1 <- max(Re(eigen(J)$values))   
    
    ### survival after pertubation
    r1 <- r*(rep(runif(n,1-pertval, 1+pertval),m)*sim_opts$rho[ii]+runif(length(r),1-pertval, 1+pertval)*(1-sim_opts$rho[ii]))
    d1 <- d*(rep(runif(n,1-pertval, 1+pertval),m)*sim_opts$rho[ii]+runif(length(d),1-pertval, 1+pertval)*(1-sim_opts$rho[ii]))
    
    # integrate dynamcis
    dyn2 <- ode(y = xs0,times = seq(1,1e7, length=5), parms = list(r=r1, A=A, Q=Q, d=d1, THRESH = 1e-9), func = LV_Q, method = "lsoda")
    xs2 <- as.numeric(dyn2[nrow(dyn2),-1])
    # get number dead
    n_dead <- n-length(unique(rep(1:n,m)[xs2>THRESH]))
    
  }
  res_mat[ii,] <- c(n0, m0, n, lam1, n_dead)
}
