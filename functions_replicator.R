library(RColorBrewer)
library(scales)
library(deSolve)
library(tidyverse)
library(reshape2)
library(gridExtra)
library(igraph)
library(gridExtra)
library(Rglpk)



#######################################
# various replicator dynamical models


# the main mutator equation, which collapses to the standard eqn when Q=I
zero_sum_mutator<- function(time, x, params){

	with(as.list(params), {

		x[x<THRESH]<-0

		x<-x/sum(x)

		dxdt<-2 * Q%*%diag(as.numeric(x))%*%H%*%(x) - x

		return(list(dxdt))
	})
}


## express in terms of "seed output" of the parent
zero_sum_phenotypes<- function(time, x, params){

	with(as.list(params), {

		x[x<THRESH]<-0

		x<-x/sum(x)

		dxdt<-2 * diag(as.numeric(Q%*%x))%*%H%*%(Q%*%x) - x

		return(list(dxdt))
	})
}

## allowing for variable fitness
zero_sum_mutator_fitness<- function(time, x, params){

	with(as.list(params), {

		x[x<THRESH]<-0

		x<-x/sum(x)


		dxdt<-sum(d*x)*(2 * Q%*%diag((f*x/sum(f*x)))%*%H%*%(f*x/sum(f*x))) - d*x

		return(list(dxdt))
	})
}

# calculate the instantaneous rate of change. used for checking the tolerance solution
point_dxdt<- function(x,H,Q){
	
	dxdt<-2*Q%*%diag(x)%*%H%*%x - x

	return(dxdt)
}

point_dxdt_seed<- function(x,d,f,H,Q){
	
	dxdt<-2 * diag(as.numeric(Q%*%x))%*%H%*%(Q%*%x) - x

	return(dxdt)
}

point_fitness<- function(x,d,f,H,Q){
	
	dxdt<-sum(d*x)*(2 * Q%*%diag((f*x/sum(f*x)))%*%H%*%(f*x/sum(f*x))) - d*x

	return(dxdt)
}

# return the average error for the optimization problem
soln_error<- function(x,H,Q,scale_v=1000){
	
	# x<-x/sum(x)*scale_v
	val<-sqrt(sum(point_dxdt(x,H,Q)^2))

	return(val)
}


# generate a random x
getx<-function(D){
	
	x<-exp(D)
	x<-x/sum(x)
	return(x)
}

# calculate the Jacobian
Jacobian_H<-function(x,Q,H){
	D<-nrow(Q)
	if(D>1){
		2*Q%*%(diag(x)%*%(H)+diag(as.numeric((H)%*%x)))-diag(D)
	}
	else{
		return(as.matrix(1))	
	}
}

Jacobian_P<-function(x,Q,H){
	P<-H-t(H)
	D<-nrow(Q)
	2*Q%*%(diag(x)%*%(P)+diag(as.numeric((P)%*%x)))+Q-diag(D)
}

######################################################
## integrating, plotting, checking dynamical results


integrate_dynamics<- function(x0=NULL,pars,func_use,maxtime=1000,maxsteps=10000,lengthtime=100,thresh=1e-7,method="ode45",rtol=NULL,atol=NULL){

	if(is.null(x0)){
		x0<-runif(nrow(pars$H))
	}
	
	x0<-x0/sum(x0)
	
	times <- seq(0, maxtime, length=lengthtime)

	pars$THRESH<-thresh
	
	if(!is.null(rtol) & !is.null(atol)){
		out <- as.matrix(ode(x0, times, func_use, pars,method=method, maxsteps = maxsteps,rtol=rtol,atol=atol))  
	}
	else{
		out <- as.matrix(ode(x0, times, func_use, pars,method=method, maxsteps = maxsteps))
	}
	
	out<-out[!is.na(rowSums(as.matrix(out[,2:ncol(out)]))),] 

	out[out<thresh]<-0
	return(out)
}

# plotting dynamics, as well as calculating the final commmunity composition
# this can also "collapse" across phenotypes of the same species to get
# species-level averages
plot_dynamics<-function(comdyn,m,xlim=NULL,ylim=NULL,collapse=FALSE,avg_prop=NULL,show_plot=T,normalize=T,return_grob=F,trans="identity"){

	n<-length(m)
	nm<-sum(m)
	
	membership<-get_membership(m)

	# group the species?
	if(collapse){

		## add the abundances for each species
		cd2<-comdyn[,1]
		comdyn<-comdyn[,-1]
		for(i in 1:n){
			cd2<-cbind(cd2,apply(as.matrix(comdyn[,membership==i]),1,sum))
		}

		comdyn<-data.frame(cd2)
		names(comdyn)<-c("time",paste("X",1:n,sep=""))
	}

	# rescale to within 0,1 to account for round-off error?
	if(normalize){
		rs<-apply(comdyn[,-1],1,sum)
		comdyn[,-1]<-sweep(comdyn[,-1],1,rs,FUN="/")
	}
	
	abund_mat<- data.frame(comdyn) %>% melt(id="time") %>% dplyr::rename(species=variable)
	
	#plot? with or without limits?
	if(show_plot){		
		if(!is.null(xlim) & !is.null(ylim)){
			print(g11<-ggplot(abund_mat, aes(x=as.numeric(time), y=value+1e-16, colour=species)) + geom_line() +
				  	xlim(xlim) +ylim(ylim)+theme(legend.position = "none")+scale_y_continuous(trans=trans,name="Time")+xlab("Relative abundance"))
		}
		else if(!is.null(xlim) & is.null(ylim)){
			print(g11<-ggplot(abund_mat, aes(x=as.numeric(time), y=value+1e-16, colour=species)) + geom_line() +
				  	xlim(xlim)+theme(legend.position = "none")+scale_y_continuous(trans=trans,name="Time")+xlab("Relative abundance"))
		}	
		else if(is.null(xlim) & !is.null(ylim)){
			print(g11<-ggplot(abund_mat, aes(x=as.numeric(time), y=value+1e-16, colour=species)) + geom_line() +
				  	ylim(ylim)+theme(legend.position = "none")+scale_y_continuous(trans=trans,name="Time")+xlab("Relative abundance"))
		}	
		else{
			print(g11<-ggplot(abund_mat, aes(x=as.numeric(time), y=value+1e-16, colour=species)) + geom_line()+
				  	theme(legend.position = "none")+scale_y_continuous(trans=trans,name="Time")+xlab("Relative abundance"))
		}	
	}
	
	Nst<-as.numeric(comdyn[nrow(comdyn),-1])
	
	if(!is.null(avg_prop)){
		get_row<-round(nrow(comdyn)*avg_prop)
		if(get_row<nrow(comdyn)){
			comdyn<-comdyn[get_row:nrow(comdyn),-1]
			Nst<-as.numeric(apply(comdyn,2,mean))
		}
		else{
			Nst<-as.numeric(comdyn[nrow(comdyn),-1])
		}
	}
	
	if(return_grob){
		return(g11)
	}
	else{
		return(Nst)
	}
}


check_zero_sum<-function(results){
	if(max(abs(apply(results[,-1],1,sum)-1))>1e-2){
		print("Not all zero")
		print(paste("Max deviation = ",max(max(abs(apply(results[,-1],1,sum)-1)))))
	}
	else{
		print("Good")
	}
}


# use the solution to the Q=I scenario to seed the models
get_starting_x<-function(m,H){
	
	B<-make_species_sum_matrix(m)
	xs<-B%*%find_optimal_strategy(H,verbose=F)
	xs[xs==0]<-max(min(xs),.01)
	xs<-(xs/sum(xs)*1/m)
	xs<-t(B)%*%xs
	xs<-xs/sum(xs)
	
	return(as.numeric(xs))
}


#############################
# constructing matrices



make_Q_payoff<-function(Q){

	diag(Q)<--(1-diag(Q))

	return(Q)
}



make.symm<-function(X){
	X[lower.tri(X)]<-t(X)[lower.tri(t(X))]
	return(X)	
}



make.comp<-function(X){
	X[lower.tri(X)]<-1-t(X)[lower.tri(t(X))]
	return(X)	
}


# this creates a matrix which sums all of values within species
make_species_sum_matrix<-function(m){
	
	nm<-sum(m)
	n<-length(m)
	
	B<-matrix(0,n,nm)
	
	B[1,1:m[1]]<-1
	end<-m[1]
	for(i in 2:n){
		start<-end+1
		end<-start+m[i]-1
		B[i,start:end]<-1
	}

	return(B)
}



#########################################
# solving for optimal solution using linear programming, assuming uncorrelated phenotypes

find_optimal_strategy<- function(H,verbose=TRUE,time_limit=5000){
    n <- dim(H)[1]
    f.obj <- rep(1, n)
    In<-diag(n)
    f.con <- H
    f.rhs <- rep(1, n)
    f.dir <-rep("<=", n)
    z<-Rglpk_solve_LP(obj=f.obj,mat=f.con,dir=f.dir,rhs=f.rhs,max=TRUE,control=list(tm_limit=time_limit,presolve=TRUE,verbose=verbose))
    return(z$solution / sum(z$solution))
}


############################################
## calculating the theoretical probability of survival

prob_m_given_k<-function(m,k,s,N){
	tsum<-0

	for(i in 0:m){
		tsum<-tsum+(-1)^i*exp(log(choose(m,i))+log(choose(s*(m-i),k))-log(choose(N*s,k))+log(choose(N,m)))
	}

	return(as.numeric(tsum))
}


prob_k_given_Nst<-function(k,s,N){

	if(k%%2==1){

		tval<-choose(s*N,k)*2^(-N*s+1)

		return(tval)
	}
	else{
		return(0)
	}

}


prob_m<-function(m,s,N){
	D<-N*s
	tsum<-0
	for(k in 1:D){
		if(k%%2==1){
			for(i in 0:m){
				
				tsum<-tsum+2^(-(N*s)+1)*(-1)^i*choose(N,m)*choose(m,i)*choose(s*(m-i),k)
			}
		
		}
	}

	return(tsum)
}


# get the theoretical expecation for the number of species that 
# should survive, given on average m_avg strains per species
get_dist_survive<-function(n,m_avg){


	val<-rep(0,n)
	for(k in 1:n){
		val[k]<-as.numeric(prob_m(k,m_avg,n))
	}

	return(data.frame(richness=1:n,expected=val))
}

get_expected_value<-function(n,m){
	
	vals<-get_dist_survive(n,m)
	
	if(abs(sum(vals$expected)-1)>1e-8){
		stop("probabilties don't add up")
	}
	
	return(sum(vals$richness*vals$expected))
}

	
get_stable_subset<-function(s,N,xm,H,Q){
	
	# get the species level average
	xm_c<-make_species_sum_matrix(s,N)%*%xm
	
	alive<-rep(xm_c>0,each=s)
	
	Q<-Q[alive,alive]
	H<-H[alive,alive]
	
	xm_c<-xm_c[xm_c>0]
	xm_c<-xm_c/sum(xm_c)
	xm<-xm[alive]
	xm<-xm/sum(xm)
	
	return(list(H=H,Q=Q,xm=xm,xm_c=xm_c,N=sum(xm_c>0)))
}
############################################################
# constructing tournament matrices
	


	 	
# generate a vector of the number of phenotypes per species. 
# mean is the target value, range is the spred (uniform) around the mean
# optionally, nm constrains the total number of phenotypes in the community
# and it will search for a solution with exactly sum(m)=nm

assign_phenotypes<-function(n,nm=NULL,mean,range=1){
	
	
	# get the min and max, truncated at 1
	min<-max(1,mean-range)
	max<-mean+range
	
	# sample values (or assign if min=max)
	if(min!=max){
		m<-sample(min:max,n,replace=T)
	}
	else{
		m<-rep(min,n)
	}
	
	# see if we want to ensure that the num. of phenotypes adds up to a fixed value
	if(!is.null(nm)){
		if(nm<(n*min)){
			stop("can't have more fewer phenotypes than min*n")
		}
		else if(nm>(n*max)){
			stop("can't have more phenotypes than max*n")
		}
		else{
			count<-0
			while(sum(m)!=nm){
				count<-count+1
				m<-sample(min:max,n,replace=T)
				if(count==10000){
					stop("something's wrong -- can't find a distribution of phenotypes that matches n*m with the given min/max values")
				}
			}
		}
	}
	
	return(m=m)
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
construct_H_matrix<-function(m,rho=0,return_zipped=F,rho_within=NULL){

	n<-length(m)
	nm<-sum(m)

	
	if(rho<0 | rho>1){
		stop("rho must be between zero and one")
	}
	
	if(is.null(rho_within)){
		rho_within<-rho
	}
	
	# H is the full phenotype matrix, Hz is the zipped-up average matrix
	H<-Hz<-NULL
	
	# adjust the resource vector so as to maintain the correlation but also have the target probability
	for(i in 1:n){
		rowmat<-zrow<-NULL
		
		# note that this calcualtes values for all entries, not just the upper triangle
		# then calls make.comp to make sure that H+t(H)=1
		for(j in 1:n){
			
			# select the average value for this block
			rho_use<-rho
			
			# first see if we are in the diagonal, if so use rho_within
			if(i==j){
				rho_use<-rho_within
				s1<-0.5
			}
			else{
				s1<-runif(1)
			}
			
			#get the maximum dist from the boundary
			mv<-max(1-s1,s1)
			
			# sample uniformly, with the bounds scaled by rho_use
			subblock<-matrix(runif(m[i]*m[j],max(0,s1-mv*(1-rho_use)),min(1,s1+mv*(1-rho_use))),m[i],m[j])

			# bind the new columns to the previous ones
			rowmat<-cbind(rowmat,subblock)
 			# get the zipped up group average
 			zrow<-cbind(zrow,matrix(mean(subblock),1,1))
		}

		# bind the new rows to the previous ones
		H<-rbind(H,rowmat)
		Hz<-rbind(Hz,zrow)			
	}

	
	# make them complementary
	H<-as.matrix(make.comp(H))
	Hz<-as.matrix(make.comp(Hz))

	diag(H)<-diag(Hz)<-.5
	
	if(return_zipped){
		return(list(H=H,Hz=Hz))
	}
	else{
		return(H)
	}

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
	
	# warming message in case a condition was missed
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

	if(!all(round(colSums(Q),12)==1)){
		stop("Q is not stochastic, something went wrong...")
	}
	
	
	return(Q)
}


### assemble a full community
build_community<-function(n,m,mean_m,range_m=0,rho=0,rho_within=NULL,p,range_p=0,range_q=0,ignore_min=F,mal_case=F,neut_case=F){
	

	if(!mal_case & !neut_case){
		if(missing(m)){
			if(missing(n) | missing(mean_m)){
				stop("Must provide both n and mean_m")
			}
			m<-assign_phenotypes(n,mean=mean_m,range=range_m) #svals$mean[i]+round(svals$mean[i]*svals$range_prop[i]))
			H<-construct_H_matrix(m,rho=rho,rho_within=rho_within)		
		}	
		else{
			if(!missing(n) | !missing(mean_m)){
				stop("Note: ignoring n and mean_m; just using m")
			}
			n<-length(m)
			H<-construct_H_matrix(m,rho=rho,rho_within=rho_within)			
		}
	}
	else{
		if(!missing(m) | !missing(mean_m) | missing(n)){
			stop("Must only give n when for maladapted and neural cases (leave m and mean_m blank)")
		}
		
		m<-rep(2,n)	
		# construct the matrices
		if(neut_case & !mal_case){
			print("NOTE: returning NEUTRAL case")
			H<-construct_H_one_neutral(n)
		}
		else if(!neut_case & mal_case){
			print("NOTE: returning MALADAPTED case")
			H<-construct_H_one_maladpated(n)	
		}
		else if(neut_case & mal_case){
			stop("Cannot have both neut_case and mal_case as TRUE)")
		}
	}
		
	nm<-sum(m)
	membership<-get_membership(m)	
	Q<-construct_Q_matrix(m,p=p,range_p = range_p,range_q = range_q,ignore_min=F)		

  return(list(H = H, Q = Q, nm = nm, n = n, m=m, membership = membership))
}



################ 
# for ode
prune_system<-function(resm,m,Q,H,membership,THRESH,show_plot=T,prop=0.5){
	
	# m<<-as.numeric(table(membership))
		
	if(show_plot){
		final<-plot_dynamics(resm,m=m,show_plot=show_plot,avg_prop = prop,normalize=F)
	}
	else{
		final<-as.numeric(resm[nrow(resm),-1])
	}
	
	keep<-final>THRESH
	
	# final.o<<-final
	# m.o<<-m
	# membership.o<<-membership
	# Q.o<<-Q
	# H.o<<-H
	
	B<-make_species_sum_matrix(m)
	
	keep.spec<-rep(B%*%final>THRESH,m)
	
	keep<-(keep & keep.spec)
	# get the new membership
	membership2<-membership[keep]

	#keep every phenotype that has at least one sibling phenotype survivng
	keep<-membership%in%membership2

	# prune
	Q<-Q[keep,keep]
	H<-H[keep,keep]
	final<-final[keep]
	membership<-membership[keep]	
	m<-as.numeric(table(membership))

	if(sum(m)!=length(final)){
		stop("m and final are different lengths")
	}
	
	return(list(Q=Q,H=Q,final=final,membership=membership,m=m))
		
}
	
	
# calculate the number of species that survive, assuming uncorrelated phenotypes
get_num_survive<-function(membership,H,verbose=FALSE,time_limit=2500){
	
	xstb<-find_optimal_strategy(H,verbose=verbose,time_limit=time_limit)

	# get the number species that survived
	surv<-membership[xstb>0]

	return(list(n=length(unique(surv)),nm=length(surv)))
}

# this takes the ode output and sums up the phenotypes to convert it to species fluctuations
convert_to_species_dynamics<-function(resm,m,THRESH){

	B<-make_species_sum_matrix(m)
	if(colnames(resm)[1]=="time"){
		timev<-resm[,1]
		resm<-resm[,-1]
	}
	else{
		timev<-NULL
	}
	
	# normalize the rows
	resm[resm<-THRESH]<-0
	resm<-diag(1/apply(resm,1,sum))%*%resm
	
	Sdyn<-NULL
	for(i in 1:nrow(resm)){
		Sdyn<-rbind(Sdyn,as.numeric(B%*%as.numeric(resm[i,])))
	}
	Sdyn<-data.frame(Sdyn)
	colnames(Sdyn)<-1:ncol(Sdyn)
	if(!is.null(timev)){
		Sdyn<-data.frame(time=timev,Sdyn)
	}
	
	return(Sdyn)
}
# this function calculates the observed species diversity and thoeretical diversity, and optionally plots
compare_obs_expected<-function(m,n,rho,div,ptsize=2,print_plot=F){
	
	# get average pheno. div
	m_avg<-round(mean(m))
	m_cv<-sd(m)/m_avg
	
	# get analytic expectation
	val<-get_dist_survive(n,m_avg)
	if(round(sum(val$expected),5)!=1){
		stop("probabilities don't sum to one")
	}
	
	# get observed values
	val<-val%>%mutate(observed=rep(0,nrow(val)))
	for(i in 1:nrow(val)){
		val$observed[i]<-sum(div>(val$richness[i]-.5) & div<(val$richness[i]+.5))/length(div)
	}

	# double check
	if(round(sum(val$expected),5)!=1){
		stop("probabilities don't sum to one")
	}	
	
	# reshape
	val<- val %>% mutate(expected_sp=get_dist_survive(n,1)$expected) %>% melt(id="richness",value.name="probability")
	
	# make plot
	legend.title<-""
	g1<-ggplot(val,aes(x=richness,y=probability,color=variable,shape=variable))+
		geom_point(size=ptsize)+geom_line()+guides(color = guide_legend(legend.title), shape = guide_legend(legend.title))+
		labs(title=paste("rho=",rho,"m avg",m_avg,", m var",m_cv,sep=""))
	
	# optional plot
	if(print_plot){
		print(g1)
	}

	# add in covariates
	val<-val %>% mutate(n=n,m_avg=m_avg,m_cv=m_cv,rho=rho)
	return(val)

}


################################
# phenotype dynamics and colors
###################################

construct_H_one_neutral <- function(n){

	# set up vector of membership, number of m, number of species
	m<-rep(2,n)
	nm <- sum(m)

	membership <- numeric(0)
	for (i in 1:n) membership <- c(membership, rep(i, m[i]))
	
	# build matrix H with random uniform values
	H <- matrix(runif(nm * nm), nm, nm)
	diag(H) <- 1/2
	H[lower.tri(H)] <- 0
	tH <- t(H)
	tH[upper.tri(H)] <- 1 - H[upper.tri(H)]
	H <- tH
	
	# now fill in the blocks
	for (i in 1:n) H[membership == i, membership == i] <- 1/2

	
	# Make one of the phenotypes neutral
	H[(1:nm)%% 2 == 0,] <- 1/2
	H[,(1:nm)%% 2 == 0] <- 1/2

	return(H)
}

construct_H_one_maladpated <- function(n,block_diag=F){

	# set up vector of membership, number of m, number of species
	m<-rep(2,n)
	nm <- sum(m)

	membership <- numeric(0)
	for (i in 1:n) membership <- c(membership, rep(i, m[i]))
	
	# build matrix H with random uniform values
	H <- matrix(runif(nm * nm), nm, nm)
	diag(H) <- 1/2
	H[lower.tri(H)] <- 0
	tH <- t(H)
	tH[upper.tri(H)] <- 1 - H[upper.tri(H)]
	H <- tH
	
	welladapted <- (1:nm) %% 2 == 1
	maladapted <- (1:nm) %% 2 == 0
	# now fill in the blocks 
	H[maladapted, welladapted] <- 0
	H[welladapted, maladapted] <- 1
	H[maladapted,maladapted] <- 1/2
	
	if(block_diag){
		for (i in 1:n) {
			H[membership == i, membership == i] <- 1/2
		}
	}
	return(H)
}


convert_to_hex<-function(x){
	rgb(red=x[1], green=x[2], blue=x[3], alpha=0.5, maxColorValue=255)
}

choose_phenotype_colors <- function(membership, nspp, pl=0.5, pd=0.5){
	# Choose colors
	# max 6 species, 9 phenotypes each
	if(nspp>6) stop("can't have more than 6 species -- not enough colors")
	
	mycolors <- character(0)
	sppcolors <- character(0)
	sppdarkcolors <- character(0)
	palettenames <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f") #c("Reds", "Blues", "Greens", "Purples", "Oranges","Greys")
	for (i in 1:nspp){
		
		npheno <- sum(membership == i)
		mid<-col2rgb(palettenames[i])
		light<-convert_to_hex(round(sqrt((pl*mid^2+(1-pl)*c(255,255,255)^2))))
		dark<-convert_to_hex(round(sqrt((pd*mid^2+(1-pd)*c(0,0,0)^2))))
		cr<-colorRampPalette(c(light,dark))
		sppdarkcolors<-c(sppdarkcolors,dark)
		tmpcol<-palettenames[i]
		if(npheno>1) tmpcol<-cr(npheno)
		mycolors <- c(mycolors, tmpcol)
		
	}
	return(list(mycolors = mycolors, sppcolors = palettenames[1:nspp], sppdarkcolors = sppdarkcolors))
}


plot_competition_graph_phenotypes <- function(
  H, # matrix where Hij = prob i beats j
  membership, # vector where membership[i] is the species phenotype i belongs to 
  n, # number of species (i.e., max(membership))
  pl=0.8, # blend of color vs. white for phenotypes
  pd=0.8, # scale of color vs. black for phenotypes
  write_svg=FALSE,
  width=4,
  height=4,
  bottom_pad=0,
  file_out
){
  # get colors
  tmp <- choose_phenotype_colors(membership, n, pl=pl, pd=pd)
  mycolors <- tmp$mycolors
  sppcolors <- tmp$sppcolors
  # Build graph
  Pay <- H - t(H)
  Pay[Pay <= 0] <- 0
  g <- graph_from_adjacency_matrix(Pay, weighted = TRUE, mode = "directed")
  V(g)$color <- mycolors
  coords <- layout_in_circle(g, order = V(g))

  g1<-plot(g, mark.groups = split(1:vcount(g), membership), layout = coords, mark.col = sppcolors, mark.border = sppcolors,
       edge.width= 3 * E(g)$weight, edge.color = "black",vertex.label=NA)
  
  if(write_svg){
  	svg(file_out,width=width,height=height)
  	par(mar=c(bottom_pad,0,0,0)+.1)
  	plot(g, mark.groups = split(1:vcount(g), membership), layout = coords, mark.col = sppcolors, mark.border = sppcolors,
       edge.width= 3 * E(g)$weight, edge.color = "black",vertex.label=NA)
  	dev.off()
  }
  return(NULL)
  
}



get_phenotype_dynamics <- function(pars, x0=NULL, random_seed = 100, tot_time = 100, step_time = 0.1,THRESH=0,nsteps=NULL,method="lsoda",pd=0.8,pl=0.8,alpha=0.8,func_use=zero_sum_phenotypes){
	# set initial conditions
	set.seed(random_seed)
	if(is.null(x0)){
		x0 <- runif(pars$nm)
		x0 <- x0 / sum(x0)
	}
	else{
		x0<-x0/sum(x0)
	}
	
	if(is.null(nsteps)){
		nsteps=round(tot_time/step_time)
	}
	# a) integrate dynamics and get original output
	
	out<-as.data.frame(integrate_dynamics(x0=x0,pars=pars,func_use=func_use,maxtime=tot_time,lengthtime=nsteps,method=method,thresh=THRESH))	# out <- as.data.frame(ode(y = x0, times = seq(0, tot_time, by = step_time), func = dynamics, parms = pars, method = "ode45"))
	# b) produce tidy output
	if(check_zero_sum(out)!="Good"){
		stop("Dynamics were not zero sum")
	}
	
	pcolors <- choose_phenotype_colors(pars$membership, pars$n,pl=pl,pd=pd)
	tidy_out <- out %>% gather(phenotype_n, density, -time) %>% as.tbl() %>% mutate(phenotype_n = as.integer(phenotype_n)) 
	
	# table to use for join 
	# get species names, phenotype number within species, color
	tmp <- data.frame(phenotype_n = as.integer(1:pars$nm), species = as.integer(pars$membership), aaa = 1L)
	tmp <- tmp %>% group_by(species) %>% mutate(phenotype = as.integer(cumsum(aaa))) %>% dplyr::select(-aaa)
	tidy_out <- inner_join(tidy_out, tmp, by = "phenotype_n") %>% mutate(label = (paste0("x[",species,"]^(", phenotype,")")))	%>% 
					mutate(spp_color=scales::alpha(pcolors$sppdarkcolors[match(species,1:length(unique(species)))],alpha))
	## Now do the same for species
	tmp <- data.frame(species = unique(as.integer(pars$membership)))
	tmp$color <- scales::alpha(pcolors$sppdarkcolors,alpha)
	tmp$label <- paste0("x[", tmp$species, "]")
	tidy_species <- tidy_out %>% group_by(time, species) %>% dplyr::summarise(density = sum(density)) %>% inner_join(tmp, by = "species")
	return(list(ode_out = out,
							tidy_out = tidy_out,
							tidy_species = tidy_species,
							pars = pars,
							x0 = x0,
							tot_time = tot_time,
							step_time = step_time))
}

plot_species_dynamics <- function(output,xlab="Generations",ylab="Relative abundance",
								  axis_size=12,label_size=14,lwd=2.5,xbreaks=3,ybreaks=4,
								  x_space=25,y_space=25,right_mar=20){
	df <- output$tidy_species
	pl <- ggplot(df, aes(x = time, y = density, colour = label)) + 
			geom_line(size = lwd) + 
			theme_bw() + 
			scale_color_manual("Species", values = unique(df$color), label = 1:output$pars$n) +#lapply(unique(df$label), parse, file = NULL, n = NULL))
			scale_y_continuous(name=ylab,labels=function(x) sprintf("%.1f", x),breaks=pretty_breaks(ybreaks))+
			scale_x_continuous(name=xlab,labels=function(x) sprintf("%.0f", x),breaks=pretty_breaks(xbreaks))+
			theme(  axis.text=element_text(size=axis_size),
				  	axis.title.x=element_text(size=label_size,margin = unit(c(x_space, 0, 0, 0), "pt")),
				  	axis.title.y=element_text(size=label_size,margin = unit(c(0, y_space, 0, 0), "pt")),
					plot.margin = unit(c(5.5,right_mar,5.5,5.5), "pt"))	
		return(pl)
}



plot_phenotype_dynamics<- function(output,pl,pd,lwd=2){
	
	if(max(output$pars$m)>5) stop("don't have enough line types for more than 5 phenotypes")
	
	line_types<-c("dashed", "dotdash", "longdash", "twodash","dotted")

	#get the ordering of phenotypes within species
	ltp<-NULL
	mem<-1
	val<-0
	for(i in 1:output$pars$nm){
		mem2<-output$par$membership[i]
		if(mem2==mem){
			val<-val+1
			ltp<-c(ltp,val)
		}
		else{
			val<-1
			ltp<-c(ltp,val)
		}
		mem<-mem2
	}

	
	df <- output$tidy_out 
	
	labs<-lapply(unique(df$label), parse, file = NULL, n = NULL)

	labs<-paste0(output$pars$membership,"(",ltp,")")
	
	pl <- ggplot(df, aes(x = time, y = density,colour=interaction(phenotype,species),linetype=interaction(phenotype,species))) + 
		geom_line(size=lwd) + theme_bw() + 
		scale_color_manual("Phenotype", values = rep(unique(df$spp_color),out$pars$m),labels=labs)+
		scale_linetype_manual("Phenotype",values=rep(1,out$pars$nm),labels=labs)+
		guides(colour = guide_legend(keywidth=2.5))
	
	return(pl)
}


plot_phenotype_separate<- function(output,xlab="Generations",trim=0,ylab="Relative abundance",scales="free_y",
								   lwd=2,axis_size=12,label_size=14,ncol=1,pl=pl,pd=pd,xlim=NULL,
								  	xbreaks=3,ybreaks=3, x_space=15,y_space=15,right_mar=20){
	
	nspp<-out$pars$n
	
	if(max(output$pars$m)>5) stop("don't have enough line types for more than 5 phenotypes")
	
	line_types<-c("dashed", "dotdash", "longdash", "twodash","dotted")

	#get the ordering of phenotypes within species
	ltp<-NULL
	mem<-1
	val<-0
	for(i in 1:output$pars$nm){
		mem2<-output$par$membership[i]
		if(mem2==mem){
			val<-val+1
			ltp<-c(ltp,val)
		}
		else{
			val<-1
			ltp<-c(ltp,val)
		}
		mem<-mem2
	}

	
	df <- output$tidy_out 
	
	labs<-lapply(unique(df$label), parse, file = NULL, n = NULL)

	plot_list<-list()
	labs<-paste0("(",output$pars$membership,",",ltp,")")
	pcolors <- choose_phenotype_colors(out$pars$membership, out$pars$n,pl=pl,pd=pd)	
	
	nrow<-ceiling(nspp/ncol)

	grid_vals<-data.frame(species=1:(ncol*nrow),expand.grid(col=1:ncol,row=1:nrow))
	
	df <- left_join(df,grid_vals,"species")
	
	
	if(is.null(xlim)){
		xlim<-c(0,max(df$time))
	}
	gpl<-ggplot(df %>% filter(time>=trim),aes(x = time, y = density,colour=factor(phenotype_n))) + 
		geom_line(size=lwd) + 
		theme_bw() + 
		scale_y_continuous(name=ylab,labels=function(x) sprintf("%.1f", x),breaks=pretty_breaks(ybreaks))+
		scale_x_continuous(name=xlab,labels=function(x) sprintf("%.0f", x),breaks=pretty_breaks(xbreaks),limits = xlim)+
		scale_color_manual(values=pcolors$mycolors)+
		theme(axis.text=element_text(size=axis_size),strip.background = element_blank(), 
			  strip.text.x = element_blank(),strip.text.y = element_blank(),
				axis.title.x=element_text(size=label_size,margin = unit(c(x_space, 0, 0, 0), "pt")),
				axis.title.y=element_text(size=label_size,margin = unit(c(0, y_space, 0, 0), "pt")),
				plot.margin = unit(c(5.5,right_mar,5.5,-20), "pt"))+
		guides(color=FALSE) +facet_grid(row~col,scales=scales)
	


	return(gpl)
}





