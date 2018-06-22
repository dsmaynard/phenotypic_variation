
rm(list=ls())

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
library(numDeriv)

source("functions_replicator.R")
# source("phenotypes_parameters.R")
# source("phenotype_dynamics.R")
# source("plot_phenotype_graph.R")




outfile<-paste0("../../data/simulations/simululations_var_fitness_fixedrho.csv")

# conservative threshold to avoid pruning during transiet dynamics
THRESH_sim<-1e-10

# strict thresh once the sims are done
THRESH_prune<-1e-6

# number of species
n<-10 #c(5,10,15)
m<-rep(5,n)


rhoseq<-seq(0,1,length=15)
pseq<-c(0.5,0.75,0.95) #seq(0.2,0.95,length=10) #seq(0.05,0.19,by=0.02916667)


range_q<-0
range_p<-0

nsim<-50000
max_count<-80
max_count2<-30
max_count_null<-5

fit_var<-seq(0.001,0.2,length=20)
#i<-695
maxed<-maxed2<-0

#svals<-expand.grid(n=nseq,p=pseq,mprop=mprop,rho=rhoseq,range_m=range_m,range_p=range_p,range_q=range_q,nsim=1)
svals<-expand.grid(p=pseq,rho=rhoseq,fit_var=fit_var)

# shuffle the ordering
svals<-svals[sample(1:nrow(svals),nrow(svals),replace=F),]
print(kn<-nrow(svals))

trunc<-1000

for(j in 1:nsim){
	
	res1<-NULL
	
	for(i in 1:nrow(svals)){
	
		# for(j in 1:1000){			
	
		print(paste0("simulation ",i," of ",kn,"; rep ",j," of ",nsim))
		randv<-sample.int(1e9,1)
		set.seed(randv)
		
		# get phenotypes for each
		# m<-assign_phenotypes(svals$n[i],mean=svals$mprop[i],range=svals$range_m[i]) #svals$mean[i]+round(svals$mean[i]*svals$range_prop[i]))
		membership<-get_membership(m)
				

		# construct the matrices
		H<-H0<-construct_H_matrix(m,rho=svals$rho[i])
		Q<-Q0<-construct_Q_matrix(m,p=svals$p[i],range_p = 0,range_q = 0,ignore_min = F)
		P<-H-t(H)
		B<-make_species_sum_matrix(m)			
	
		# get the starting conditions
		xst<-find_optimal_strategy(H,verbose=F)+0.01
		xst<-xst/sum(xst)
		
		d<-runif(sum(m),1-svals$fit_var[i],1+svals$fit_var[i])
		f<-runif(sum(m),1-svals$fit_var[i],1+svals$fit_var[i])
		
		# do the initial integration to get past transients
		resm<-integrate_dynamics(x0=xst,pars=list(H=H,Q=Q,f=f,d=d),maxtime=1000,lengthtime=1000,maxsteps=1e6,method="lsoda",thresh=THRESH_sim,func_use=zero_sum_mutator_fitness)
		soln<-resm[nrow(resm),-1]
		soln[soln<THRESH_prune]<-0
		#now integrate for a bit to get oscillations
		resm<-integrate_dynamics(x0=soln,pars=list(H=H,Q=Q,f=f,d=d),maxtime=20000,lengthtime=500,maxsteps=1e6,method="lsoda",thresh=THRESH_sim,func_use=zero_sum_mutator_fitness)		
		soln<-resm[nrow(resm),-1]
		soln[soln<THRESH_prune]<-0
		resm<-integrate_dynamics(x0=soln,pars=list(H=H,Q=Q,f=f,d=d),maxtime=1000,lengthtime=1000,maxsteps=1e6,method="lsoda",thresh=THRESH_sim,func_use=zero_sum_mutator_fitness)		
		sdval<-apply(resm[,-1],2,sd)
		meanval<-apply(resm[,-1],2,mean)
		cvval<-sdval/meanval
		
		# update<-prune_system(resm=resm,m=m,Q=Q,H=H,membership = membership,THRESH=0,prop=1,show_plot=T)
		
		alive<-soln>THRESH_prune
		alive_tab<-table(membership[alive])
		
		div<-length(alive_tab)
	
			
		# print every 200
		res1<-rbind(res1,data.frame(n=n,div=div,mean_cv=mean(cvval,na.rm=T),mean_sd=mean(sdval),x_est=mean(meanval),seed=randv))
		
	}
	
	res1<-data.frame(res1,svals)
	# if it exists, just append
	if(file.exists(outfile)){
		write_csv(res1,outfile,append=T)
	}
	else{# otherwise, create
		write_csv(res1,outfile)
	}

}
