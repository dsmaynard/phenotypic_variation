
rm(list=ls())

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)


source("functions_replicator.R")
# source("phenotypes_parameters.R")
# source("phenotype_dynamics.R")
# source("plot_phenotype_graph.R")




outfile<-paste0("../../data/simulations/simululations_rho_vs_p_uniform.csv")

# conservative threshold to avoid pruning during transiet dynamics
THRESH_sim<-1e-10

# strict thresh once the sims are done
THRESH_prune<-1e-6

# number of species
nseq<-1 #c(5,10,15)

rhoseq<-seq(0,0.95,length=25)
# mmean<-c(2,3,5,10)
mprop<-1 #c(2,3,5) #seq(0.4,1,length=4)
pseq<-0.999999 #seq(0.95,0.99,length=25) #seq(0.05,0.19,by=0.02916667)
range_m<-1 #c(0,1,3,5)

range_q<-0.1
range_p<-0.1

nsim<-2000
max_count<-80
max_count2<-30
	
#i<-695
maxed<-maxed2<-0

#svals<-expand.grid(n=nseq,p=pseq,mprop=mprop,rho=rhoseq,range_m=range_m,range_p=range_p,range_q=range_q,nsim=1)
svals<-expand.grid(p=pseq,rho=rhoseq,nsim=1)

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
		
		n<-sample(5:15,1)
		m<-sample(1:5,n,replace=T)
		

		# get phenotypes for each
		# m<-assign_phenotypes(svals$n[i],mean=svals$mprop[i],range=svals$range_m[i]) #svals$mean[i]+round(svals$mean[i]*svals$range_prop[i]))
		membership<-get_membership(m)
				

		# construct the matrices
		H<-H0<-construct_H_matrix(m,rho=svals$rho[i])
		Q<-Q0<-construct_Q_matrix(m,p=svals$p[i],range_p = 0.15,range_q = 0.15,ignore_min = F)
		P<-H-t(H)
		B<-make_species_sum_matrix(m)			
	
		# get the starting conditions
		xst<-find_optimal_strategy(H,verbose=F)+0.01
		xst<-xst/sum(xst)

		# do the initial integration to get past transients
		resm<-integrate_dynamics(x0=xst,pars=list(H=H,Q=Q),maxtime=500,lengthtime=500,maxsteps=1e6,method="lsoda",thresh=THRESH_sim,func_use=zero_sum_mutator)
	
		#now integrate for a bit to get oscillations
		resm<-integrate_dynamics(x0=resm[nrow(resm),-1],pars=list(H=H,Q=Q),maxtime=5000,lengthtime=1000,maxsteps=1e6,method="lsoda",thresh=THRESH_sim,func_use=zero_sum_mutator)		
		soln<-apply(resm,2,mean)[-1]
	
		# continually integrate and average until we find the solution			
		count<-0
		skip<-FALSE
		while(any(point_dxdt(soln,H,Q)>1e-7)){
			if(count==max_count){
				break()
			}		
			count<-count+1
			resm<-integrate_dynamics(x0=soln,pars=list(H=H,Q=Q),maxtime=1000,lengthtime=500,maxsteps=1e6,method="lsoda",thresh=THRESH_sim,func_use=zero_sum_mutator)
			soln<-apply(resm,2,mean)[-1]			
		}
			
		# plotting, for debugging
		# update<-prune_system(resm=resm,m=m,Q=Q,H=H,membership = membership,THRESH=0,prop=1,show_plot=T)

		# make sure we aren't pruning partial species
		alive<-soln>THRESH_prune
		alive_tab<-table(membership[alive])

		count2<-0
		if(count<max_count){
			while(!all(m[as.numeric(names(alive_tab))]==alive_tab)){
				if(count2==max_count2){
					# can't get it to converge
					break()
				}
				resm<-integrate_dynamics(x0=soln,pars=list(H=H,Q=Q),maxtime=1000,lengthtime=500,maxsteps=1e6,method="lsoda",thresh=THRESH_sim,func_use=zero_sum_mutator)
				soln<-apply(resm,2,mean)[-1]
				alive<-soln>THRESH_prune
				alive_tab<-table(membership[alive])	
				count2<-count2+1
			}
		}
		


		
		check<-FALSE
		div<-resil<-NA		

		if(count<max_count & count2<max_count2){
			
			# prune the system and get the eigenvalue
			H1<-H[alive,alive]
			Q1<-Q[alive,alive]
			soln1<-soln[alive]
			soln1<-soln1/sum(soln1)
			resm1<-integrate_dynamics(x0=soln1,pars=list(H=as.matrix(H1),Q=as.matrix(Q1)),maxtime=100,lengthtime=100,maxsteps=1e6,method="lsoda",thresh=THRESH_prune,func_use=zero_sum_mutator)
			soln1<-as.numeric(resm1[nrow(resm1),-1])
			ev1<-Re(eigen(Jacobian_H(soln1,as.matrix(Q1),as.matrix(H1)))$values[-1])	

			# calculate the metrics
			if(sum(alive)==1){
				div<-1
				resil<-NA
			}
			else if(all(ev1<0)){
				resil<-max(ev1)
				soln[soln<THRESH_prune]<-0
				div<-sum((B%*%soln)>THRESH_prune)
				# if(abs(resil)>0.01 & div>1){
				# 	stop("found one")
				# }
			}
			else{
				# stop("have positive eigenvalues")		
				div<-resil<-NA	
				check<-TRUE
			}
		}
		else if(count<max_count & count2==max_count2){
			maxed2<-maxed2+1
			# stop()
			print(c(maxed,maxed2))
		}
		else if(count==max_count & count2<max_count2){
			maxed<-maxed+1
			print(c(maxed,maxed2))
		}
		
		# print every 200
		res1<-rbind(res1,data.frame(n=n,div=div,resil=resil,seed=randv,check=check,count=count,count2=count2,
									mean_m=mean(m),min_m=min(m),median_m=median(m),max_m=max(m),var_m=var(m),
									mean_p=mean(diag(Q0)),min_p=min(diag(Q0)),median_p=median(diag(Q0)),max_p=max(diag(Q0)),var_p=var(diag(Q0))))
		
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
