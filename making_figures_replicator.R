library(tidyverse)
library(gridExtra)
library(scales)
rm(list=ls())

source("functions_replicator.R")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### plotting params

label_size<-15
main_size<-16
tic_size<-13
offset<- -0.3
offsetv<- 1
main_face<-"bold"
kn<-23

colpal_yellow<-colorRampPalette(c("white","white","#ffffb2","#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8","#253494","#033456","#033456","#033456"))

sim_label<-expression(paste("Phenotypic similarity (",tau,")"))
mem_label<-"Phenotypic memory (p)"
prop_label<-"Diversity"
avg_label<-"Number of phenotypes"
stab_label<-"Return\ntime"

###########################################################
############ stability and return time

res_orig<-read_csv("../../data/simulations/simululations_rho_vs_p_uniform.csv") 
glist<-read.csv("../../data/simulations/n_survive_by1_spec_back.csv",header=T) %>% filter(!is.na(m_avg),rho<1,div>0,!is.na(div),rho<0.95)


 res0<- res_orig %>% mutate(n_prop=div/n,prop_scaled=(n_prop-0.5)/0.5) %>% filter(p<0.951 | p>0.9999,p>0.20)



g1<-ggplot(res0, aes(x = p, y = rho, z = n_prop)) + 
    stat_summary_2d(bins =  c(23,25),color = NA, fun = mean) +		
 	labs(title = "", x = mem_label, y = sim_label, fill = prop_label) +	
    theme_classic()+
	theme(panel.grid = element_blank(), plot.title = element_text(size=main_size, face=main_face,hjust = offset), 
		  axis.title=element_text(size=label_size),axis.text=element_text(size=tic_size),
		  legend.title = element_text(size=label_size-1),legend.text=element_text(size=tic_size-1)) +
	scale_fill_gradientn(colours=colpal_yellow(100),values=seq(0,1,length=100)^(2),limits=c(0.5,0.85),na.value=colpal_yellow(1))
	# scale_fill_gradient2(low="red",mid = "yellow",high="blue",midpoint=7.5)


g2<-ggplot(res0, aes(x = p, y = rho, z = log(-1/resil))) +
    stat_summary_2d(bins = c(23,25), color = NA, fun = mean) + 
 	labs(title = "",  x = mem_label, y = sim_label, fill = stab_label)+ #expression(paste("log(-Re(λ"[1]*")"^-1*")"))) +	
    theme_classic()+
	theme(panel.grid = element_blank(), plot.title = element_text(size=main_size, face=main_face,hjust = offset), 
		  axis.title=element_text(size=label_size),axis.text=element_text(size=tic_size),
		  legend.title = element_text(size=label_size-1),legend.text=element_text(size=tic_size-1)) +
	scale_fill_gradientn(colours=colpal_yellow(100),values=seq(0,1,length=100)^(1.5),limits=c(4.7,11),na.value=colpal_yellow(1))





##################################
# richness as a function of m vs n



set.seed(10)
glist<-glist[sample(1:nrow(glist),500000),]


# colpal_yellow<-colorRampPalette(c("white","#ffffdd","#ffffb2","#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8","#253494","#033456"))


g3<-ggplot(glist %>% mutate(prop_m=m/n,prop_div=div/n) %>% filter(m_avg<=10,m_avg>=1), aes(y = n,x = m_avg, z = prop_div))+xlim(c(2,10))+
    stat_summary_2d(bins=c(25,17),color = NA, fun = mean) + 
 	labs(title = "", y = "Number of species", x = avg_label, fill = prop_label) +
    theme_classic()+
	theme(panel.grid = element_blank(), plot.title = element_text(size=main_size, face=main_face,hjust = offset), 
		  axis.title=element_text(size=label_size),axis.text=element_text(size=tic_size-1),
		  legend.title = element_text(size=label_size-1),legend.text=element_text(size=tic_size)) +
	scale_fill_gradientn(colours=colpal_yellow(100))

g4<-ggplot(glist %>% mutate(mround=round(m_avg),prop_m=m/n,prop_div=div/n,sdv=sqrt(m_var),cv=sdv/m_avg) %>% 
		   	filter(n%in%c(3,9,20),m_avg<=10,m_avg>=2), aes(x = m_avg, y = rho, z = prop_div))+
    stat_summary_2d(bins=c(35,23),color = NA, fun = mean) + xlim(c(2,10))+
 	labs(title = "", x = avg_label, y = sim_label, fill = prop_label) +
    theme_classic()+
	theme(panel.grid = element_blank(), plot.title = element_text(size=main_size, face=main_face,hjust = offset), 
		  axis.title=element_text(size=label_size),axis.text=element_text(size=tic_size-1),
		  legend.title = element_text(size=label_size-1),legend.text=element_text(size=tic_size)) +
	scale_fill_gradientn(colours=colpal_yellow(100),values=seq(0,1,length=100)^(1.5/2))

library(cowplot)

g1<-ggplotGrob(g1)
g2<-ggplotGrob(g2)
g3<-ggplotGrob(g3)
g4<-ggplotGrob(g4)
g4$widths<-g3$widths<-g2$widths<-g1$widths

svg("../../figures/diversity_stability_stacked.svg",width=9,height=6.5)
plot_grid(g3, g4, g1, g2, labels = "AUTO",label_size=16,hjust=0,vjust=1.5)
dev.off()




rm(list = ls())

source("functions_replicator.R")

# Example with general model
light_dyn<-1
dark_dyn<-.6

light_net<-0.6
dark_net<-0.6

# for(i in 1:10000){
# seed<-seed+1
set.seed(681)

scales_val<-"free_y"

pars <- build_community(m = c(1, 2, 3, 4), p = 0.7,rho=0,rho_within=1,range_p=0,range_q=0)
out <- get_phenotype_dynamics(pars, random_seed=7,tot_time =350, step_time =1,alpha=.9,pl=light_dyn,pd=dark_dyn)
g0<-plot_competition_graph_phenotypes(pars$H, pars$membership, pars$n,pl=light_net,pd=dark_net,write_svg=T,height=6.5*.7,width=5.5*.7,file_out="../../figures/network.svg")
g1<-plot_species_dynamics(out,axis_size=30,label_size=35,lwd=2.5,x_space=20,y_space=20,right_mar=15)+theme(legend.position = "none")
g2<-plot_phenotype_separate(out,pl=light_dyn,pd=dark_dyn,axis_size=25,label_size=35,trim=0,scales=scales_val,lwd=2,ylab="",y_space=0,xlim=c(0,300))



svg("../../figures/species_dynamics.svg",width=5.2,height=7*1.1)
g1
dev.off()

svg("../../figures/phenotype_dynamics.svg",width=5.2,height=7*1.1)
g2
dev.off()



######################################
################### H & Q matrices

grid_size<-2
axis_size<-70
# seed<-seed+1
set.seed(10)
pars <- build_community(m = c(2, 3, 1, 4), p = 0.95,rho=0.8,range_p=0.4,range_q=0.2,rho_within=.7)

Q<-pars$Q
Q[1,2]<-0.033
Q[2,1]<-0.016
Q[1,1]<-0.9

Qm<-melt(Q)  %>% mutate(Var1=max(Var1)-Var1,value=ifelse(value==0,NA,log10(value))) %>% 
	mutate(value=(value-min(value,na.rm=T))/max(value-min(value,na.rm=T),na.rm=T)) %>% 
	mutate(value=ifelse(value>.6,value^3,value^(1/2)))

Qm$value[Qm$value==0 & !is.na(Qm$value)]<-min(Qm$value[Qm$value>0 & !is.na(Qm$value)],na.rm=T)

high_col<-"black"
						
high_col<-"#1c5d8c" ##00455e" #002d5b"
low_col<-"#a33030" #red4"

g1<-ggplot(Qm,aes(x = Var2, y = Var1)) +
	geom_tile(aes(fill=value),color=NA) + ylab("Q = ")+xlab("")+
  	scale_fill_gradient2(low="white",high=high_col,na.value = "white",limits=c(0,1)) +
  	theme_bw() + 
	theme(axis.title.y=element_text(size=axis_size,angle=0,vjust=.5,face = "italic"),
		  axis.text=element_blank(),axis.ticks = element_blank(),
		  plot.margin = unit(c(5.5,5.5,50,5.5), "pt"))+
	guides(fill=FALSE)+
	scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) #+
	# geom_vline(xintercept=0.5,size=grid_size)+geom_hline(yintercept= sum(pars$m)-0.5,size=grid_size)

g2<-ggplot(melt(pars$H) %>% mutate(Var1=max(Var1)-Var1),aes(x = Var2, y = Var1)) + 
	geom_tile(aes(fill=value),color=NA) + ylab("H = ")+xlab("")+
  	scale_fill_gradient2(low=low_col,high=high_col,midpoint = 0.5) +
  	theme_bw() + 
	theme(axis.title.y =element_text(size=axis_size,angle=0,vjust=.5,face = "italic"),
		  axis.title.x =element_blank(),		  
		  axis.text=element_blank(),axis.ticks = element_blank(),
		  plot.margin = unit(c(5.5,5.5,50,5.5), "pt"))+
	guides(fill=FALSE)+
	scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+


svg("../../figures/Q_H_stacked.svg",width=5.7,height=9.1)
grid.arrange(g1,g2,nrow=2)
dev.off()




#####################################################
#  three examples
#########################################################


######################################################
### Example 1: rps, one maladapted


set.seed(6)
pars <- build_community(m=c(2,1,1), p = 1)
pars$H<-rbind(c(.5,0,0,0),c(1,.5,0,1),c(1,1,.5,0),c(1,0,1,.5))
p<-0.9
tm<-matrix(1-p,2,2)
diag(tm)<-p
pars$Q[1:2,1:2]<-tm
pars$Q<-pars$Q[c(3,1,2,4),c(3,1,2,4)]
pars$H<-pars$H[c(3,1,2,4),c(3,1,2,4)]
pars$m<-c(1,2,1)
pars$membership<-c(1,2,2,3)
out <- get_phenotype_dynamics(pars, x0=c(.5,1,1,1), tot_time = 4000, step_time = 75, random_seed = 8,nsteps=45,alpha=.9,pl=light_dyn,pd=dark_dyn)

g0<-plot_competition_graph_phenotypes(pars$H, pars$membership, pars$n,pl=light_net,pd=dark_net,write_svg=T,height=3,width=2.5,bottom_pad=2.5,file_out="../../figures/rps_network.svg")
g1<-plot_species_dynamics(out,axis_size=30,label_size=35,lwd=3,x_space=20,y_space=20,right_mar=40,xbreaks=2)+theme(legend.position = "none")
g2<-plot_phenotype_separate(out,pl=light_dyn,pd=dark_dyn,axis_size=25,label_size=35,trim=0,scales=scales_val,lwd=3,ylab="",y_space=0,xbreaks=1)



g1
g2


svg("../../figures/rps_species_dynamics.svg",width=6,height=7*1.1)
g1
dev.off()

svg("../../figures/rps_phenotype_dynamics.svg",width=4,height=7*1.1)
g2
dev.off()


######################################################
# Example 2: two  phenotypes, one maladapted


set.seed(6)
seed<-1
seed<-seed+1
set.seed(3)

pars <- build_community(n=3,mal_case = T, p = 0.75)
out <- get_phenotype_dynamics(pars, x0=c(2,1,2,1,2,2), tot_time = 250000, step_time = 1000, nsteps=100,random_seed = 9,alpha=.9,pl=light_dyn,pd=dark_dyn)


g0<-plot_competition_graph_phenotypes(pars$H, pars$membership, pars$n,pl=light_net,pd=dark_net,write_svg=T,height=3,width=2.5,bottom_pad=2.5,file_out="../../figures/mal_network.svg")
g1<-plot_species_dynamics(out,axis_size=30,label_size=35,lwd=3,x_space=20,y_space=20,right_mar=40,xbreaks=2)+theme(legend.position = "none")
g2<-plot_phenotype_separate(out,pl=light_dyn,pd=dark_dyn,axis_size=25,label_size=35,trim=0,scales=scales_val,lwd=3,ylab="",y_space=0,xbreaks=1)

g1
g2


svg("../../figures/mal_species_dynamics.svg",width=6,height=7*1.1)
g1
dev.off()

svg("../../figures/mal_phenotype_dynamics.svg",width=4,height=7*1.1)
g2
dev.off()




######################################################
# Example 3: two phenotypes, one neutral



set.seed(14)
pars <- build_community(n=3,neut_case = T, p = 0.95)
out <- get_phenotype_dynamics(pars, tot_time = 2500, step_time = 2.5,alpha=.9,pl=light_dyn,pd=dark_dyn)

g0<-plot_competition_graph_phenotypes(pars$H, pars$membership, pars$n,pl=light_net,pd=dark_net,write_svg=T,height=3,width=2.5,bottom_pad=3,file_out="../../figures/neutral_network.svg")
g1<-plot_species_dynamics(out,axis_size=30,label_size=35,lwd=3,x_space=20,y_space=20,right_mar=40,xbreaks=2)+theme(legend.position = "none")
g2<-plot_phenotype_separate(out,pl=light_dyn,pd=dark_dyn,axis_size=25,label_size=35,trim=0,scales=scales_val,lwd=3,ylab="",y_space=0,xbreaks=1)


g1
g2

svg("../../figures/neutral_species_dynamics.svg",width=6,height=7*1.1)
g1
dev.off()

svg("../../figures/neutral_phenotype_dynamics.svg",width=4,height=7*1.1)
g2
dev.off()










######################################################################################
################################ SUPPLEMENTAL FIGURES ################################




###########################
# expected survival plot


set.seed(10)

nsim<-500

nseq<-c(3,5,10,20)
rho_seq<-c(0) #seq(0,1,length=4)
m_mean<-c(1,3,6,12)
m_range<-c(0) #,1,2,3)

simvals<-expand.grid(n=nseq,rho=rho_seq,m_mean=m_mean,m_range=m_range)

kk<-nrow(simvals)

glist<-NULL
for(i in 1:kk){
	
	div<-div_pheno<-NULL
	for(j in 1:nsim){
		# print(paste("scenario",i,"of",kk,";",j,"of",nsim))
		
		#get the community
		m<-assign_phenotypes(simvals$n[i],mean=simvals$m_mean[i],range=simvals$m_range[i])
		membership<-get_membership(m)
		H<-construct_H_matrix(m,rho=simvals$rho[i],return_zipped=F)
		
		# get the number of species and phenotypes that survive
		tdiv<-get_num_survive(membership,H)		
		div<-c(div,tdiv$n)
		div_pheno<-c(div_pheno,tdiv$nm)
		
	}

	# bind the results
	glist<-rbind(glist,compare_obs_expected(m=m,n=simvals$n[i],rho=simvals$rho[i],div=div))
	
}


legend.title<-""
glist_filt<-glist %>% filter(variable!="expected_sp",rho==0,n>3) %>% 
	rename(Richness=richness,Probability=probability) %>% 
	mutate(variable=ifelse(variable=="expected","Expected","Observed")) %>% 
	mutate(n_name=factor(paste("n =",n), levels = paste("n =",sort(unique(n))))) %>% 
	mutate(m_name=factor(paste("m =",m_avg), levels = paste("m =",sort(unique(m_avg)))))
	

gs1<-ggplot(glist_filt,aes(x=Richness,y=Probability,color=variable,shape=variable))+
	geom_point(size=2)+geom_line()+guides(color = guide_legend(legend.title), shape = guide_legend(legend.title))+
	facet_grid(m_name~n_name,scales="free") + theme_bw() +
	theme(strip.text = element_text(size=13),axis.title = element_text(size=13),legend.text = element_text(size=13),axis.text = element_text(size=11))

# ggsave("../../figures/survival_distribution_analytical.svg",plot=g2,device="svg",width=8,height=7)
svg("../../figures/survival_distribution_analytical.svg",width=6,height=5.5)
gs1
dev.off()


##########################################
# distribution of surival probabilities


mseq<-seq(1,20,by=1)
nseq<-c(3,5,10,20)
nm<-data.frame(expand.grid(n=nseq,m=mseq))

nm$val<-0
for(i in 1:nrow(nm)){
	nm$val[i]<-get_expected_value(n=nm$n[i],m=nm$m[i])
}

nm_filt<-nm %>% mutate(n=factor(n,levels=nseq)) %>%  rename(`Number of phenotypes`=m, `Number of species`=n, Richness=val)


gs2<-ggplot(nm_filt,aes(x=`Number of phenotypes`,y=Richness,col=`Number of species`)) +geom_line(size=1.3)+theme_bw()+theme(axis.title = element_text(size=13),legend.title = element_text(size=13),legend.text = element_text(size=12),axis.text = element_text(size=11))

# ggsave("../../figures/exdiversity.svg",plot=g1,device="svg",width=6,height=4)
svg("../../figures/exdiversity.svg",width=6,height=4)
gs2
dev.off()



#######################################
# expected survival as rho increases
label_size<-16

glist<-read.csv("../../data/simulations/n_survive_by1_spec_back.csv",header=T) %>% filter(!is.na(m_avg),rho<1,div>0,!is.na(div),rho<0.95)

set.seed(10)
glist<-glist[sample(1:nrow(glist),500000),]

kmat<-data.frame(expand.grid(nv=1:max(glist$n),mv=1:max(glist$m_avg)))
expect<-NULL
for(i in 1:nrow(kmat)){
	
	val<-get_dist_survive(kmat$nv[i],kmat$mv[i])
	expect<-c(expect,sum(val[,1]*val[,2]))
}
kmat<- kmat %>% mutate(expected=expect)
km2<-dcast(kmat,nv~mv)
kmat$nm<-paste0(kmat$nv,"-",kmat$mv)


# filter and summarize
glist <- glist %>% mutate(nm=paste0(n,"-",round(m_avg))) 
gdat<-left_join(glist,kmat,"nm")%>% mutate(diff_o=div-n/2,prop_o=div/n,diff_e=expected-n/2,prop_e=expected/n) %>% 
	filter(m_avg<10,m_var<3) %>% mutate(cv=round(sqrt(m_var)/m_avg,1)) %>% mutate(rho=round(rho*10,0)/10) %>% 
	mutate(mround=round(m_avg)) %>% group_by(n,rho,mround) %>% 
	dplyr::summarize(mean_cv=mean(cv),mean_div=mean(div),var_m=mean(m_var),expected=mean(expected),diff_o=mean(diff_o),prop_o=mean(prop_o),diff_e=mean(diff_e),prop_e=mean(prop_e)) %>% ungroup() %>%  gather(key="key",value="value","prop_o","prop_e") %>% 
	filter(rho%in%c(0,0.3,0.6,0.9),n%in%c(4,8,16)) %>% mutate(print_tau=factor(paste("τ =",rho)))


# plot

exp_line_size<-1.2
exp_line_col<-"gray35"

g1<-ggplot(gdat %>% dplyr::filter(key!="prop_e") ,aes(x=mround,y=value,color=factor(n),linetype=key))+
	geom_smooth(size=1.5,span=3)+
	xlim(c(2,10))+ylim(c(0.46,1.05))+
	geom_abline(color=exp_line_col,size=exp_line_size,slope=0,intercept=0.5)+
	labs(title = "", x = "Average number of phenotypes per species", y = "Proportion of species", color="Number of\nSpecies") +
    theme_classic()+theme(plot.title = element_text(size=main_size,face="bold",hjust = offset), 
    					  axis.title.x=element_text(size=label_size+1),axis.title.y=element_text(size=label_size+1),
    					  axis.text=element_text(size=tic_size-2))+
	scale_color_manual(values=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f"))+facet_grid(~print_tau) +
	geom_line(aes(x=mround,y=value),data=gdat %>% filter(key=="prop_e"),color=exp_line_col,size=exp_line_size) +
	scale_linetype_manual(guide=FALSE,values=c("solid","dashed"))+
	theme(strip.text.x = element_text(size = 12,face="bold"),legend.position = c(0.15, 0.4),legend.background=element_rect(color="black"),
    		legend.title=element_text(size=12), 
    		legend.text=element_text(size=12))

g1


svg("../../figures/rho_vs_m.svg",width=9,height=3.6)
g1
dev.off()



######################### variable fitness
res_orig<-read_csv("../../data/simulations/simululations_var_fitness_fixedrho.csv")  

res0<- res_orig %>%	mutate(n_prop=div/n,prop_scaled=(n_prop-0.5)/0.5)  

label_size<-16
main_size<-16
tic_size<-13
offset<- -0.3
offsetv<- 1
main_face<-"bold"
kn<-23


colpal_yellow<-colorRampPalette(c("white","white","#ffffb2","#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8","#253494","#033456","#033456","#033456"))

sim_label<-expression(paste("Phenotypic similarity (",tau,")"))
mem_label<-"Phenotypic memory (p)"
prop_label<-"Diversity"
	

g1<-ggplot(res0 %>% mutate(pp=paste("Memory (p) =",p)), aes(y = rho, x = fit_var, z = n_prop)) + 
    stat_summary_2d(bins =  c(18,16), fun = mean) +		
 	labs(title = "", y = sim_label, x = "Variation in demographic rates", fill = prop_label) +	
    theme_classic()+theme(plot.title = element_text(size=main_size, face=main_face,hjust = offset), 
    					  axis.title=element_text(size=label_size), axis.title.x=element_blank(),axis.text=element_text(size=tic_size),
    					  strip.text.x = element_text(size = 13)) +
	scale_fill_gradientn(colours=colpal_yellow(100))+ facet_grid(~pp)


g2<-ggplot(res0 %>% mutate(pp=paste("Memory (p) =",p)), aes(y = rho, x = fit_var, z = (mean_sd)^(1/10))) + 
    stat_summary_2d(bins =  c(18,16), color = NA, fun = mean) + 
 	labs(title = "",  y = sim_label, x = "Variation in demographic rates", fill = "Magnitude of\n oscillations") +	
    theme_classic()+theme(plot.title = element_text(size=main_size,face=main_face,hjust = offset), 
    					  axis.title=element_text(size=label_size),axis.text=element_text(size=tic_size),
    					  strip.text.x = element_text(size = 13)) +
	scale_fill_gradientn(colours=colpal_yellow(100),na.value=colpal_yellow(1))+	facet_grid(~pp)

svg("../../figures/SI_variable_fitness.svg",width=8.5,height=6.5)
grid.arrange(g1,g2)
dev.off()


