rm(list=ls())
require(here)
setwd("../output/")
prez<-read.csv('../data/prez_clean_gdk.csv') # invasible host range (from FIA)
prez2<-read.csv('../data/prez2_clean_gdk.csv') # pest presences for all species
prez2[,1]<-readRDS('presences_time_eab.rds')[[1]][,5] # pest presences in 2020
L<-rep(0,64) # size of each pest's host range
V_i<-read.csv('../../eab_mgmt/data/streettrees_grid.csv')[,20]
L<-length(unique(c(which(V_i>0), prez[,1], prez2[,1])))-1
prez[,1]<-c(unique(c(which(V_i!=0), prez[which(prez[,1]!=0),1], prez2[which(prez2[,1]!=0),1])), rep(0, 3372-L))

# 
mgmt_itme<-read.csv('M_0.9_0.5.csv') # examine gurobi Management scenario, change file name to different efficiency scenaris

mgmt<-list()
for (time in 6:11)
{
  mgmt[[time]]<-which(mgmt_itme[(L+1):nrow(mgmt_itme),time-5]==1)
}
qin<-unlist(apply(mgmt_itme[((L+1):(2*L)),],2, function(x){length(which(x==1))}))#
qout<-unlist(apply(mgmt_itme[((2*L+1):(3*L)),],2, function(x){length(which(x==1))}))
bio<-unlist(apply(mgmt_itme[((3*L+1):(4*L)),],2, function(x){length(which(x==1))}))
cost_each<-matrix(c(unlist(lapply(qin,function(x){(x)*(15421)})), unlist(lapply(qout,function(x){(x)*(15421)})),unlist(lapply(bio,function(x){(x)*(50000)}))),nrow=3,ncol=5,byrow=T)
cost_each<-t(cost_each)
colnames(cost_each)<-c("Quar_in", "Quar_out", "Biocontrol")
row.names(cost_each)<-seq(2025, 2045, by=5)
cost_each<-cbind(cost_each,as.matrix((cost_each[,1]+cost_each[,2])/rowSums(cost_each)))
cost_each<-cbind(cost_each,as.matrix(cost_each[,3]/rowSums(cost_each)))
colnames(cost_each)[4:5]<-c("spread_frac", "site_frac")
colMeans(cost_each)
mean(cost_each[,1]/cost_each[,2])
cost_each
colMeans(cost_each[,1:3]/rowSums(cost_each[,1:3]))
sum(cost_each[,1:3])/(1650000*5)