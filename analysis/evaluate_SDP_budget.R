require(here)
setwd(paste0(here(), "/../output/"))
budget_scen<-data.frame(site_bud=seq(0,1, length.out=11), spread_bud=seq(1,0,length.out=11))
obj<-data.frame(frac_site=0,frac_spread=0,q_in=0,qbio=0, time1=0,time2=0,time3=0,time4=0, time5=0 )
V_i<-read.csv('../data/streettrees_grid.csv')[,20]
prez<-read.csv('../data/prez_clean_gdk.csv') # invasible host range (from FIA)
qz<-c(0.3,0.6,0.9)
bios<-c(0.1,0.3,0.5)
for (q_in in qz)
{
  
  for (qbio in bios)
  {
    for (scen in c(1,3,5,7,9,11))
    {
      frac_site=budget_scen$site_bud[scen]
      frac_spread=budget_scen$spread_bud[scen]
      d4prime<-read.csv(paste("../../eab_mgmt/output/vecptime",frac_spread,q_in,qbio, "bud.csv", sep="_"))[,2:6]/1000
      obj<-rbind(obj, setNames(c(frac_site,frac_spread,q_in,qbio,colSums(sweep(as.matrix(d4prime),MARGIN=1,as.vector(V_i[prez[,1]]+1),"*"))),names(obj)))
    }
  }
}
obj<-obj[2:nrow(obj),]
library(viridis)
plot((obj$time1)~obj$frac_site, col=viridis(9)[as.factor(paste0(obj$q_in,obj$q_out, obj$qbio))], xlab="Site-focused budget proportion", ylab="Exposed ash street trees", pch="1", ylim=c(5e+03, 2e+05))
points((obj$time2)~obj$frac_site, col=viridis(9)[as.factor(paste0(obj$q_in,obj$q_out, obj$qbio))], xlab="Site-focused budget proportion", ylab="Exposed ash street trees", pch="2")
points((obj$time3)~obj$frac_site, col=viridis(9)[as.factor(paste0(obj$q_in,obj$q_out, obj$qbio))], xlab="Site-focused budget proportion", ylab="Exposed ash street trees", pch="3")
points((obj$time4)~obj$frac_site, col=viridis(9)[as.factor(paste0(obj$q_in,obj$q_out, obj$qbio))], xlab="Site-focused budget proportion", ylab="Exposed ash street trees", pch="4")
points((obj$time5)~obj$frac_site, col=viridis(9)[as.factor(paste0(obj$q_in,obj$q_out, obj$qbio))], xlab="Site-focused budget proportion", ylab="Exposed ash street trees", pch="5")

# 
 mgmt_itme<-read.csv('../../eab_mgmt/output/M_0.9_0.5.csv', header=F) # examine gurobi Management scenario
 #d<-read.csv('../../eab_mgmt/analysis/python/d_0.3_0.3_0.1.csv', header=F) # examine gurobi pest density output
# 
mgmt<-list()
for (time in 6:11)
{
  mgmt[[time]]<-which(mgmt_itme[1800:nrow(mgmt_itme),time-5]==1)
}
qin<-unlist(apply(mgmt_itme[((1799+1):(2*1799)),],2, function(x){length(which(x==1))}))#
qout<-unlist(apply(mgmt_itme[((2*1799+1):(3*1799)),],2, function(x){length(which(x==1))}))
bio<-unlist(apply(mgmt_itme[((3*1799+1):(4*1799)),],2, function(x){length(which(x==1))}))
cost_each<-matrix(c(unlist(lapply(qin,function(x){(x)*(646863/309)})), unlist(lapply(qout,function(x){(x)*(646863/309)})),unlist(lapply(bio,function(x){(x)*(50000)}))),nrow=3,ncol=5,byrow=T)
cost_each<-t(coxxst_each)
colnames(cost_each)<-c("Quar_in", "Quar_out", "Biocontrol")
row.names(cost_each)<-seq(2025, 2045, by=5)
cost_each<-cbind(cost_each,as.matrix((cost_each[,1]+cost_each[,2])/rowSums(cost_each)))
cost_each<-cbind(cost_each,as.matrix(cost_each[,3]/rowSums(cost_each)))
colnames(cost_each)[4:5]<-c("spread_frac", "site_frac")
colMeans(cost_each)
mean(cost_each[,1]/cost_each[,2])
cost_each
