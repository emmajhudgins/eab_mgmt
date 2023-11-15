require(here)
library(ggplot2)
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
      d4prime<-read.csv(paste("vecptime",frac_spread,q_in,qbio, "bud.csv", sep="_"))[,2:6]/1000
      obj<-rbind(obj, setNames(c(frac_site,frac_spread,q_in,qbio,colSums(sweep(as.matrix(d4prime),MARGIN=1,as.vector(V_i[prez[,1]]+1),"*"))),names(obj)))
    }
  }
}
obj<-obj[2:nrow(obj),]
library(viridis)
plot((obj$time1)~obj$frac_site, col=viridis(9)[as.factor(paste0(obj$q_in,obj$q_out, obj$qbio))], xlab="Site-focused budget proportion", ylab="Exposed ash street trees", pch="1", ylim=c(0, 2e+05))
points((obj$time2)~obj$frac_site, col=viridis(9)[as.factor(paste0(obj$q_in,obj$q_out, obj$qbio))], xlab="Site-focused budget proportion", ylab="Exposed ash street trees", pch="2")
points((obj$time3)~obj$frac_site, col=viridis(9)[as.factor(paste0(obj$q_in,obj$q_out, obj$qbio))], xlab="Site-focused budget proportion", ylab="Exposed ash street trees", pch="3")
points((obj$time4)~obj$frac_site, col=viridis(9)[as.factor(paste0(obj$q_in,obj$q_out, obj$qbio))], xlab="Site-focused budget proportion", ylab="Exposed ash street trees", pch="4")
points((obj$time5)~obj$frac_site, col=viridis(9)[as.factor(paste0(obj$q_in,obj$q_out, obj$qbio))], xlab="Site-focused budget proportion", ylab="Exposed ash street trees", pch="5")


obj<-data.frame(frac_site=0,frac_spread=0,q_in=0,q_out=0,qbio=0, obj=0)
names(obj)

for (q_in in qz)
{
  for (q_out in qz)
  {
    for (qbio in bios)
    {
      for (scen in c(1,3,5,7,9,11))
      {
        frac_site=budget_scen$site_bud[scen]
        frac_spread=budget_scen$spread_bud[scen]
        d4prime<-read.csv(paste("vecptime",frac_spread,q_in,qbio, "bud.csv", sep="_"), header=F)[,2:8]/1000
        obj<-rbind(obj, setNames(c(frac_site,frac_spread,q_in,q_out,qbio,sum(sweep(as.matrix(d4prime),MARGIN=1,as.vector(V_i[prez[,1]]+1),"*"))),names(obj)))
      }
    }
  }
}
obj<-obj[2:nrow(obj),]
plot(y=obj$obj,x=(obj$frac_site), col=alpha(viridis(9)[as.factor(paste0(obj$q_in, obj$qbio))],0.75), xlab="Site-focused budget proportion", ylab="Exposed ash street trees", pch=19, xlim=c(-0.05,1.05), ylim=c(100000,1400000))
library(readxl)
dat<-read.csv('results_table.csv')
points(y=dat$Exposed.Ash.Street.Trees..thousands.*1000, x=dat$bio_prop/100, lty=2,col=viridis(9), pch=5)

mgmt_itme<-read.csv('M_0.6_0.1.csv', header=F) # examine gurobi Management scenario
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
cost_each<-t(cost_each)
colnames(cost_each)<-c("Quar_in", "Quar_out", "Biocontrol")
row.names(cost_each)<-seq(2025, 2045, by=5)
cost_each<-cbind(cost_each,as.matrix((cost_each[,1]+cost_each[,2])/rowSums(cost_each)))
cost_each<-cbind(cost_each,as.matrix(cost_each[,3]/rowSums(cost_each)))
colnames(cost_each)[4:5]<-c("spread_frac", "site_frac")
colMeans(cost_each)
mean(cost_each[,1]/(cost_each[,1]+cost_each[,2]))
cost_each

colMeans(cost_each/rowSums(cost_each[,1:3]))

sum(cost_each[,1:3])/(5*1650000)
obj<-data.frame(q_in=0,q_out=0,qbio=0, obj=0)

for (q_in in qz)
{
  q_out=q_in
  for (qbio in bios)
  {
    d4prime<-read.csv(paste("vecptime",q_in,qbio, "fg.csv", sep="_"), header=F)[,2:8]/1000
    obj<-rbind(obj, setNames(c(q_in,q_out,qbio,sum(sweep(as.matrix(d4prime),MARGIN=1,as.vector(V_i[prez[,1]]+1),"*"))),names(obj)))
  }
}
obj<-obj[2:nrow(obj),]

layout(matrix(c(1,2), ncol=2, byrow=TRUE), heights=c(1,1), widths=c(0.7,0.3))
par(mar=c(4,4,2,0))
par(oma=c(0,0,0,0))
par(xpd = FALSE)


plot(y=obj$obj,x=(obj$q_in), col=alpha(viridis(9)[as.factor(paste0(obj$q_in, obj$qbio))],0.75),xlab="Quarantine efficiency", ylab="Exposed ash street trees", pch=19, xlim=c(0,0.95), ylim=c(0,900000), main="Fast Growth Scenario", cex.main=0.75, axes=F)
axis(1,labels=c("0","30%", "60%", '90%'), at=c(0,0.3,0.6,0.9), outer=F)
axis(2, labels=c('0',"200000", "400000", '600000', '800000', '1000000'), at=c(0,200000,400000,600000,800000, 1000000), outer=F, las=2, cex.axis=0.75, hadj=0.8)
dat<-read.csv('~/Downloads/postdocdat.csv')
points(y=dat$Exposed.Ash.Street.Trees..thousands.*1000, x=dat$Quarantine.efficiency+0.05, lty=2,col=viridis(9), pch=5)
par(mai=c(0.4,0.5,0.1,0.6))

image(1,1:9,t(matrix(1:9)), col=viridis(9), axes=FALSE, ann=F)
axis(2,labels=c("30%","30%","30%","60%","60%","60%", "90%","90%","90%"),at=c(1:9), cex.axis=0.75, padj=1)
axis(4,at=c(1:9), labels=rep("", 9),cex.axis=0.75, padj=2)
par(xpd = TRUE) #Draw outside plot area
corners<-par("usr")
text("50%         30%        10%",y=c(8.1), x=corners[2]+.2, cex=0.75, srt=270)
text("50%         30%        10%",y=c(5.1), x=corners[2]+.2, cex=0.75, srt=270)
text("50%         30%        10%",y=c(2.1), x=corners[2]+.2, cex=0.75, srt=270)

mtext(side=2, "Quarantine efficiency",line=1.25, cex=0.75)
corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
text(x = corners[2]+0.4, y = mean(corners[3:4]),"Biological Control efficiency",cex=0.7, srt = 270)

obj<-data.frame(q_in=0,q_out=0,qbio=0, obj=0)

for (q_in in qz)
{
  q_out=q_in
  for (qbio in bios)
  {
    d4prime<-read.csv(paste("vecptime",q_in,qbio, "fdfg.csv", sep="_"), header=F)[,2:8]/1000
    obj<-rbind(obj, setNames(c(q_in,q_out,qbio,sum(sweep(as.matrix(d4prime),MARGIN=1,as.vector(V_i[prez[,1]]+1),"*"))),names(obj)))
  }
}
obj<-obj[2:nrow(obj),]

layout(matrix(c(1,2), ncol=2, byrow=TRUE), heights=c(1,1), widths=c(0.7,0.3))
par(mar=c(4,4,2,0))
par(oma=c(0,0,0,0))
par(xpd = FALSE)


plot(y=obj$obj,x=(obj$q_in), col=alpha(viridis(9)[as.factor(paste0(obj$q_in, obj$qbio))],0.75),xlab="Quarantine efficiency", ylab="Exposed ash street trees", pch=19, xlim=c(0,0.95), ylim=c(0,900000), main="Fast Dispersal and Growth Scenario", cex.main=0.75, axes=F)
axis(1,labels=c("0","30%", "60%", '90%'), at=c(0,0.3,0.6,0.9), outer=F)
axis(2, labels=c('0',"200000", "400000", '600000', '800000', '1000000'), at=c(0,200000,400000,600000,800000, 1000000), outer=F, las=2, cex.axis=0.75, hadj=0.8)
dat<-read.csv('~/Downloads/postdocdat.csv')
points(y=dat$Exposed.Ash.Street.Trees..thousands.*1000, x=dat$Quarantine.efficiency+0.05, lty=2,col=viridis(9), pch=5)
par(mai=c(0.4,0.5,0.1,0.6))


image(1,1:9,t(matrix(1:9)), col=viridis(9), axes=FALSE, ann=F)
axis(2,labels=c("30%","30%","30%","60%","60%","60%", "90%","90%","90%"),at=c(1:9), cex.axis=0.75, padj=1)
axis(4,at=c(1:9), labels=rep("", 9),cex.axis=0.75, padj=2)
par(xpd = TRUE) #Draw outside plot area
corners<-par("usr")
text("50%         30%        10%",y=c(8.1), x=corners[2]+.2, cex=0.75, srt=270)
text("50%         30%        10%",y=c(5.1), x=corners[2]+.2, cex=0.75, srt=270)
text("50%         30%        10%",y=c(2.1), x=corners[2]+.2, cex=0.75, srt=270)

mtext(side=2, "Quarantine efficiency",line=1.25, cex=0.75)
corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
text(x = corners[2]+0.4, y = mean(corners[3:4]),"Biological Control efficiency",cex=0.7, srt = 270)

obj<-data.frame(q_in=0,q_out=0,qbio=0, obj=0)

for (q_in in qz)
{
  q_out=q_in
  for (qbio in bios)
  {
    d4prime<-read.csv(paste("vecptime",q_in,qbio, "nothresh.csv", sep="_"), header=F)[,2:8]/1000
    obj<-rbind(obj, setNames(c(q_in,q_out,qbio,sum(sweep(as.matrix(d4prime),MARGIN=1,as.vector(V_i[prez[,1]]+1),"*"))),names(obj)))
  }
}
obj<-obj[2:nrow(obj),]

layout(matrix(c(1,2), ncol=2, byrow=TRUE), heights=c(1,1), widths=c(0.7,0.3))
par(mar=c(4,4,2,0))
par(oma=c(0,0,0,0))
par(xpd = FALSE)


plot(y=obj$obj,x=(obj$q_in), col=alpha(viridis(9)[as.factor(paste0(obj$q_in, obj$qbio))],0.75),xlab="Quarantine efficiency", ylab="Exposed ash street trees", pch=19, xlim=c(0,0.95), ylim=c(0,900000), main="No Threshold Scenario", cex.main=0.75, axes=F)
axis(1,labels=c("0","30%", "60%", '90%'), at=c(0,0.3,0.6,0.9), outer=F)
axis(2, labels=c('0',"200000", "400000", '600000', '800000', '1000000'), at=c(0,200000,400000,600000,800000, 1000000), outer=F, las=2, cex.axis=0.75, hadj=0.8)
dat<-read.csv('~/Downloads/postdocdat.csv')
points(y=dat$Exposed.Ash.Street.Trees..thousands.*1000, x=dat$Quarantine.efficiency+0.05, lty=2,col=viridis(9), pch=5)
par(mai=c(0.4,0.5,0.1,0.6))


image(1,1:9,t(matrix(1:9)), col=viridis(9), axes=FALSE, ann=F)
axis(2,labels=c("30%","30%","30%","60%","60%","60%", "90%","90%","90%"),at=c(1:9), cex.axis=0.75, padj=1)
axis(4,at=c(1:9), labels=rep("", 9),cex.axis=0.75, padj=2)
par(xpd = TRUE) #Draw outside plot area
corners<-par("usr")
text("50%         30%        10%",y=c(8.1), x=corners[2]+.2, cex=0.75, srt=270)
text("50%         30%        10%",y=c(5.1), x=corners[2]+.2, cex=0.75, srt=270)
text("50%         30%        10%",y=c(2.1), x=corners[2]+.2, cex=0.75, srt=270)

mtext(side=2, "Quarantine efficiency",line=1.25, cex=0.75)
corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
text(x = corners[2]+0.4, y = mean(corners[3:4]),"Biological Control efficiency",cex=0.7, srt = 270)
