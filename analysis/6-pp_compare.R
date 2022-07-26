#script to compare propagule pressure in each site across management actions. This creates Fig. S3 and the ANOVA+Tukey test results.

#Written by Emma J. Hudgins, emma.hudgins@carleton.ca

rm(list=ls()) 
library(dplyr)
library(here)
library(ggplot2)

setwd(paste0(here(), '/../output/'))

mgmt_sites<-tibble(action=character(), pp=numeric(), qin=numeric(), qbio=numeric())

for (q_in in c(0.3,0.6,0.9)){
  for (qbio in c(0.1,0.3,0.5)){
    vecP_time<-read.csv(paste("vecptime_",q_in,"_",qbio, ".csv", sep=""), header=F)[,1:7]/1000
    vecP_time<-cbind(rep(0,1799),rep(0,1799),rep(0,1799),rep(0,1799),rep(0,1799),vecP_time)
    mgmt_itme<-read.csv(paste("M_",q_in,"_",qbio, ".csv", sep=""), header=F)
    mgmt_sites<-rbind(mgmt_sites,setNames(data.frame(cbind(rep("bio", length(c(vecP_time[which(mgmt_itme[(3*1799+1):(4*1799),1]==1),6], vecP_time[which(mgmt_itme[(3*1799+1):(4*1799),2]==1),7], vecP_time[which(mgmt_itme[(3*1799+1):(4*1799),3]==1),8], vecP_time[which(mgmt_itme[(3*1799+1):(4*1799),4]==1),9], vecP_time[which(mgmt_itme[(3*1799+1):(4*1799),5]==1),10]))), c(vecP_time[which(mgmt_itme[(3*1799+1):(4*1799),1]==1),6], vecP_time[which(mgmt_itme[(3*1799+1):(4*1799),2]==1),7], vecP_time[which(mgmt_itme[(3*1799+1):(4*1799),3]==1),8], vecP_time[which(mgmt_itme[(3*1799+1):(4*1799),4]==1),9], vecP_time[which(mgmt_itme[(3*1799+1):(4*1799),5]==1),10]), rep(q_in, length(c(vecP_time[which(mgmt_itme[(3*1799+1):(4*1799),1]==1),6], vecP_time[which(mgmt_itme[(3*1799+1):(4*1799),2]==1),7], vecP_time[which(mgmt_itme[(3*1799+1):(4*1799),3]==1),8], vecP_time[which(mgmt_itme[(3*1799+1):(4*1799),4]==1),9], vecP_time[which(mgmt_itme[(3*1799+1):(4*1799),5]==1),10]))), rep(qbio, length(c(vecP_time[which(mgmt_itme[(3*1799+1):(4*1799),1]==1),6], vecP_time[which(mgmt_itme[(3*1799+1):(4*1799),2]==1),7], vecP_time[which(mgmt_itme[(3*1799+1):(4*1799),3]==1),8], vecP_time[which(mgmt_itme[(3*1799+1):(4*1799),4]==1),9], vecP_time[which(mgmt_itme[(3*1799+1):(4*1799),5]==1),10]))))) , names(mgmt_sites)))
    mgmt_sites<-bind_rows(mgmt_sites, setNames(data.frame(cbind(rep("qin", length(c(vecP_time[which(mgmt_itme[(1*1799+1):(2*1799),1]==1),6], vecP_time[which(mgmt_itme[(1*1799+1):(2*1799),2]==1),7], vecP_time[which(mgmt_itme[(1*1799+1):(2*1799),3]==1),8], vecP_time[which(mgmt_itme[(1*1799+1):(2*1799),4]==1),9], vecP_time[which(mgmt_itme[(1*1799+1):(2*1799),5]==1),10]))),c(vecP_time[which(mgmt_itme[(1*1799+1):(2*1799),1]==1),6], vecP_time[which(mgmt_itme[(1*1799+1):(2*1799),2]==1),7], vecP_time[which(mgmt_itme[(1*1799+1):(2*1799),3]==1),8], vecP_time[which(mgmt_itme[(1*1799+1):(2*1799),4]==1),9], vecP_time[which(mgmt_itme[(1*1799+1):(2*1799),5]==1),10]),rep(q_in, length(c(vecP_time[which(mgmt_itme[(1*1799+1):(2*1799),1]==1),6], vecP_time[which(mgmt_itme[(1*1799+1):(2*1799),2]==1),7], vecP_time[which(mgmt_itme[(1*1799+1):(2*1799),3]==1),8], vecP_time[which(mgmt_itme[(1*1799+1):(2*1799),4]==1),9], vecP_time[which(mgmt_itme[(1*1799+1):(2*1799),5]==1),10]))),rep(qbio, length(c(vecP_time[which(mgmt_itme[(1*1799+1):(2*1799),1]==1),6], vecP_time[which(mgmt_itme[(1*1799+1):(2*1799),2]==1),7], vecP_time[which(mgmt_itme[(1*1799+1):(2*1799),3]==1),8], vecP_time[which(mgmt_itme[(1*1799+1):(2*1799),4]==1),9], vecP_time[which(mgmt_itme[(1*1799+1):(2*1799),5]==1),10]))))), names(mgmt_sites)))
    mgmt_sites<-bind_rows(mgmt_sites, setNames(data.frame(cbind(rep("qout", length(c(vecP_time[which(mgmt_itme[(2*1799+1):(3*1799),1]==1),6], vecP_time[which(mgmt_itme[(2*1799+1):(3*1799),2]==1),7], vecP_time[which(mgmt_itme[(2*1799+1):(3*1799),3]==1),8], vecP_time[which(mgmt_itme[(2*1799+1):(3*1799),4]==1),9], vecP_time[which(mgmt_itme[(2*1799+1):(3*1799),5]==1),10]))),c(vecP_time[which(mgmt_itme[(2*1799+1):(3*1799),1]==1),6], vecP_time[which(mgmt_itme[(2*1799+1):(3*1799),2]==1),7], vecP_time[which(mgmt_itme[(2*1799+1):(3*1799),3]==1),8], vecP_time[which(mgmt_itme[(2*1799+1):(3*1799),4]==1),9], vecP_time[which(mgmt_itme[(2*1799+1):(3*1799),5]==1),10]), rep(q_in, length(c(vecP_time[which(mgmt_itme[(2*1799+1):(3*1799),1]==1),6], vecP_time[which(mgmt_itme[(2*1799+1):(3*1799),2]==1),7], vecP_time[which(mgmt_itme[(2*1799+1):(3*1799),3]==1),8], vecP_time[which(mgmt_itme[(2*1799+1):(3*1799),4]==1),9], vecP_time[which(mgmt_itme[(2*1799+1):(3*1799),5]==1),10]))), rep(qbio, length(c(vecP_time[which(mgmt_itme[(2*1799+1):(3*1799),1]==1),6], vecP_time[which(mgmt_itme[(2*1799+1):(3*1799),2]==1),7], vecP_time[which(mgmt_itme[(2*1799+1):(3*1799),3]==1),8], vecP_time[which(mgmt_itme[(2*1799+1):(3*1799),4]==1),9], vecP_time[which(mgmt_itme[(2*1799+1):(3*1799),5]==1),10]))))), names(mgmt_sites)))
    mgmt_sites<-bind_rows(mgmt_sites, setNames(data.frame(cbind(rep("none", length(c(vecP_time[which(mgmt_itme[(1):(1799),1]==1),6], vecP_time[which(mgmt_itme[(1):(1799),2]==1),7], vecP_time[which(mgmt_itme[(1):(1799),3]==1),8], vecP_time[which(mgmt_itme[(1):(1799),4]==1),9], vecP_time[which(mgmt_itme[(1):(1799),5]==1),10]))),c(vecP_time[which(mgmt_itme[(1):(1799),1]==1),6], vecP_time[which(mgmt_itme[(1):(1799),2]==1),7], vecP_time[which(mgmt_itme[(1):(1799),3]==1),8], vecP_time[which(mgmt_itme[(1):(1799),4]==1),9], vecP_time[which(mgmt_itme[(1):(1799),5]==1),10]),rep(q_in, length(c(vecP_time[which(mgmt_itme[(1):(1799),1]==1),6], vecP_time[which(mgmt_itme[(1):(1799),2]==1),7], vecP_time[which(mgmt_itme[(1):(1799),3]==1),8], vecP_time[which(mgmt_itme[(1):(1799),4]==1),9], vecP_time[which(mgmt_itme[(1):(1799),5]==1),10]))), rep(qbio, length(c(vecP_time[which(mgmt_itme[(1):(1799),1]==1),6], vecP_time[which(mgmt_itme[(1):(1799),2]==1),7], vecP_time[which(mgmt_itme[(1):(1799),3]==1),8], vecP_time[which(mgmt_itme[(1):(1799),4]==1),9], vecP_time[which(mgmt_itme[(1):(1799),5]==1),10]))))), names(mgmt_sites)))
  }}
mgmt_sites$qin<-as.numeric(mgmt_sites$qin)
mgmt_sites$qbio<-as.numeric(mgmt_sites$qbio)
mgmt_sites$action<-as.factor(mgmt_sites$action)
mgmt_sites$pp<-as.numeric(mgmt_sites$pp)
m<-lm(log(pp+1)~action+qin+qbio, mgmt_sites)
library(multcomp)
mult<-glht(m, linfct=mcp(action="Tukey"))
summary(mult)
summary<-mgmt_sites%>%group_by(action, qin)%>%summarize(mean(pp), sqrt(var(pp)/length(pp)))
write.csv(summary, row.names=F, file="summarybyaction.csv")

mgmt_sites$eff<-as.factor(mgmt_sites$qin)

plot<-ggplot(aes(x=action, y=pp, fill=eff), data=mgmt_sites)+geom_boxplot(aes(y=pp,x=action, fill=eff))+scale_y_log10(limits=c(0.0001,1), breaks=c(0.0001,0.001,0.01,0.1,1), labels=c("0.0001", "0.001", "0.01", "0.1",1))+scale_fill_viridis(discrete=T, name="Quarantine effectiveness",labels=c('30%', '60%','90%'))+xlab('Action type')+scale_x_discrete(labels=c('Biological control', 'None', "Quarantine in", "Quarantine Out"))+ylab('Relative propagule pressure')+theme_classic()
plot


