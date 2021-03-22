rm(list=ls()) 
n_spp=66

library(pdist)
correction=FALSE
setwd("~/Desktop/OneDrive - McGill University/Grad/scripts/")
#Read in Data
data<-read.csv('data_minus5_july.csv', stringsAsFactors = FALSE)
data2<-read.csv('datanorm.csv', stringsAsFactors = FALSE)
gen<-matrix(0,3372, 78)
for (sppp in 1:75)
{gen[,sppp]<-c(which(data[,paste(data2$COLUMN_NAM[sppp])]>0), rep(0,(3372-length(which(data[,paste(data2$COLUMN_NAM[sppp])]>0)))))}
rr<-which(data2$YEAR>=10)
data2<-data2[rr,]
gen<-gen[,rr]
FIA<-read.csv('FIAcodes_notypos.csv', stringsAsFactors = FALSE)
FIA<-FIA[rr,]
FIA2<-read.csv('FIA_march.csv', stringsAsFactors = FALSE)
fia<-list()
FIA$FIA<-as.character(FIA$FIA)
fia<-strsplit(FIA$FIA, split=", ")
currpopden<-as.matrix(read.csv("currpopden_5.csv", stringsAsFactors = FALSE))
currpopden<-as.matrix(read.csv("currpopden_5.csv", stringsAsFactors = FALSE))
currpopden2<-as.matrix(read.csv("future_scaled_pop2.csv"))
twenty05<-rowMeans(cbind(currpopden[,47], currpopden2[,1]))
twenty15<-rowMeans(cbind(currpopden2[,1], currpopden2[,2]))
twenty25<-rowMeans(cbind(currpopden2[,2], currpopden2[,3]))
currpopden<-cbind(currpopden, twenty05, currpopden2[,1], twenty15, currpopden2[,2], twenty25, currpopden2[,3])
sources<-as.list(read.csv('Psources_notypos.csv')[,1])
body<-read.csv('body_pred.csv', stringsAsFactors = F)
body$size[is.na(body$Source)==FALSE]<-scale(body$size[is.na(body$Source)==FALSE], center=TRUE)
mfive <- function(x){ 
 5*floor(x/5) 
} 
L<-rep(0,n_spp)
prez<-matrix(0,3372,n_spp)
to_manage<-readRDS('to_manage.RDS')
  for (spp in 1:n_spp) 
  { 
    cols<-paste("FIA", fia[[spp]], sep="_")
    pres<-rep(0, 3372)
    for (q in 1:length(cols))
    {
      if (cols[q] %in% colnames(FIA2))
      {
        nn<-which(FIA2[,cols[q]]!=0)
        pres[nn]=pres[nn]+1
      }
    }
    prez[,spp]<-c(which(pres!=0),rep(0, length(which(pres==0))))
    L[spp]<-length(which(pres!=0))
  }
  good<-which(L!=0)
  prez<-prez[,good]
  L<-L[good]
  data2<-data2[good,]
  gen<-gen[,good]
  host.density2<-read.csv("hostvol_notypos.csv", stringsAsFactors = FALSE)
  hostvol<-read.csv('hostvol_sdp.csv') #need to redo for 2 spp
  hostvol<-hostvol[,rr]
  hostvol<-hostvol[,good]
  publand<-read.csv('publand_sdp.csv')
  harvestland<-read.csv('harvest_sdp.csv')
  host_harv<-matrix(0,3372,75)
  for (spp in good) 
  { 
    cols<-paste("FIA", fia[[spp]], sep="_")
    pres<-rep(0, 3372)
    for (q in 1:length(cols))
    {
      if (cols[q] %in% colnames(FIA2))
      {nn<-which(FIA2[,cols[q]]!=0)
      host_harv[,spp]<-host_harv[,spp]+FIA2[,cols[q]]
      pres[nn]=pres[nn]+1}
    }
  }
  host_harv<-host_harv[,rr]
  host_harv<-host_harv[,good]
  fia<-fia[good]
  n_spp=length(data2[,1])
  prez2<-matrix(0,3372,64)
  L2<-rep(0,64)
  L3<-rep(0,64)
  for (sppp in 1:64)
  {prez2[,sppp]<-c(intersect(prez[which(prez[,sppp]!=0),sppp], gen[which(gen[,sppp]!=0),sppp]), rep(0,(3372-length(intersect(prez[which(prez[,sppp]!=0),sppp], gen[which(gen[,sppp]!=0),sppp])))))}
  # for (sppp in 1:64)
  # {L2[sppp]<-length(which(prez2[,sppp]!=0))}
  # for (sppp in 1:64)
  # {L3[sppp]<-length(which(gen[,sppp]!=0))}
  # Tr1<-function(x)
  # {
  #   sqrt((data$X_coord-data$X_coord[x])^2+(data$Y_coord-data$Y_coord[x])^2)
  # }
  # dists<-sapply(1:3372, Tr1)
  # T1<-exp(-dists/50000)
  # rm(dists)
  # YEARS<-data2$YEAR
  # pp_last<-matrix(0,3372,64) 
  # risk_pp<-matrix(0,3372,3372)
  # spp=1  
  #     par<-rep(0,23)
  #     startpt<-read.csv('startpt.csv')[,1]
  #     if (spp %in% startpt==F)
  #     {
  #       par[1]<-read.table(paste("./new_adj_inter/par",spp,"adj_inter",spp, sep="."))[,1]
  #     }
  #     if (spp %in% startpt==T)
  #     {
  #       par[1]<-read.csv(paste("./adj_inter_start/par_GDKic",spp,"csv", sep="."))[,2]
  #     }
  #     par[c(1,21,22,4,18,20,8)]<-as.numeric(c(par[1],c(0.000538410692229749, 0.299034706404549, -0.525670755351726, 15.6132848183217,-0.163552592351765, 0.323831382884772)))
  # 
  #     par[22]<-abs(par[22])+1
  # 
  #     YEAR<-YEARS[spp]
  #    #Pest Parameters
  #   total_time<-YEAR/5 +1
  #   Pfull<-matrix(0, 3372, total_time+5)
  #   Pfull_good<-matrix(0, 3372, total_time+5)
  #   Pfull_time<-Pfull
  #   vecP_time=d2prime=d3prime=d4prime=dprime=matrix(0,L[spp], total_time+5)
  # 
  #       constpD=rep(0,64)
  #       constpD=matrix(rep(constpD),3372,64, byrow=TRUE)
  #       constpD2<-matrix(rep(par[9]*data[,18]+par[10]*data[,16]+par[16]*data[,19]+par[17]*data[,20]+par[18]*data[,21]),3372,64)+par[19]*host.density2
  #       constpD<-as.numeric(constpD)+constpD2
  #       constpD3<-matrix(rep(par[4]*data[,19]+par[12]*data[,20]+par[3]*data[,21]+par[13]*data[,18]+par[15]*data[,16]),3372,64)+par[5]*host.density2
  #         
  #         #Pest Parameters
  #         Psource=sources[[spp]]
  #         Discovery<-2009-YEAR
  #         #Pfull<-matrix(0, 3372, total_time)
  #         rem<-mfive(Discovery)
  #         
  #         T2<-T1[prez[1:L[spp],spp],prez[1:L[spp],spp]]
  #         vecP<-rep(0,L[spp])
  #         for (rrr in 1:length(Psource))
  #         {vecP[which(prez[,spp]==Psource[rrr])]=1}
  # 
  #         r0<-par[22]
  # 
  #         for (time in 1:(total_time+5))
  #         {
  #           
  #           vecP[which(prez[,spp]==Psource)]=1
  #           dprime[,time]<-vecP
  #           d2prime[,time]<-vecP
  #           Pnext<-rep(0,L[spp])
  #           qq<-0
  #           column<-(((Discovery+5*(time-1))-1790)/5)+1
  # 
  #           qq<-matrix(rep(constpD[prez[which(vecP>=par[21]),spp],spp]+par[8]*currpopden[prez[which(vecP>=par[21]),spp],column], L[spp]), nrow=length(which(vecP>=par[21])), ncol=L[spp])
  #           zzz<-matrix(rep(constpD3[prez[1:L[spp],spp],spp]+par[20]*currpopden[prez[1:L[spp],spp],column], L[spp]), nrow=L[spp], ncol=L[spp], byrow=TRUE)
  #           
  #           qq<-(2*par[1]*exp((zzz[which(vecP>=par[21]),]+qq)))/(1+exp(zzz[which(vecP>=par[21]),]+qq))
  #           qq<-T2[which(vecP>=par[21]),]^qq
  #           if (length(which(vecP>=par[21]))>1){qq<-qq/rowSums(qq)}
  #           if (length(which(vecP>=par[21]))==1){qq<-qq/sum(qq)}
  #           qq[which(qq<0.001)]=0
  #           qq2<-matrix(0,L[spp], L[spp])
  #           qq2[which(vecP>=par[21]),]<-qq
  #           write.csv(qq2, file=paste("transmat_", spp,time, ".csv", sep=""), row.names=F) 
  #           Pnext=(vecP[which(vecP>=par[21])])%*%(qq)
  #           d3prime[,time]<-Pnext
  #           Pnext[which(prez[,spp]==Psource)]=1
  #           Pfull[,time]<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) 
  #           Pfull_time[,time]<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) #threshold at 'discoverable' value
  #           reset_yrs<-YEAR/5
  #             
  #           if (time>1)
  #           {
  # 
  #             dddd<-which(!(Pfull_time[1:length(which(Pfull_time[,time-1]!=0)),time-1]%in%Pfull_time[1:length(which(Pfull_time[,time]!=0)),time]))
  #             ffff<-which(prez[1:length(which(prez[,spp]!=0)), spp]%in%Pfull_time[dddd,time-1])
  #       
  #             Pnext[ffff]<-par[21]
  #             Pfull[,time]<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) 
  #             Pfull_time[,time]<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) #threshold at 'discoverable' value
  #           }
  #           if (time==floor(YEAR/5))
  #             {
  #               dddd<-which(prez[1:length(which(prez[,spp]!=0)),spp]%in%prez2[1:length(which(prez2[,spp]!=0)),spp])
  #               cccc<-which(!(prez[1:length(which(prez[,spp]!=0)),spp]%in%prez2[1:length(which(prez2[,spp]!=0)),spp]))
  #               eeee<-which(Pnext[dddd]<par[21])
  #               Pnext[dddd[eeee]]<-par[21] #testing out setting to true presences
  #               Pnext[cccc]<-0 #remove this when not doing sdp
  #               Pfull_time[,time]<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) 
  #           }
  #           
  #           Pnext[which(prez[,spp]==Psource)]=1
  #           Pnext[which(Pnext>=par[21])]=Pnext[which(Pnext>=par[21])]*r0
  #           d4prime[,time]<-Pnext
  #           Pnext[which(Pnext>=1)]<-1
  # 
  #           vecP=Pnext
  #           vecP[which(prez[,spp]==Psource)]=1
  #           vecP_time[,time]<-vecP
  #       }
  # 
  #   
# c_4=c_5=c_6=c_7=matrix(0,L[spp], ncol(vecP_time))
# for (i in 1:ncol(vecP_time))
# {
#   c_4[which(d2prime[,i]>0.000538410692229749),i]=1
#   d2prime[which(d2prime[,i]<=0.000538410692229749),i]<-0
#   c_5[which(d3prime[,i]>=0.000538410692229749),i]=1
#   c_6[which(d4prime[,i]<1),i]=1
#   c_7[which(d4prime[,i]>=1),i]=1
# }
# write.csv(dprime[,4:9], file=paste("dprime", ".csv", sep=""), row.names=F)
# write.csv(d2prime[,4:9], file=paste("d2prime", ".csv", sep=""), row.names=F)
# write.csv(d3prime[,4:9], file=paste("d3prime", ".csv", sep=""), row.names=F)
# write.csv(d4prime[,4:9], file=paste("d4prime", ".csv", sep=""), row.names=F)
# write.csv(vecP_time[,3:9], file=paste("vecptime",".csv", sep=""), row.names=F)
# write.csv(c_4[,4:9], file=paste("c_4",".csv", sep=""), row.names=F)
# write.csv(c_5[,4:9], file=paste("c_5",".csv", sep=""), row.names=F)
# write.csv(c_6[,4:9], file=paste("c_6",".csv", sep=""), row.names=F)
# write.csv(c_7[,4:9], file=paste("c_7",".csv", sep=""), row.names=F)
# qq<-as.matrix(read.csv('transmat_18.csv'))
# Pnext2=(d2prime[which(d2prime[,8]>=par[21]),8])%*%(qq[which(d2prime[,8]>=par[21]),])
# Pnext2[2]
# colSums(results[1:1799,]) # do nothing
# colSums(results[1800:(2*1799),]) # quar in
# colSums(results[(2*1799+1):(3*1799),]) # quar out
# colSums(results[(3*1799+1):(4*1799),]) #biological control
# colSums(results[(4*1799+1):(5*1799),]) # eradication
# figure out objective for rules of thumb
# maybe make some slides?
# ash_streettrees<-rowSums(readRDS('bestguess_genus.RDS')[,45,])
# write.csv(ash_streettrees, file="streettreevol_EAB.csv", row.names=F)

# # Create two panels side by side
par(mar=c(2,0,0,4))
par(oma=c(0,0,0,0))
library(sp)
library(maptools)
library(maps)
m<-SpatialPointsDataFrame(coords=cbind(data$X_coord, data$Y_coord), data=data)

USA_merged<-map('usa', fill=TRUE, plot=F)
USAsstates<-map('state', fill=F, plot=T)

IDs <- sapply(strsplit(USA_merged$names, ":"), function(x) x[1])
sp_map_usa <- map2SpatialPolygons(USA_merged, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))
transform_usa<-spTransform(sp_map_usa, CRS("+proj=eqdc +lat_0=39 +lon_0=-96 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))

IDsS <- sapply(strsplit(USAsstates$names, ":"), function(x) x[1])
sp_map_states <- map2SpatialPolygons(USAsstates, IDs=IDsS, proj4string=CRS("+proj=longlat +datum=WGS84"))
transform_states<-spTransform(sp_map_states, CRS("+proj=eqdc +lat_0=39 +lon_0=-96 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))


time<-seq(1:5)
quar<-c(0.3,0.6,0.9)
bio<-c(0.1,0.3,0.5)
# for (i in quar)
# {
#   for (j in quar){
#     for (k in bio){
#         qin=i
#         qout=j
#         qbio=k
qin=qout=0.9
qbio=0.5
results<-read.csv(paste0('./RoT/management_test_1_0_',qin,'_',qout,'_',qbio,'_.csv'))
pdf(paste("eab_mgmt_rot_site", qin, qout, qbio,".pdf", sep="_"))
layout(matrix(c(1,2,3,4,5,6), ncol=2, byrow=TRUE), heights=c(1,1,1))
par(mai=c(0.25,0,0.25,0))
years<-seq(2025,2045, by=5)
for (time in 1:5)
{
plot(transform_usa, lwd=0.5, main=years[time], col="grey")
points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5)
points(cbind(data$X_coord[prez[which(results[((1799+1):(2*1799)), time]==1),1]], data$Y_coord[prez[which(results[((1799+1):(2*1799)), time]==1),1]]), pch=15, cex=0.5, col="yellow")
points(cbind(data$X_coord[prez[which(results[((2*1799+1):(3*1799)), time]==1),1]], data$Y_coord[prez[which(results[((2*1799+1):(3*1799)), time]==1),1]]), pch=15, cex=0.5, col="orange")
points(cbind(data$X_coord[prez[which(results[((3*1799+1):(4*1799)), time]==1),1]], data$Y_coord[prez[which(results[((3*1799+1):(4*1799)), time]==1),1]]), pch=15, cex=0.5, col="red")
points(cbind(data$X_coord[prez[which(results[((4*1799+1):(5*1799)), time]==1),1]], data$Y_coord[prez[which(results[((4*1799+1):(5*1799)), time]==1),1]]), pch=15, cex=0.5, col="darkred")
}
plot.new()
par(xpd=FALSE)
legend('center', c("Quarantine In", "Quarantine Out", "Biocontrol","Eradicate"), col=c('yellow', 'orange', 'red', 'darkred'), pch=15, cex=2)
dev.off()

#results<-read.csv('pestden3_nomgmt_new.csv')

results<-read.csv(paste0('pestden3_',qin,'_',qout,'_',qbio,'_new.csv'))
#results[,4]<-results[,4]-results[,3]
#results[,3]<-results[,3]-results[,2]
#results[,2]<-results[,2]-results[,1]
#results[,2][which(results[,1]>=0.000538)]<-0
#results[,3][which(results[,2]>=0.000538)]<-0
#results[,4][which(results[,3]>=0.000538)]<-0
#results[,4][which(results[,4]<0.000538)]<-0
#results[,3][which(results[,3]<0.000538)]<-0
#results[,2][which(results[,2]<0.000538)]<-0

time<-seq(1:6)
years<-seq(2025,2050, by=5)
pdf(paste("eab_dens", qin, qout, qbio,".pdf", sep="_"))
layout(matrix(c(1,2,3,4,5,6,7,7), ncol=2, byrow=TRUE), heights=c(0.3,0.3,0.3,0.1))
par(mai=c(0.25,0,0.25,0))
years<-seq(2025,2050, by=5)
 bins<-seq(5,0,length.out=50)
 bins<-10^-bins
 library(viridis)
library(ggplot2)
for (time in 1:6)
{
col<-viridis(51)[findInterval(results[,time], bins)+1]
plot(transform_usa, lwd=0.5, main=years[time], col="grey")
points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5)
points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5, col=col)
#plot(transform_states, lwd=0.5, main=years[time], fill=FALSE, add=T,border='white')
par(xpd=T)
}
 par(xpd=FALSE)
 par(mai=c(0.5,0.25,0,0.25))
 image(1:51,1,matrix(1:51), col=viridis(51), axes=FALSE, ann=F)
 axis(1,labels=c("0",expression(10^{-5}), expression(10^{-3}),expression(10^{-1}),expression(10^{0})),at=c(seq(0,51,length.out=5)), cex.axis=1)
 mtext(side=1, "EAB propagule pressure", line=2)
dev.off()

pdf(paste("eab_exp", qin, qout, qbio,".pdf", sep="_"))
layout(matrix(c(1,2,3,4,5,6,7,7), ncol=2, byrow=TRUE), heights=c(0.3,0.3,0.3,0.1))
par(mai=c(0.25,0,0.25,0))
years<-seq(2025,2050, by=5)
bins<-seq(-5,5,length.out=50)
bins<-10^bins
ash<-read.csv('~/Desktop/OneDrive - McGill University/Grad/scripts/streettrees_grid.csv')
for (time in 1:6)
{
  col<-viridis(50)[findInterval(results[,time]*ash[prez[1:1799,1],20], bins)+1]
  plot(transform_usa, lwd=0.5, main=years[time], col="grey")
  points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5)
  points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5, col=col)
  par(xpd=T)
}
par(xpd=FALSE)
par(mai=c(0.5,0.25,0,0.25))
image(1:50,1,matrix(1:50), col=viridis(50), axes=FALSE, ann=F)
axis(1,labels=c("0",expression(10^{2}), expression(10^{4}),expression(10^{6})),at=c(seq(0,50,length.out=4)), cex.axis=1)
mtext(side=1, "Exposed Street Ash", line=2)
dev.off()
  #  }}}
#pdf("totalash.pdf")
layout(matrix(c(1,2), ncol=1, byrow=TRUE), heights=c(0.75,0.25))
bins<-seq(0,6,length.out=50)
bins<-10^bins
ash<-read.csv('~/Desktop/OneDrive - McGill University/Grad/scripts/streettrees_grid.csv')
col<-viridis(50)[findInterval(ash[prez[1:1799,1],20], bins)+1]
par(mai=c(0,0,0,0))
par(mar=c(0,0,0,0))
par(oma=c(0,0,0,0))
  plot(transform_usa, lwd=0.5, main=NULL, col="grey")
  points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5)
  points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5, col=col)
  par(mai=c(0.75,0.25,0,0.25))
    image(1:50,1,matrix(1:50), col=viridis(50), axes=FALSE, ann=F)
  axis(1,labels=c("0",expression(10^{2}), expression(10^{4}),expression(10^{6})),at=c(seq(0,50,length.out=4)), cex.axis=1)
  mtext(side=1, "Street Ash", line=2)
#dev.off()

  layout(matrix(c(1,2), ncol=1, byrow=TRUE), heights=c(0.75,0.25))
  bins<-seq(0,8,length.out=50)
  bins<-10^bins
  ash<-read.csv('~/Desktop/OneDrive - McGill University/Grad/scripts/hostvol_sdp2020.csv')
  col<-viridis(50)[findInterval(ash[prez[1:1799,1],1], bins)+1]
  par(mai=c(0,0,0,0))
  par(mar=c(0,0,0,0))
  par(oma=c(0,0,0,0))
  plot(transform_usa, lwd=0.5, main=NULL, col="grey")
  points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5)
  points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5, col=col)
  par(mai=c(0.75,0.25,0,0.25))
  image(1:50,1,matrix(1:50), col=viridis(50), axes=FALSE, ann=F)
  axis(1,labels=c("0",expression(10^{2}), expression(10^{4}),expression(10^{6}),expression(10^{8})),at=c(seq(0,50,length.out=5)), cex.axis=1)
  mtext(side=1, "Forest Ash Volume", line=2)
dat<-read.csv('~/Downloads/postdocdat.csv')
