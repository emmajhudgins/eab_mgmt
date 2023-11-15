rm(list=ls()) 
presences_time<-list()
missing<-to_rm<-list()
library(pdist)
#Read in Data
data<-read.csv('../data/countydatanorm_march.csv', stringsAsFactors = FALSE) # spatial data
data2<-read.csv('../data/datanorm.csv', stringsAsFactors = FALSE) #species-level data (including life history traits)
prez2<-prez<-matrix(0,3372, 64)
prez2[,1]<-c(which(data[,paste(data2$COLUMN_NAM[1])]>0), rep(0,(3372-length(which(data[,paste(data2$COLUMN_NAM[1])]>0))))) #determine where species are present
host.density2<-read.csv('../data/hostden_clean_gdk.csv') # tree host density by pest species
new_eab<-read.csv('../data/eab_each_2020.csv') # more recent detections
L<-rep(0,64) # size of each pest's host range
ic=F
V_i<-read.csv('../data/streettrees_grid.csv')[,20]
forecast=F
prez<-read.csv('../data/prez_clean_gdk.csv')
L<-length(unique(c(which(V_i>0), prez[,1], prez2[,1])))-1
prez[,1]<-c(unique(c(which(V_i!=0), prez[which(prez[,1]!=0),1], prez2[which(prez2[,1]!=0),1])), rep(0, 3372-L))
i<-1
#sources<-as.list(read.csv('Psources_notypos.csv')[,1])
currpopden<-as.matrix(read.csv("../data/currpopden_5.csv", stringsAsFactors = FALSE))
currpopden2<-as.matrix(read.csv("../data/future_scaled_pop2.csv"))
twenty35<-rowMeans(cbind(currpopden[,47], currpopden2[,4]))
twenty45<-rowMeans(cbind(currpopden2[,4], currpopden2[,5]))
twenty55<-rowMeans(cbind(currpopden2[,5], currpopden2[,6]))
currpopden<-cbind(currpopden, twenty35, currpopden2[,4], twenty45, currpopden2[,5], twenty55, currpopden2[,6])
#host.density2<-read.csv("hostden_clean_gdk.csv", stringsAsFactors = FALSE)

Tr1<-function(x)
{
  sqrt((data$X_coord-data$X_coord[x])^2+(data$Y_coord-data$Y_coord[x])^2)
}
dists<-sapply(1:3372, Tr1)
T1<-exp(-dists/50000)
rm(dists)
YEARS<-data2$YEAR
v_scale<-scale(V_i, center=T)
LLfit=function(par)
{
  pars<-rep(0,32)
  pars[1]<-par[1]
  pars[23:24]<-par[2:3]
  pars[c(21,22,4,18,20,8)]<-c(0.000538410692229749, 0.299034706404549, -0.525670755351726, 15.6132848183217,-0.163552592351765, 0.323831382884772) #best fit for GDKu on 5-year timescale
  par<-pars
  par[21]<-abs(par[21])
  par[22]<-abs(par[22])+1
  spp=i
  
  
  current_temp<-current_temp2<-0
  hum=0
  
  #Pest Parameters
  Pfull<<-matrix(0, 3372,64)
  Pfull_time<<-Pfull
  constpD=rep(0,64)
  constpD=matrix(rep(constpD),3372,64, byrow=TRUE)
  constpD2<-matrix(rep(par[9]*data[,18]+par[10]*data[,16]+par[16]*data[,19]+par[17]*data[,20]+par[18]*data[,21]+par[23]*v_scale[,1], 64),3372,64, byrow=F)+par[19]*host.density2
  constpD<-as.numeric(constpD)+constpD2
  constpD3<-matrix(rep(par[4]*data[,19]+par[12]*data[,20]+par[3]*data[,21]+par[13]*data[,18]+par[15]*data[,16]+par[24]*v_scale[,1], 64),3372,64, byrow=F)+par[5]*host.density2
  
  #Pest Parameters

  #sources<-as.list(read.csv('Psources_notypos.csv')[,1])
  #load("Psources_closest_march.Rdata")
  
  
  Psource=3347
  YEAR<-YEARS[spp]
  Discovery<-2009-YEAR
  Pfull<<-matrix(0, 3372, 64)
  T2<-T1[prez[1:L[spp],spp],prez[1:L[spp],spp]]
  vecP<-rep(0,L[spp])
  for (rrr in 1:length(Psource))
  {vecP[which(prez[,spp]==Psource[rrr])]=1}
  r0<-par[22]
  vecptime<<-matrix(0,L,(floor(YEAR/5)+9))
  for (time in 1:(floor(YEAR/5)+9))
  {
    vecP[which(prez[,spp]==Psource)]=1
    Pnext<-rep(0,L[spp]) # vecP for next timestep
    qq<-0 # dispersal kernel
    column<-(((Discovery+5*(time-1))-1790)/5)+1 # which year of human population density to consider
    qq<-matrix(rep(constpD[prez[which(vecP>=par[21]),spp],spp]+par[8]*currpopden[prez[which(vecP>=par[21]),spp],column], L[spp]), nrow=length(which(vecP>=par[21])), ncol=L[spp])
    zzz<-matrix(rep(constpD3[prez[1:L[spp],spp],spp]+par[20]*currpopden[prez[1:L[spp],spp],column], L[spp]), nrow=L[spp], ncol=L[spp], byrow=TRUE) # add in current human populations
    qq<-(2*par[1]*exp((zzz[which(vecP>=par[21]),]+qq)))/(1+exp(zzz[which(vecP>=par[21]),]+qq)) 
    qq<-T2[which(vecP>=par[21]),]^qq # add in distance
    #scale dispersal kernel
    if (length(which(vecP>=par[21]))>1){qq<-qq/rowSums(qq)}
    if (length(which(vecP>=par[21]))==1){qq<-qq/sum(qq)}
    qq[which(qq<0.001)]=0
    Pnext=(vecP[which(vecP>=par[21])])%*%(qq) # dispersal into and out of all sites
    
    Pnext[which(prez[,spp]==Psource)]=1
    Pfull_time[,time]<<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) #threshold at 'discoverable' value
    if (time>1)
    {
      dddd<-which(!(Pfull_time[1:length(which(Pfull_time[,time-1]!=0)),time-1]%in%Pfull_time[1:length(which(Pfull_time[,time]!=0)),time]))
      ffff<-which(prez[1:length(which(prez[,spp]!=0)), spp]%in%Pfull_time[dddd,time-1])
      Pnext[ffff]<-par[21]
    }
    if (forecast==T){
    if (time<floor(YEAR/5))
    {
      cccc<-which(!(prez[1:length(which(prez[,spp]!=0)),spp]%in%new_eab[1:length(which(new_eab[,4]!=0)),4]))
      Pnext[cccc]<-0  
    }
    if (time>=floor(YEAR/5) &time<=(floor(YEAR/5)+3))
    {
    dddd<-which(prez[1:length(which(prez[,spp]!=0)),spp]%in%new_eab[1:length(which(new_eab[,((time-2)*5)-1]!=0)),((time-2)*5)-1])
    cccc<-which(!(prez[1:length(which(prez[,spp]!=0)),spp]%in%new_eab[1:length(which(new_eab[,((time-2)*5)-1]!=0)),((time-2)*5)-1]))
    eeee<-which(Pnext[dddd]<par[21])
    Pnext[dddd[eeee]]<-par[21]
    Pnext[cccc]<-0    }
    }
    current_temp<-current_temp2
    Pnext[which(prez[,spp]==Psource)]=1
    Pnext[which(Pnext>=par[21])]=Pnext[which(Pnext>=par[21])]*r0
    Pnext[which(Pnext>=1)]<-1
    vecP=Pnext
    vecP[which(prez[,spp]==Psource)]=1
    Pfull_time[,time]<<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) #threshold at 'discoverable' value
    if (time==floor(YEAR/5))
   {
      Pfull[,spp]<<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) #threshold at 'discoverable' value  
    }
    vecptime[,time]<<-vecP
    
  }
  MET=function(spp) # MET calculation
  {
    dii<-2*sum(dist(cbind(data$X_coord[Pfull[1:length(which(Pfull[,spp]!=0)),spp]], data$Y_coord[Pfull[1:length(which(Pfull[,spp]!=0)),spp]]), upper=F))
    dij<-sum(pdist(cbind(data$X_coord[Pfull[1:length(which(Pfull[,spp]!=0)),spp]], data$Y_coord[Pfull[1:length(which(Pfull[,spp]!=0)),spp]]), cbind(data$X_coord[prez2[1:length(which(prez2[,spp]!=0)),spp]], data$Y_coord[prez2[1:length(which(prez2[,spp]!=0)),spp]]))@dist)
    djj<-2*sum(dist(cbind(data$X_coord[prez2[1:length(which(prez2[,spp]!=0)),spp]], data$Y_coord[prez2[1:length(which(prez2[,spp]!=0)),spp]])))
    return(((1/(length(which(Pfull[,spp]!=0))*length(which(prez2[,spp]!=0))))*dij)-((1/(2*length(which(Pfull[,spp]!=0))^2))*dii)-((1/(2*length(which(prez2[,spp]!=0))^2))*djj))
  }
  return (sum(tapply(spp, spp, MET)))
}

pars<-read.table("../data/par.1.adj_inter.1")[,1]
model<-optim(par=c(3.41, -0.421, -0.205), fn=LLfit, control=list(trace=100, maxit=1000))
#MET for full model=7.81
#MET for pars 1,2=23.34
#MET for pars 1,3=25.08
#MET for par 1=59.14
yy<-model$value # continue refitting until output is the same as input to avoid local minima
xx<-1000000
while (yy!=xx)
{
  yy<-model$value
  model<-optim(par=model$par, fn=LLfit, control=list(trace=100, maxit=1000))
  xx<-model$value
}
pars<-model$par
forecast=T
LLfit(c(pars))
presences_time[[i]]<-Pfull_time
saveRDS(model, '../output/full_model.RDS')
saveRDS(presences_time, file="../output/presences_time_eab.rds")


### plot model fit
forecast=F
pars<-readRDS('../output/full_model.rds')$par
LLfit(c(pars))
(length(which(Pfull[Pfull[,1]!=0,1]%in%prez2[,1]==T))+length((which(c(1:3060)%in%prez2[,1]==F)%in%Pfull[Pfull[,1]!=0,1])==F))/3060
#97.1% accuracy

# # Create two panels side by side
par(mar=c(2,0,0,4))
par(oma=c(0,0,0,0))
library(sp)
library(maptools)
library(maps)
m<-SpatialPointsDataFrame(coords=cbind(data$X_coord, data$Y_coord), data=data)

USA_merged<-map('usa', fill=TRUE, plot=F)
USAsstates<-map('state', fill=TRUE, plot=F, res=0)

IDs <- sapply(strsplit(USA_merged$names, ":"), function(x) x[1])
sp_map_usa <- map2SpatialPolygons(USA_merged, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))
transform_usa<-spTransform(sp_map_usa, CRS("+proj=eqdc +lat_0=39 +lon_0=-96 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))

IDsS <- sapply(strsplit(USAsstates$names, ":"), function(x) x[1])
sp_map_states <- map2SpatialPolygons(USAsstates, IDs=IDsS, proj4string=CRS("+proj=longlat +datum=WGS84"))
transform_states<-spTransform(sp_map_states, CRS("+proj=eqdc +lat_0=39 +lon_0=-96 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))

quar_bound<-readShapeLines('../plots/quarantine_boundary.shp')
proj4string(quar_bound)<-"+proj=longlat +datum=WGS84"
quar_bound<-spTransform(quar_bound,"+proj=eqdc +lat_0=39 +lon_0=-96 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")



pdf(paste("../plots/predvsobs_eab.pdf", sep="_"))
layout(matrix(c(1,2), ncol=2))
par(mai=c(0.25,0,0.25,0))
bins<-seq(-4,0,length.out=50)
bins<-10^bins
ash<-read.csv('streettrees_grid.csv')
colfunc <- colorRampPalette(c(alpha("grey",0), "darkred"))
for (time in c(6))
{
  col<-colfunc(50)[findInterval(vecptime[,time], bins)+1]
  plot(transform_usa, lwd=0.5, main="Predicted 2005 Density", col="grey")
  points(cbind(data$X_coord[prez[,1][which(vecptime[,time]>0)]], data$Y_coord[prez[,1][which(vecptime[,time]>0)]]), pch=15, cex=0.5, col=col[which(vecptime[,time]>0)])
  plot(transform_usa, lwd=0.5, fill=FALSE, add=T)
  par(xpd=T)
}
for (time in c(6))
{
  plot(transform_usa, lwd=0.5, main="Observed 2005 Detections",col="grey")
   points(cbind(data$X_coord[prez2[,1]], data$Y_coord[prez2[,1]]), pch=15, cex=0.5, col=colfunc(50)[25])
   plot(transform_usa, fill=FALSE, add=T)
   par(xpd=T)
}
dev.off()

