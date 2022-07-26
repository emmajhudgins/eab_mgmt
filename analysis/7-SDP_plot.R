rm(list=ls()) 
n_spp=66
library(sp)
library(pdist)
library(gpplot2)
library(here)
setwd(paste0(here(),'../data/'))

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
for (i in 1:66)
{
  fia[[i]]<-gsub("^", "FIA_", fia[[i]])
  fia[[i]]<-gsub("_(\\d{3}$)", "_0\\1", fia[[i]])
  fia[[i]]<-gsub("_(\\d{2}$)", "_00\\1", fia[[i]])
}
currpopden<-as.matrix(read.csv("currpopden_5.csv", stringsAsFactors = FALSE))
currpopden2<-as.matrix(read.csv("future_scaled_pop2.csv"))
twenty05<-rowMeans(cbind(currpopden[,47], currpopden2[,1]))
twenty15<-rowMeans(cbind(currpopden2[,1], currpopden2[,2]))
twenty25<-rowMeans(cbind(currpopden2[,2], currpopden2[,3]))
currpopden<-cbind(currpopden, twenty05, currpopden2[,1], twenty15, currpopden2[,2], twenty25, currpopden2[,3])
sources<-as.list(read.csv('Psources_notypos.csv')[,1])
L<-rep(0,n_spp)
prez<-matrix(0,3372,n_spp)
for (spp in 1:n_spp) 
{ 
  cols<-paste(fia[[spp]])
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
fia<-fia[good]
n_spp=length(data2[,1])
prez2<-matrix(0,3372,64)
L2<-rep(0,64)
L3<-rep(0,64)
for (sppp in 1:64)
{prez2[,sppp]<-c(intersect(prez[which(prez[,sppp]!=0),sppp], gen[which(gen[,sppp]!=0),sppp]), rep(0,(3372-length(intersect(prez[which(prez[,sppp]!=0),sppp], gen[which(gen[,sppp]!=0),sppp])))))}

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

quar_bound<-readShapeLines('~/Desktop/OneDrive - McGill University/GitHub/eab_mgmt/plots/quarantine_boundary.shp')
proj4string(quar_bound)<-"+proj=longlat +datum=WGS84"
quar_bound<-spTransform(quar_bound,"+proj=eqdc +lat_0=39 +lon_0=-96 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")


time<-seq(1:5)
quar<-c(0.3,0.6,0.9)
bio<-c(0.1,0.3,0.5)

qin=qout=0.3
qbio=0.1
frac_site=0.2
frac_spread=0.8
results<-read.csv(paste0('~/Desktop/OneDrive - McGill University/GitHub/eab_mgmt/output/M_',qin,'_',qbio,'.csv'), header=F)

pdf(paste("~/Desktop/OneDrive - McGill University/GitHub/eab_mgmt/plots/eab_mgmt_site",qin, qbio,".pdf", sep="_"))
layout(matrix(c(1,2,3,4,5,6,7,7,7), ncol=3, byrow=T), heights=c(0.5,0.5,0.4))
par(mai=c(0.25,0,0.25,0))
years<-seq(2025,2045, by=5)
for (time in c(1,3,5))
{
  plot(transform_usa, lwd=0.5, main=years[time], col="white")
  points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5, col="grey")
  plot(transform_usa, lwd=0.5, main=years[time], col=alpha("grey",0), add=T)
  
  lines(quar_bound, col=alpha("deeppink4",1), lwd=2)
  points(cbind(data$X_coord[prez[which(results[((1799+1):(2*1799)), time]==1),1]], data$Y_coord[prez[which(results[((1799+1):(2*1799)), time]==1),1]]), pch=15, cex=0.5, col="yellow")
  points(cbind(data$X_coord[prez[which(results[((2*1799+1):(3*1799)), time]==1),1]], data$Y_coord[prez[which(results[((2*1799+1):(3*1799)), time]==1),1]]), pch=15, cex=0.5, col="orange")
  points(cbind(data$X_coord[prez[which(results[((3*1799+1):(4*1799)), time]==1),1]], data$Y_coord[prez[which(results[((3*1799+1):(4*1799)), time]==1),1]]), pch=15, cex=0.5, col="red")

  }
results<-read.csv(paste0('~/Desktop/OneDrive - McGill University/GitHub/eab_mgmt/output/management_test_0.2_0.8_',qin,'_',qbio,'_.csv'))
for (time in c(1,3,5))
{
  plot(transform_usa, lwd=0.5, main=years[time], col="white")
  points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5, col="grey")
  lines(quar_bound, col=alpha("deeppink4",1), lwd=2)
  plot(transform_usa, lwd=0.5, main=years[time], col=alpha("grey",0), add=T)
  points(cbind(data$X_coord[prez[which(results[((1799+1):(2*1799)), time]==1),1]], data$Y_coord[prez[which(results[((1799+1):(2*1799)), time]==1),1]]), pch=15, cex=0.5, col="yellow")
  points(cbind(data$X_coord[prez[which(results[((2*1799+1):(3*1799)), time]==1),1]], data$Y_coord[prez[which(results[((2*1799+1):(3*1799)), time]==1),1]]), pch=15, cex=0.5, col="orange")
  points(cbind(data$X_coord[prez[which(results[((3*1799+1):(4*1799)), time]==1),1]], data$Y_coord[prez[which(results[((3*1799+1):(4*1799)), time]==1),1]]), pch=15, cex=0.5, col="red")

  }
plot.new()
par(xpd=FALSE)
legend('center', c("Quarantine In", "Quarantine Out", "Biocontrol", "Previous Quarantine Boundary"), col=c('yellow', 'orange', 'red', "deeppink4"), pch=15, cex=2)
dev.off()

#results<-read.csv('pestden3_nomgmt_new.csv')

#results[,4]<-results[,4]-results[,3]
#results[,3]<-results[,3]-results[,2]
#results[,2]<-results[,2]-results[,1]
#results[,2][which(results[,1]>=0.000538)]<-0
#results[,3][which(results[,2]>=0.000538)]<-0
#results[,4][which(results[,3]>=0.000538)]<-0
#results[,4][which(results[,4]<0.000538)]<-0
#results[,3][which(results[,3]<0.000538)]<-0
#results[,2][which(results[,2]<0.000538)]<-0
qin=qout=0.3
qbio=0.1
frac_site=0.2
frac_spread=0.8
results<-read.csv(paste('~/Desktop/OneDrive - McGill University/GitHub/eab_mgmt/output/vecptime_0.9_0.5.csv'), header=F)/1000
time<-seq(1:6)
years<-seq(2020,2045, by=5)
pdf(paste("~/Desktop/OneDrive - McGill University/GitHub/eab_mgmt/plots/RoT_eab_dens_",qin, qout, qbio,".pdf", sep="_"))
years<-seq(2020,2045, by=5)
bins<-seq(5,0,length.out=50)
bins<-10^-bins
library(viridis)
library(ggplot2)
layout(matrix(c(1,2,3,4,5,6,7,7,7), ncol=3, byrow=T), heights=c(0.5,0.5,0.1))
par(mai=c(0,0,0.25,0))
for (time in c(2,4,6))
{
  col<-viridis(51)[findInterval(results[,time], bins)+1]
  plot(transform_usa, lwd=0.5, main=years[time], col="grey")
  points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5)
  points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5, col=col)
  #plot(transform_states, lwd=0.5, main=years[time], fill=FALSE, add=T,border='white')
  par(xpd=T)
}
results<-read.csv(paste0("~/Desktop/OneDrive - McGill University/GitHub/eab_mgmt/output/RoT_vecptime_", frac_site,"_",frac_spread, "_", qin,'_',qbio,'_.csv'))
for (time in c(2,4,6))
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

bins<-seq(-5,5,length.out=50)
bins<-10^bins


results<-read.csv(paste('~/Desktop/OneDrive - McGill University/GitHub/eab_mgmt/output/vecptime_0.3_0.1.csv'), header=F)/1000
pdf(paste("~/Desktop/OneDrive - McGill University/GitHub/eab_mgmt/plots/eab_exp",qin,qbio,".pdf", sep="_"))
layout(matrix(c(1,2,3,4,5,6,7,7,7), ncol=3, byrow=T), heights=c(0.5,0.5,0.1))
par(mai=c(0.25,0,0.25,0))
years<-seq(2020,2045, by=5)
bins<-seq(-5,5,length.out=50)
bins<-10^bins
ash<-read.csv('~/Desktop/OneDrive - McGill University/GitHub/eab_mgmt/data/streettrees_grid.csv')
 for (time in c(2,4,6))
  {
   col<-viridis(50)[findInterval(results[,time]*ash[prez[1:1799,1],20], bins)+1]
   plot(transform_usa, lwd=0.5, main=years[time], col="grey")
    points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5)
    points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5, col=col)
    #plot(transform_states, lwd=0.5, main=years[time], fill=FALSE, add=T,border='white')
    par(xpd=T)
  }
  results<-read.csv(paste0("~/Desktop/OneDrive - McGill University/GitHub/eab_mgmt/output/RoT_vecptime_", frac_site,"_",frac_spread, "_", qin,'_',qbio,'_.csv'))
  for (time in c(2,4,6))
  {
    col<-viridis(50)[findInterval(results[,time]*ash[prez[1:1799,1],20], bins)+1]
    plot(transform_usa, lwd=0.5, main=years[time], col="grey")
    points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5)
    points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5, col=col)
    #plot(transform_states, lwd=0.5, main=years[time], fill=FALSE, add=T,border='white')
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

qz<-c(0.3,0.6,0.9)
bios<-c(0.1,0.3,0.5)
budget_scen<-data.frame(site_bud=seq(0,1, length.out=11), spread_bud=seq(1,0,length.out=11))
obj<-data.frame(frac_site=0,frac_spread=0,q_in=0,q_out=0,qbio=0, obj=0)
names(obj)
V_i<-read.csv('../../GitHub/eab_mgmt/data/streettrees_grid.csv')[,20]

for (q_in in qz)
{
q_out=q_in
    for (qbio in bios)
    {
      for (scen in c(1,3,5,7,9,11))
      {
        frac_site=budget_scen$site_bud[scen]
        frac_spread=budget_scen$spread_bud[scen]
        d4prime<-read.csv(paste("~/Desktop/OneDrive - McGill University/GitHub/eab_mgmt/output/RoT_vecptime",frac_site,frac_spread,q_in,qbio, ".csv", sep="_"))[,1:7]
        obj<-rbind(obj, setNames(c(frac_site,frac_spread,q_in,q_out,qbio,sum(sweep(as.matrix(d4prime),MARGIN=1,as.vector(V_i[prez[,1]]+1),"*"))),names(obj)))
      }
    }
  }

obj<-obj[2:nrow(obj),]
library(viridis)
library(ggplot2)

layout(matrix(c(1,2), ncol=2, byrow=TRUE), heights=c(1,1), widths=c(0.7,0.3))
par(mar=c(4,4,2,0))
par(oma=c(0,0,0,0))
par(xpd = FALSE)

plot(y=obj$obj,x=(obj$frac_site), col=alpha(viridis(9)[as.factor(paste0(obj$q_in, obj$qbio))],0.75),xlab="Biological control proportion", ylab="Exposed ash street trees", pch=19, xlim=c(-0.05,1.05), ylim=c(100000,1300000), axes=F)
axis(1,  at=c(0,0.2,0.4,0.6,0.8,1.0))
axis(2, labels=c('0',"200000", "400000", '600000', '800000', '1000000', '1200000', '1400000'), at=c(0,200000,400000,600000,800000, 1000000,1200000, 1400000), outer=F, las=2, cex.axis=0.5)

library(readxl)
dat<-read.csv('~/Downloads/postdocdat.csv')
points(y=dat$Exposed.Ash.Street.Trees..thousands.*1000, x=dat$bio_prop/100, lty=2,col=viridis(9), pch=5)
par(mai=c(0.4,0.5,0.1,0.6))

image(1,1:9,t(matrix(1:9)), col=viridis(9), axes=FALSE, ann=F)
axis(2,labels=c("30%","30%","30%","60%","60%","60%", "90%","90%","90%"),at=c(1:9), cex.axis=0.5, padj=2)
axis(4,at=c(1:9), labels=rep("", 9),cex.axis=0.5, padj=2)
par(xpd = TRUE) #Draw outside plot area
corners<-par("usr")
text("50%             30%             10%",y=c(8), x=corners[2]+.2, cex=0.5, srt=270)
text("50%             30%             10%",y=c(5), x=corners[2]+.2, cex=0.5, srt=270)
text("50%             30%             10%",y=c(2), x=corners[2]+.2, cex=0.5, srt=270)

mtext(side=2, "Quarantine efficiency",line=1, cex=0.75)
corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
text(x = corners[2]+0.4, y = mean(corners[3:4]),"Biological Control efficiency",cex=0.7, srt = 270)




qz<-c(0.3,0.6,0.9)
bios<-c(0.1,0.3,0.5)
obj<-data.frame(q_in=0,qbio=0, obj=0)
names(obj)

for (q_in in qz)
{
    for (qbio in bios)
    {
    d4prime<-read.csv(paste("~/Desktop/OneDrive - McGill University/GitHub/eab_mgmt/output/vecptime_",q_in,"_",qbio, ".csv", sep=""), header=F)[,2:8]/1000
        obj<-rbind(obj, setNames(c(q_in,qbio,sum(sweep(as.matrix(d4prime),MARGIN=1,as.vector(V_i[prez[,1]]+1),"*"))),names(obj)))
      }
  }
layout(1,1)
obj<-obj[2:nrow(obj),]
library(viridis)
library(ggplot2)
plot(y=obj$obj,x=obj$q_in, col=alpha(viridis(27)[as.factor(paste0(obj$q_in,obj$qbio))],0.75), xlab="Quarantine efficiency", ylab="Exposed ash street trees", pch=19, ylim=c(100000,1200000))
par(xpd = FALSE)
par(mar=c(4,4,2,2))
plot(dat$Exposed.Ash.Street.Trees..thousands.*1000~dat$Biocontrol.efficiency,col=viridis(27)[c(1,1,1,13,13,13,25,25,25,25)], pch=19, xlab="Biocontrol Efficiency", ylab="Exposed ash street trees", ylim=c(100000,1500000), axes=F) 
axis(1,  at=c(0,0.1,0.3,0.5))
axis(2, labels=c('0',"200000", "400000", '600000', '800000', '1000000', '1200000', '1400000'), at=c(0,200000,400000,600000,800000, 1000000,1200000, 1400000), outer=F, las=2, cex.axis=0.5)

points(x=c(0.1,0.3,0.5), y=c(1239440,1050990,767630),col="red", pch=24) 
legend('topright', legend=c("30% Quarantine Efficiency","60% Quarantine Efficiency", "90% Quarantine Efficiency", "Biocontrol Only"), pch=c(19,19,19,24), col=c(viridis(27)[c(1,13,25)], "red"), cex=0.75)
# points(y=rep(978764.3,5),x=c(0.5403754, 0.3209943, 0.0000000, 0.3491804, 0.0000000 ), col=viridis(27)[1], pch=c("1","2","3","4","5"))

# abline(h=946951, lty=2,col=viridis(27)[4])
# abline(h=946951, lty=2,col=viridis(27)[7])
# abline(h=887046.1, lty=2,col=viridis(27)[16])
# points(y=rep(887046.1,5),x=c(0.3692045, 0.2444020, 0.2881563, 0.2918899, 0.0000000 ), col=viridis(27)[16], pch=c("1","2","3","4","5")) 
# # abline(h=946951, lty=2,col=viridis(27)[13])
# # abline(h=946951, lty=2,col=viridis(27)[16])
# abline(h=, lty=2,col=viridis(27)[27])
# points(y=rep(797804,5),x=c(5.403754e-01, 3.209943e-01, 2.425140e-07, 3.491806e-01, 1.734203e-07), col=viridis(27)[27], pch=c("1","2","3","4","5"))

# abline(h=946951, lty=2,col=viridis(27)[22])
# abline(h=946951, lty=2,col=viridis(27)[25])

#points(y=rep(946951,5),x=cost_each[,6], col=viridis(27)[1], pch=8)
#legend('bottomright', legend=c("Site first", "Spread First", "Optimality"), pch=c(19,23,8), col=viridis(1))
#legend('bottomright', legend=c("Rule of Thumb", "Optimality"), pch=c(19,NA), lty=c(NA,2), col=viridis(27)[13], cex=0.5)
mgmt_itme<-read.csv('M3_0.9_0.9_0.1.csv', header=F)
mgmt_itme<-round(mgmt_itme)
d_time<-read.csv('pestden3_0.6_0.6_0.3_new.csv')
# d_time<-read.csv('../vecptime_0.9_0.1_0.3_0.3_0.1_

qin=qout=0.3
qbio=0.1
mgmt<-read.csv(paste('~/Desktop/OneDrive - McGill University/GitHub/eab_mgmt/output/M_',qin,"_",qbio,'.csv', sep=""), header=F)
results<-(read.csv(paste('~/Desktop/OneDrive - McGill University/GitHub/eab_mgmt/output/vecptime_',qin,"_",qbio,'.csv', sep=""), header=F)/1000)[,1:6]
no_action<-(as.matrix(results)[,2])[which(mgmt[1:1799,2]==1)]
q_in<-(as.matrix(results[,2]))[which(mgmt[1800:(2*1799),2]==1)]
q_out<-(as.matrix(results[,2]))[which(mgmt[(2*1799+1):(3*1799),2]==1)]
bio<-(as.matrix(results))[which(mgmt[(3*1799+1):(4*1799),]==1)]
df <- data.frame(values = c(no_action,q_in,q_out,bio),
                 vars = rep(c("no_action","q_in", "q_out","bio"), times = c(length(no_action), length(q_in), length(q_out), length(bio))))
par(mar=c(4,4,2,2))
boxplot(values~vars, df)
library(vioplot)
vioplot(values~vars, df,horizontal=T)


qz<-c(0.3,0.6,0.9)
bios<-c(0.1,0.3,0.5)
obj<-data.frame(q_in=0,qbio=0, obj=0)
names(obj)

for (q_in in qz)
{
  for (qbio in bios)
  {
    d4prime<-read.csv(paste("~/Desktop/OneDrive - McGill University/GitHub/eab_mgmt/output/vecptime_",q_in,"_",qbio, "_fdfg.csv", sep=""), header=F)[,2:8]/1000
    obj<-rbind(obj, setNames(c(q_in,qbio,sum(sweep(as.matrix(d4prime),MARGIN=1,as.vector(V_i[prez[,1]]+1),"*"))),names(obj)))
  }
}

obj<-obj[2:nrow(obj),]


qz<-c(0.3,0.6,0.9)
bios<-c(0.1,0.3,0.5)
budget_scen<-data.frame(site_bud=seq(0,1, length.out=6), spread_bud=seq(1,0,length.out=6))
obj<-data.frame(frac_site=0,q_in=0,qbio=0, obj=0)
names(obj)
V_i<-read.csv('../../GitHub/eab_mgmt/data/streettrees_grid.csv')[,20]

for (q_in in qz)
{
    for (qbio in bios)
    {
      for (scen in 1:6)
      {
        frac_site=budget_scen$site_bud[scen]
        frac_spread=budget_scen$spread_bud[scen]
    
        d4prime<-read.csv(paste("~/Desktop/OneDrive - McGill University/GitHub/eab_mgmt/output/vecptime",frac_spread,q_in,qbio, "bud.csv", sep="_"), header=F)[,2:8]/1000
        obj<-rbind(obj, setNames(c(frac_site, q_in,qbio,sum(sweep(as.matrix(d4prime),MARGIN=1,as.vector(V_i[prez[,1]]+1),"*"))),names(obj)))
      }
    }
  }

obj<-obj[2:nrow(obj),]

layout(matrix(c(1,2), ncol=2, byrow=TRUE), heights=c(1,1), widths=c(0.7,0.3))
par(mar=c(4,4,2,0))
par(oma=c(0,0,0,0))
par(xpd = FALSE)


plot(y=obj$obj,x=(obj$frac_site), col=alpha(viridis(9)[as.factor(paste0(obj$q_in, obj$qbio))],0.75), xlab="Biological control proportion", ylab="Exposed ash street trees",  pch=19, xlim=c(-0.05,1.05), ylim=c(0,1400000), cex.main=0.75, axes=F)
axis(1,labels=c("0","20%", "40%","60%", '80%', "100%"), at=c(0,0.2,0.4,0.6,0.8,1.0), outer=F)
axis(2, labels=c('0',"200000", "400000", '600000', '800000', '1000000','1200000','1400000'), at=c(0,200000,400000,600000,800000, 1000000,1200000,1400000), outer=F, las=2,cex.axis=0.5)
library(readxl)
dat<-read.csv('~/Downloads/postdocdat.csv')
points(y=dat$Exposed.Ash.Street.Trees..thousands.*1000, x=dat$bio_prop/100, lty=2,col=viridis(9), pch=5)
par(mai=c(0.4,0.5,0.1,0.6))

image(1,1:9,t(matrix(1:9)), col=viridis(9), axes=FALSE, ann=F)
axis(2,labels=c("30%","30%","30%","60%","60%","60%", "90%","90%","90%"),at=c(1:9), cex.axis=0.5, padj=2)
axis(4,at=c(1:9), labels=rep("", 9),cex.axis=0.5, padj=2)
par(xpd = TRUE) #Draw outside plot area
corners<-par("usr")
text("50%             30%             10%",y=c(8), x=corners[2]+.2, cex=0.5, srt=270)
text("50%             30%             10%",y=c(5), x=corners[2]+.2, cex=0.5, srt=270)
text("50%             30%             10%",y=c(2), x=corners[2]+.2, cex=0.5, srt=270)

mtext(side=2, "Quarantine efficiency",line=1, cex=0.75)
corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
text(x = corners[2]+0.4, y = mean(corners[3:4]),"Biological Control efficiency",cex=0.7, srt = 270)



layout(matrix(c(1,2), ncol=2, byrow=TRUE), heights=c(1,1), widths=c(0.7,0.3))
par(mar=c(4,4,2,0))
par(oma=c(0,0,0,0))
par(xpd = FALSE)

plot(y=obj$obj,x=(obj$frac_site), col=alpha(viridis(9)[as.factor(paste0(obj$q_in, obj$qbio))],0.75),xlab="Biological control proportion", ylab="Exposed ash street trees", pch=19, xlim=c(-0.05,1.05), ylim=c(100000,1200000))

library(readxl)
dat<-read.csv('~/Downloads/postdocdat.csv')
points(y=dat$Exposed.Ash.Street.Trees..thousands.*1000, x=dat$bio_prop/100, lty=2,col=viridis(9), pch=5)
par(mai=c(0.4,0.5,0.1,0.6))
image(1,1:9,t(matrix(1:9)), col=viridis(9), axes=FALSE, ann=F)
axis(2,labels=c("30%","60%", "90%"),at=c(seq(1,9,length.out=3)), cex.axis=0.5, padj=2)
axis(4,at=c(1:9), labels=rep("", 9),cex.axis=0.5, padj=2)
par(xpd = TRUE) #Draw outside plot area
corners<-par("usr")
text("50%             30%             10%",y=c(8), x=corners[2]+.65, cex=0.5, srt=270)
text("50%             30%             10%",y=c(5), x=corners[2]+.65, cex=0.5, srt=270)
text("50%             30%             10%",y=c(2), x=corners[2]+.65, cex=0.5, srt=270)

mtext(side=2, "Quarantine efficiency",line=1, cex=0.5)
corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
text(x = corners[2]+1.25, y = mean(corners[3:4]),"Biological Control efficiency",cex=0.5, srt = 270)



plot(dat$Exposed.Ash.Street.Trees..thousands.*1000~dat$Biocontrol.efficiency,col=viridis(27)[c(1,1,1,13,13,13,25,25,25,25)], pch=19, xlab="Biocontrol Efficiency", ylab="Exposed Ash Trees", ylim=c(100000,1500000)) 
points(x=c(0.1,0.3,0.5), y=c(1233200,1035330, 729520),col="red", pch=23) 
legend('topright', legend=c("All Management", "Biocontrol Only"), pch=c(19,23), col=c(viridis(27)[c(13)], "red"), cex=0.75)


bins<-seq(-5,5,length.out=50)
bins<-10^bins


results<-read.csv(paste('~/Desktop/OneDrive - McGill University/GitHub/eab_mgmt/data/vecptime_0_1_0.3_0.1_.csv'))
pdf(paste("~/Desktop/OneDrive - McGill University/GitHub/eab_mgmt/plots/ash",qin,qbio,".pdf", sep="_"))
layout(matrix(c(1,2,3,3), ncol=2, byrow=T),heights=c(0.8,0.2))
par(mai=c(0.25,0,0.25,0))
bins<-seq(-5,5,length.out=50)
bins<-10^bins
ash<-read.csv('~/Desktop/OneDrive - McGill University/GitHub/eab_mgmt/data/streettrees_grid.csv')
for (time in c(6))
{
  col<-viridis(50)[findInterval(ash[prez[1:1799,1],20], bins)+1]
  plot(transform_usa, lwd=0.5, main="Total", col="grey")
  points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5)
  points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5, col=col)
  #plot(transform_states, lwd=0.5, main=years[time], fill=FALSE, add=T,border='white')
  par(xpd=T)
}
results<-read.csv("~/Desktop/OneDrive - McGill University/GitHub/eab_mgmt/data/vecP_noaction.csv")
for (time in c(6))
{
  col<-viridis(50)[findInterval((ifelse(rowSums(results[,1:time])>1,1,rowSums(results[,1:time])))*ash[prez[1:1799,1],20], bins)+1]
  plot(transform_usa, lwd=0.5, main="Exposed by 2045",col="grey")
  points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5)
  points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5, col=col)
  #plot(transform_states, lwd=0.5, main=years[time], fill=FALSE, add=T,border='white')
  par(xpd=T)
}
par(xpd=FALSE)
par(mai=c(0.5,0.25,0,0.25))
image(1:50,1,matrix(1:50), col=viridis(50), axes=FALSE, ann=F)
axis(1,labels=c("0",expression(10^{2}), expression(10^{4}),expression(10^{6})),at=c(seq(0,50,length.out=4)), cex.axis=1)
mtext(side=1, "Street Ash", line=2)
dev.off()
