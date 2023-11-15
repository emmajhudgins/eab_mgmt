rm(list=ls()) 
n_spp=66
library(sp)
library(ggplot2)
library(viridis)

setwd('../data/')

#Read in Data
data<-read.csv('data_minus5_july.csv', stringsAsFactors = FALSE)
data2<-read.csv('datanorm.csv', stringsAsFactors = FALSE)
prez2<-prez<-matrix(0,3372, 64)
prez2[,1]<-c(which(data[,paste(data2$COLUMN_NAM[1])]>0), rep(0,(3372-length(which(data[,paste(data2$COLUMN_NAM[1])]>0))))) #determine where species are present
host.density2<-read.csv('hostden_clean_gdk.csv') # tree host density by pest species
new_eab<-read.csv('eab_each_2020.csv')
L<-rep(0,64) # size of each pest's host range
ic=F
V_i<-read.csv('streettrees_grid.csv')[,20]
forecast=F
prez<-read.csv('prez_clean_gdk.csv')
L<-length(unique(c(which(V_i>0), prez[,1], prez2[,1])))-1
prez[,1]<-c(unique(c(which(V_i!=0), prez[which(prez[,1]!=0),1], prez2[which(prez2[,1]!=0),1])), rep(0, 3372-L))

currpopden<-as.matrix(read.csv("currpopden_5.csv", stringsAsFactors = FALSE))
currpopden2<-as.matrix(read.csv("future_scaled_pop2.csv"))
twenty05<-rowMeans(cbind(currpopden[,47], currpopden2[,1]))
twenty15<-rowMeans(cbind(currpopden2[,1], currpopden2[,2]))
twenty25<-rowMeans(cbind(currpopden2[,2], currpopden2[,3]))
currpopden<-cbind(currpopden, twenty05, currpopden2[,1], twenty15, currpopden2[,2], twenty25, currpopden2[,3])

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


time<-seq(1:5)
quar<-c(0.3,0.6,0.9)
bio<-c(0.1,0.3,0.5)

qin=qout=0.9
qbio=0.5
frac_site=0.2
frac_spread=0.8
results<-read.csv(paste0('../output/M_',qin,'_',qbio,'.csv'), header=F)

pdf(paste("../plots/NEW_eab_mgmt_site_street",qin, qbio,".pdf", sep="_"))
layout(matrix(c(1,2,3,4,5,6,7,7,7), ncol=3, byrow=T), heights=c(0.5,0.5,0.4))
par(mai=c(0.25,0,0.25,0))
years<-seq(2025,2045, by=5)
for (time in c(1,3,5))
{
  plot(transform_usa, lwd=0.5, main=years[time], col="white")
  points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5, col="grey")
  plot(transform_usa, lwd=0.5, main=years[time], col=alpha("grey",0), add=T)
  
  lines(quar_bound, col=alpha("deeppink4",1), lwd=2)
  points(cbind(data$X_coord[prez[which(results[((L+1):(2*L)), time]==1),1]], data$Y_coord[prez[which(results[((L+1):(2*L)), time]==1),1]]), pch=15, cex=0.5, col="yellow")
  points(cbind(data$X_coord[prez[which(results[((2*L+1):(3*L)), time]==1),1]], data$Y_coord[prez[which(results[((2*L+1):(3*L)), time]==1),1]]), pch=15, cex=0.5, col="orange")
  points(cbind(data$X_coord[prez[which(results[((3*L+1):(4*L)), time]==1),1]], data$Y_coord[prez[which(results[((3*L+1):(4*L)), time]==1),1]]), pch=15, cex=0.5, col="red")

  }
results<-read.csv(paste0('../output/management_test_0.2_0.8_',qin,'_',qbio,'_.csv'))
for (time in c(1,3,5))
{
  plot(transform_usa, lwd=0.5, main=years[time], col="white")
  points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5, col="grey")
  lines(quar_bound, col=alpha("deeppink4",1), lwd=2)
  plot(transform_usa, lwd=0.5, main=years[time], col=alpha("grey",0), add=T)
  points(cbind(data$X_coord[prez[which(results[((L+1):(2*L)), time]==1),1]], data$Y_coord[prez[which(results[((L+1):(2*L)), time]==1),1]]), pch=15, cex=0.5, col="yellow")
  points(cbind(data$X_coord[prez[which(results[((2*L+1):(3*L)), time]==1),1]], data$Y_coord[prez[which(results[((2*L+1):(3*L)), time]==1),1]]), pch=15, cex=0.5, col="orange")
  points(cbind(data$X_coord[prez[which(results[((3*L+1):(4*L)), time]==1),1]], data$Y_coord[prez[which(results[((3*L+1):(4*L)), time]==1),1]]), pch=15, cex=0.5, col="red")

  }
plot.new()
par(xpd=FALSE)
legend('center', c("Quarantine In", "Quarantine Out", "Biocontrol", "Previous Quarantine Boundary"), col=c('yellow', 'orange', 'red', "deeppink4"), pch=15, cex=2)
dev.off()



time<-seq(1:5)
quar<-c(0.3,0.6,0.9)
bio<-c(0.1,0.3,0.5)

qin=qout=0.9
qbio=0.5
frac_site=0.2
frac_spread=0.8
results<-read.csv(paste0('../output/M_',qin,'_',qbio,'.csv'), header=F)
pdf(paste("../plots/NEW_eab_mgmt_site_street_mult",qin, qbio,".pdf", sep="_"))
layout(matrix(c(1,2,3,4,5,6,7,7,7), ncol=3, byrow=T), heights=c(0.5,0.5,0.4))
par(mai=c(0.25,0,0.25,0))
years<-seq(2025,2045, by=5)
for (time in c(1,3,5))
{
  plot(transform_usa, lwd=0.5, main=years[time], col="white")
  points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5, col="grey")
  plot(transform_usa, lwd=0.5, main=years[time], col=alpha("grey",0), add=T)
  
  lines(quar_bound, col=alpha("deeppink4",1), lwd=2)
  points(cbind(data$X_coord[prez[which(results[((L+1):(2*L)), time]==1),1]], data$Y_coord[prez[which(results[((L+1):(2*L)), time]==1),1]]), pch=15, cex=0.5, col="yellow")
  points(cbind(data$X_coord[prez[which(results[((2*L+1):(3*L)), time]==1),1]], data$Y_coord[prez[which(results[((2*L+1):(3*L)), time]==1),1]]), pch=15, cex=0.5, col="orange")
  points(cbind(data$X_coord[prez[which(results[((3*L+1):(4*L)), time]==1),1]], data$Y_coord[prez[which(results[((3*L+1):(4*L)), time]==1),1]]), pch=15, cex=0.5, col="red")
  
}
results<-read.csv(paste0('../output/M_',qin,'_',qbio,'_mult.csv'), header=F)

#results<-read.csv(paste0('../output/management_test_0.2_0.8_',qin,'_',qbio,'_.csv'))
for (time in c(1,3,5))
{
  plot(transform_usa, lwd=0.5, main=years[time], col="white")
  points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5, col="grey")
  lines(quar_bound, col=alpha("deeppink4",1), lwd=2)
  plot(transform_usa, lwd=0.5, main=years[time], col=alpha("grey",0), add=T)
  points(cbind(data$X_coord[prez[which(results[((L+1):(2*L)), time]==1),1]], data$Y_coord[prez[which(results[((L+1):(2*L)), time]==1),1]]), pch=15, cex=0.5, col="yellow")
  points(cbind(data$X_coord[prez[which(results[((2*L+1):(3*L)), time]==1),1]], data$Y_coord[prez[which(results[((2*L+1):(3*L)), time]==1),1]]), pch=15, cex=0.5, col="orange")
  points(cbind(data$X_coord[prez[which(results[((3*L+1):(4*L)), time]==1),1]], data$Y_coord[prez[which(results[((3*L+1):(4*L)), time]==1),1]]), pch=15, cex=0.5, col="red")
  points(cbind(data$X_coord[prez[which(results[((3*L+1):(4*L)), time]==1 & results[((2*L+1):(3*L)), time]==1),1]], data$Y_coord[prez[which(results[((3*L+1):(4*L)), time]==1 & results[((2*L+1):(3*L)), time]==1),1]]), pch=15, cex=0.5, col="blue")
  points(cbind(data$X_coord[prez[which(results[((L+1):(2*L)), time]==1&results[((2*L+1):(3*L)), time]==1),1]], data$Y_coord[prez[which(results[((L+1):(2*L)), time]==1&results[((2*L+1):(3*L)), time]==1),1]]), pch=15, cex=0.5, col="purple")
  points(cbind(data$X_coord[prez[which(results[((L+1):(2*L)), time]==1&results[((3*L+1):(4*L)), time]==1),1]], data$Y_coord[prez[which(results[((L+1):(2*L)), time]==1&results[((3*L+1):(4*L)), time]==1),1]]), pch=15, cex=0.5, col="magenta")
  points(cbind(data$X_coord[prez[which(results[((L+1):(2*L)), time]==1&results[((3*L+1):(4*L)), time]==1& results[((2*L+1):(3*L)), time]==1),1]], data$Y_coord[prez[which(results[((L+1):(2*L)), time]==1&results[((3*L+1):(4*L)), time]==1&results[((2*L+1):(3*L)), time]==1),1]]), pch=15, cex=0.5, col="darkred")
  
}
plot.new()
par(xpd=FALSE)
legend('center', c("Quarantine In", "Quarantine Out", "Biocontrol","Quarantine out and Biocontrol", "Quarantine in and Quarantine out", "Quarantine in and Biocontrol", "All 3 Actions","Previous Quarantine Boundary"), col=c('yellow', 'orange', 'red','blue','purple','magenta', 'darkred',"deeppink4"), pch=15, cex=1.25)
dev.off()

qin=qout=0.3
qbio=0.1
frac_site=0.2
frac_spread=0.8
results<-read.csv(paste('../output/vecptime_0.3_0.1.csv'), header=F)/1000
time<-seq(1:6)
years<-seq(2020,2045, by=5)
pdf(paste("../plots/NEW_RoT_eab_dens_",qin, qout, qbio,".pdf", sep="_"))
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
results<-read.csv(paste0("../output/RoT_vecptime_", frac_site,"_",frac_spread, "_", qin,'_',qbio,'_.csv'))
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


results<-read.csv(paste('../output/vecptime_0.3_0.1.csv'), header=F)/1000
pdf(paste("../plots/NEW_eab_exp",qin,qbio,".pdf", sep="_"))
layout(matrix(c(1,2,3,4,5,6,7,7,7), ncol=3, byrow=T), heights=c(0.5,0.5,0.1))
par(mai=c(0.25,0,0.25,0))
years<-seq(2020,2045, by=5)
bins<-seq(-5,5,length.out=50)
bins<-10^bins
ash<-read.csv('streettrees_grid.csv')
 for (time in c(2,4,6))
  {
   col<-viridis(50)[findInterval(results[,time]*ash[prez[1:L,1],20], bins)+1]
   plot(transform_usa, lwd=0.5, main=years[time], col="grey")
    points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5)
    points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5, col=col)
    #plot(transform_states, lwd=0.5, main=years[time], fill=FALSE, add=T,border='white')
    par(xpd=T)
  }
  results<-read.csv(paste0("../output/RoT_vecptime_", frac_site,"_",frac_spread, "_", qin,'_',qbio,'_.csv'))
  for (time in c(2,4,6))
  {
    col<-viridis(50)[findInterval(results[,time]*ash[prez[1:L,1],20], bins)+1]
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

layout(matrix(c(1,2), ncol=1, byrow=TRUE), heights=c(0.75,0.25))
bins<-seq(0,6,length.out=50)
bins<-10^bins
col<-viridis(50)[findInterval(ash[prez[1:L,1],20], bins)+1]
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
ash<-read.csv('hostvol_sdp2020.csv')
col<-viridis(50)[findInterval(ash[prez[1:L,1],1], bins)+1]
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

dat<-read.csv('../output/postdocdat_street_debug.csv')

qz<-c(0.3,0.6,0.9)
bios<-c(0.1,0.3,0.5)
budget_scen<-data.frame(site_bud=seq(0,1, length.out=11), spread_bud=seq(1,0,length.out=11))
obj<-data.frame(frac_site=0,frac_spread=0,q_in=0,q_out=0,qbio=0, obj=0)
names(obj)
V_i<-read.csv('streettrees_grid.csv')[,20]

for (q_in in qz)
{
q_out=q_in
    for (qbio in bios)
    {
      for (scen in c(1,3,5,7,9,11))
      {
        frac_site=budget_scen$site_bud[scen]
        frac_spread=budget_scen$spread_bud[scen]
        d4prime<-read.csv(paste("../output/RoT_vecptime",frac_site,frac_spread,q_in,qbio, ".csv", sep="_"))[,2:8]
        obj<-rbind(obj, setNames(c(frac_site,frac_spread,q_in,q_out,qbio,sum(sweep(as.matrix(d4prime),MARGIN=1,as.vector(V_i[prez[,1]]+1),"*"))),names(obj)))
      }
    }
  }

obj<-obj[2:nrow(obj),]

layout(matrix(c(1,2), ncol=2, byrow=TRUE), heights=c(1,1), widths=c(0.8,0.3))
par(mar=c(4,4,2,0))
par(oma=c(0,0,0,0))
par(xpd = FALSE)

plot(y=obj$obj,x=(obj$frac_site), col=alpha(viridis(9)[as.factor(paste0(obj$q_in, obj$qbio))],0.75),xlab="Biological control proportion", ylab="Exposed ash street trees", pch=19, xlim=c(-0.05,1.05), ylim=c(100000,1500000), axes=F)
axis(1,  at=c(0,0.2,0.4,0.6,0.8,1.0))
axis(2, labels=c('0',"200000", "400000", '600000', '800000', '1000000', '1200000', '1400000', '1600000'), at=c(0,200000,400000,600000,800000, 1000000,1200000, 1400000,1600000), outer=F, las=2, cex.axis=0.75, hadj=0.8)

points(y=dat$Exposed.Ash.Street.Trees..thousands.*1000, x=dat$bio_prop/100, lty=2,col=viridis(9), pch=5)
par(mai=c(0.4,0.5,0.1,0.6))

image(1,1:9,t(matrix(1:9)), col=viridis(9), axes=FALSE, ann=F)
axis(2,labels=c("30%","30%","30%","60%","60%","60%", "90%","90%","90%"),at=c(1:9), cex.axis=0.75, padj=1)
axis(4,at=c(1:9), labels=rep("", 9),cex.axis=0.75, padj=2)
par(xpd = TRUE) #Draw outside plot area
corners<-par("usr")
text("50%      30%     10%",y=c(7.9), x=corners[2]+.2, cex=0.75, srt=270)
text("50%      30%     10%",y=c(4.9), x=corners[2]+.2, cex=0.75, srt=270)
text("50%      30%     10%",y=c(1.9), x=corners[2]+.2, cex=0.75, srt=270)

mtext(side=2, "Quarantine efficiency",line=1.25, cex=0.75)
corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
text(x = corners[2]+0.3, y = mean(corners[3:4]),"Biological Control efficiency",cex=0.7, srt = 270)


layout(1,1)
par(xpd = FALSE)
par(mar=c(4,4,2,2))
plot(dat$Exposed.Ash.Street.Trees..thousands.*1000~dat$Biocontrol.efficiency,col=viridis(27)[c(1,1,1,13,13,13,25,25,25,25)], pch=19, xlab="Biocontrol Efficiency", ylab="Exposed ash street trees", ylim=c(100000,1500000), axes=F) 
axis(1,  at=c(0,0.1,0.3,0.5), labels=c("","10%","30%","50%"))
axis(2, labels=c('0',"200000", "400000", '600000', '800000', '1000000', '1200000'), at=c(0,200000,400000,600000,800000, 1000000,1200000), outer=F, las=2,  cex.axis=0.75, hadj=0.8)

points(x=c(0.1,0.3,0.5), y=c(dat$only.biocontrol[1:3]*1000),col="red", pch=24) 
legend('topright', legend=c("30% Quarantine Efficiency","60% Quarantine Efficiency", "90% Quarantine Efficiency", "Biocontrol Only"), pch=c(19,19,19,24), col=c(viridis(27)[c(1,13,25)], "red"), cex=0.75)



qz<-c(0.3,0.6,0.9)
bios<-c(0.1,0.3,0.5)
budget_scen<-data.frame(site_bud=seq(0,1, length.out=6), spread_bud=seq(1,0,length.out=6))
obj<-data.frame(frac_site=0,q_in=0,qbio=0, obj=0)
names(obj)
V_i<-read.csv('streettrees_grid.csv')[,20]

for (q_in in qz)
{
    for (qbio in bios)
    {
      for (scen in 1:6)
      {
        frac_site=budget_scen$site_bud[scen]
        frac_spread=budget_scen$spread_bud[scen]
    
        d4prime<-read.csv(paste("../output/vecptime",frac_spread,q_in,qbio, "bud.csv", sep="_"), header=F)[,2:8]/1000
        obj<-rbind(obj, setNames(c(frac_site, q_in,qbio,sum(sweep(as.matrix(d4prime),MARGIN=1,as.vector(V_i[prez[,1]]+1),"*"))),names(obj)))
      }
    }
  }

obj<-obj[2:nrow(obj),]

layout(matrix(c(1,2), ncol=2, byrow=TRUE), heights=c(1,1), widths=c(0.7,0.3))
par(mar=c(4,4,2,0))
par(oma=c(0,0,0,0))
par(xpd = FALSE)

plot(y=obj$obj,x=(obj$frac_site), col=alpha(viridis(9)[as.factor(paste0(obj$q_in, obj$qbio))],0.75),xlab="Biological control proportion", ylab="Exposed ash street trees", pch=19, xlim=c(-0.05,1.05), ylim=c(100000,1300000), axes=F)
axis(1,  at=c(0,0.2,0.4,0.6,0.8,1.0))
axis(2, labels=c('0',"200000", "400000", '600000', '800000', '1000000', '1200000', '1400000'), at=c(0,200000,400000,600000,800000, 1000000,1200000, 1400000), outer=F, las=2, cex.axis=0.75, hadj=0.8)

points(y=dat$Exposed.Ash.Street.Trees..thousands.*1000, x=dat$bio_prop/100, lty=2,col=viridis(9), pch=5)
par(mai=c(0.4,0.5,0.1,0.6))

image(1,1:9,t(matrix(1:9)), col=viridis(9), axes=FALSE, ann=F)
axis(2,labels=c("30%","30%","30%","60%","60%","60%", "90%","90%","90%"),at=c(1:9), cex.axis=0.75, padj=1)
axis(4,at=c(1:9), labels=rep("", 9),cex.axis=0.75, padj=2)
par(xpd = TRUE) #Draw outside plot area
corners<-par("usr")
text("50%      30%     10%",y=c(7.9), x=corners[2]+.2, cex=0.75, srt=270)
text("50%      30%     10%",y=c(4.9), x=corners[2]+.2, cex=0.75, srt=270)
text("50%      30%     10%",y=c(1.9), x=corners[2]+.2, cex=0.75, srt=270)

mtext(side=2, "Quarantine efficiency",line=1.25, cex=0.75)
corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
text(x = corners[2]+0.3, y = mean(corners[3:4]),"Biological Control efficiency",cex=0.7, srt = 270)

results<-read.csv(paste('../output/vecptime_0.3_0.1.csv'))
pdf(paste("../plots/ash",qin,qbio,".pdf", sep="_"))
layout(matrix(c(1,2,3,3), ncol=2, byrow=T),heights=c(0.8,0.2))
par(mai=c(0.25,0,0.25,0))
bins<-seq(-5,5,length.out=50)
bins<-10^bins
ash<-read.csv('streettrees_grid.csv')
for (time in c(6))
{
  col<-viridis(50)[findInterval(ash[prez[1:L,1],20], bins)+1]
  plot(transform_usa, lwd=0.5, main="Total", col="grey")
  points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5)
  points(cbind(data$X_coord[prez[,1]], data$Y_coord[prez[,1]]), pch=15, cex=0.5, col=col)
  #plot(transform_states, lwd=0.5, main=years[time], fill=FALSE, add=T,border='white')
  par(xpd=T)
}
results<-read.csv("../output/vecPtime_noaction.csv")
for (time in c(6))
{
  col<-viridis(50)[findInterval((ifelse(rowSums(results[,1:time])>1,1,rowSums(results[,1:time])))*ash[prez[1:L,1],20], bins)+1]
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
