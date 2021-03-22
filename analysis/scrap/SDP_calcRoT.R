## rewrite to get a suite of site vs. spread-focussed heuristic solutions
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
  V_i<-read.csv('streettrees_grid.csv')[,20]
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
  for (sppp in 1:64)
  {L2[sppp]<-length(which(prez2[,sppp]!=0))}
  for (sppp in 1:64)
  {L3[sppp]<-length(which(gen[,sppp]!=0))}
  Tr1<-function(x)
  {
    sqrt((data$X_coord-data$X_coord[x])^2+(data$Y_coord-data$Y_coord[x])^2)
  }
  dists<-sapply(1:3372, Tr1)
  T1<-exp(-dists/50000)
  rm(dists)
  YEARS<-data2$YEAR
  pp_last<-matrix(0,3372,64) 
  risk_pp<-matrix(0,3372,3372)
  spp=1  
      par<-rep(0,23)
      startpt<-read.csv('startpt.csv')[,1]
      if (spp %in% startpt==F)
      {
        par[1]<-read.table(paste("./new_adj_inter/par",spp,"adj_inter",spp, sep="."))[,1]
      }
      if (spp %in% startpt==T)
      {
        par[1]<-read.csv(paste("./adj_inter_start/par_GDKic",spp,"csv", sep="."))[,2]
      }
      par[c(1,21,22,4,18,20,8)]<-as.numeric(c(par[1],c(0.000538410692229749, 0.299034706404549, -0.525670755351726, 15.6132848183217,-0.163552592351765, 0.323831382884772)))
      par[22]<-abs(par[22])+1

      YEAR<-YEARS[spp]
     #Pest Parameters
    total_time<-YEAR/5 +1
    Pfull<-matrix(0, 3372, total_time+5)
    Pfull_good<-matrix(0, 3372, total_time+5)
    Pfull_time<-Pfull
    vecP_time=d2prime=d3prime=d4prime=dprime=d_out=matrix(0,L[spp], total_time+5)
fitting_yr=2
        constpD=rep(0,64)
        constpD=matrix(rep(constpD),3372,64, byrow=TRUE)
        constpD2<-matrix(rep(par[9]*data[,18]+par[10]*data[,16]+par[16]*data[,19]+par[17]*data[,20]+par[18]*data[,21]),3372,64)+par[19]*host.density2
        constpD<-as.numeric(constpD)+constpD2
        constpD3<-matrix(rep(par[4]*data[,19]+par[12]*data[,20]+par[3]*data[,21]+par[13]*data[,18]+par[15]*data[,16]),3372,64)+par[5]*host.density2
          
          #Pest Parameters
          Psource=sources[[spp]]
          Discovery<-2009-YEAR
          #Pfull<-matrix(0, 3372, total_time)
          rem<-mfive(Discovery)
          
          T2<-T1[prez[1:L[spp],spp],prez[1:L[spp],spp]]
          vecP<-rep(0,L[spp])
          for (rrr in 1:length(Psource))
          {vecP[which(prez[,spp]==Psource[rrr])]=1}
  
          r0<-par[22]
          bc_pp_out<-bc_pp_in<-pp_bio<-pp_erad<-matrix(0,1799,total_time+5)
          qbio=0.1
          qerad=0.99
          q_in=q_out=0.6
          B=963943
          frac_spread=0.5
          frac_site=0.5
          c_4=c_5=c_6=c_7=matrix(0,L[spp], ncol(vecP_time))
          mgmt<-list()
          for (time in 1:(total_time+5))
          {
            
            vecP[which(prez[,spp]==Psource)]=1
            dprime[,time]<-vecP
            d2prime[,time]<-vecP
            d_out[,time]<-d2prime[,time]
            d2prime[which(vecP<par[21]),time]<-0
            c_4[which(vecP>=par[21]),time]<-1
            
            Pnext<-rep(0,L[spp])
            qq<-0
            column<-(((Discovery+5*(time-1))-1790)/5)+1

            qq<-matrix(rep(constpD[prez[which(vecP>=par[21]),spp],spp]+par[8]*currpopden[prez[which(vecP>=par[21]),spp],column], L[spp]), nrow=length(which(vecP>=par[21])), ncol=L[spp])
            zzz<-matrix(rep(constpD3[prez[1:L[spp],spp],spp]+par[20]*currpopden[prez[1:L[spp],spp],column], L[spp]), nrow=L[spp], ncol=L[spp], byrow=TRUE)
            qq<-(2*par[1]*exp((zzz[which(vecP>=par[21]),]+qq)))/(1+exp(zzz[which(vecP>=par[21]),]+qq))
            qq<-T2[which(vecP>=par[21]),]^qq
            if (length(which(vecP>=par[21]))>1){qq<-qq/rowSums(qq)}
            if (length(which(vecP>=par[21]))==1){qq<-qq/sum(qq)}
            qq[which(qq<0.001)]=0
            
            if (time<=floor(total_time))
            {
            Pnext=(vecP[which(vecP>=par[21])])%*%(qq)
            qq2<-matrix(0,L[spp], L[spp])
            qq2[which(vecP>=par[21]),]<-qq
 #           write.csv(qq2, file=paste("transmatM_", spp,time, ".csv", sep=""), row.names=F)
                        Pnext[which(prez[,spp]==Psource)]=1
            Pfull[,time]<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) 
            Pfull_time[,time]<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) #threshold at 'discoverable' value
            d3prime[,time]<-Pnext
            
            }
            if (time>floor(total_time))
            {
            vecP[which(prez[,spp]==Psource)]=1
            qq2<-matrix(0,L[spp], L[spp])
            qq2[which(vecP>=par[21]),]<-qq
            diag(qq2)<-0
            pp<-sweep(qq2,1,vecP,'*')
            Pnext<-rep(0,L[spp])
            
            
            shortpath_in<-apply(sweep(qq2,1,vecP,'*'),1,which.max)
            shortpath_out<-apply(sweep(qq2,1,vecP,'*'),2,which.max)
            for (i in 1:1799)
            {
              # bc_sources[i]<-length(which(pp[i,abs]>par[21]))
              # bc_destinations[i]<-length(which(pp[abs,i]>par[21]))
              bc_pp_out[i,time]<-sum(pp[i,which(shortpath_out==i)]*V_i[which(shortpath_out==i)])
              bc_pp_in[i,time]<-sum(pp[which(shortpath_in==i),i]*V_i[prez[i,1]])
              pp_bio[i,time]<-(vecP[i]*V_i[prez[i,1]]*(1-qbio))/141519
              pp_erad[i,time]<-(vecP[i]*V_i[prez[i,1]]*(1-qerad))/(31000*2500*vecP[i])
             # bc[i]<-length(which(shortpath_out==i))+length(which(shortpath_in==i))
            }
            pp_erad[which((31000*2500*vecP)>B),time]<-0
            pp_erad[which(vecP<par[21])]<-0
            mgmt[[time]]<-vector()
            
            ce_spread<-rank(c(bc_pp_in[,time],bc_pp_out[,time]), ties.method="random")
            ce_spread[which(ce_spread==min(ce_spread))]<-NA
            tt=min(ce_spread,na.rm=T)
            cost3=cost4=0
            while(cost4<=B*frac_spread & sum(ce_spread[which(ce_spread>=tt)], na.rm=T)!=0)
            {
              if(which(ce_spread==tt)%in%c(mgmt[[time]], mgmt[[time]]-1799,mgmt[[time]]+1799 )==F)
              {mgmt[[time]]<-c(mgmt[[time]],which(ce_spread==tt))
              cost3=cost3+(646863/309)}
              tt=min(ce_spread[which(ce_spread>tt)],na.rm=T)
              cost4=cost3+(646863/309)
            }
            ce_site<-rank(c(pp_bio[,time],pp_erad[,time]),ties.method="random")
            cost_site<-c(rep(141519,1799),(31000*2500*vecP) )
            ce_site[which(ce_site==min(ce_site))]<-NA
            ce_site[which(cost_site>B*frac_site)]<-NA
            tt=min(ce_site,na.rm=T)
            cost=cost2=0
            while(cost2<=B*frac_site& sum(ce_site[which(ce_site>=tt)], na.rm=T)!=0)
            {
              if (which(ce_site==tt)%in%c(mgmt[[time]], mgmt[[time]]-1799,mgmt[[time]]-(2*1799),mgmt[[time]]-(3*1799),mgmt[[time]]+1799,mgmt[[time]]+(2*1799),mgmt[[time]]+(3*1799))==F)
              {mgmt[[time]]<-c(mgmt[[time]],which(ce_site==tt)+2*1799)
              cost=cost+cost_site[which(ce_site==tt)]}
              tt=min(ce_site[which(ce_site>tt)],na.rm=T)
              cost2=cost+cost_site[which(ce_site==tt)]
            }
            diag(qq2[which(vecP>=par[21]),])<-diag(qq)
            vecP[(subset(mgmt[[time]], mgmt[[time]]>(3*1799) & mgmt[[time]]<=(4*1799))-(3*1799))]<-(1-qerad)*vecP[(subset(mgmt[[time]], mgmt[[time]]>(3*1799) & mgmt[[time]]<=(4*1799))-(3*1799))]
            vecP[(subset(mgmt[[time]], mgmt[[time]]>2*1799 & mgmt[[time]]<=(3*1799))-(2*1799))]<-(1-qbio)*vecP[(subset(mgmt[[time]], mgmt[[time]]>(2*1799) & mgmt[[time]]<=(3*1799))-(2*1799))]
            
            dprime[,time]<-vecP
            d2prime[,time]<-vecP
            d2prime[which(vecP<par[21]),time]<-0
            c_4[which(vecP>=par[21]),time]<-1
             qq2[(subset(mgmt[[time]], mgmt[[time]]>(1799) & mgmt[[time]]<=(2*1799))-1799),]<-(1-q_out)*qq2[(subset(mgmt[[time]], mgmt[[time]]>(1799) & mgmt[[time]]<=(2*1799))-1799),]
              d_out[,time]<-d2prime[,time]
             d_out[(subset(mgmt[[time]], mgmt[[time]]>(1799) & mgmt[[time]]<=(2*1799))-1799),time]<-vecP[(subset(mgmt[[time]], mgmt[[time]]>(1799) & mgmt[[time]]<=(2*1799))-1799)]*(1-q_out)
             qq2[,(subset(mgmt[[time]], mgmt[[time]]<(1799)))]<-(1-q_in)*qq2[,(subset(mgmt[[time]], mgmt[[time]]<(1799)))]
             diag(qq2[which(vecP>=par[21]),])<-diag(qq)
             qq<-qq2[which(vecP>=par[21]),]
            Pnext=(vecP[which(vecP>=par[21])])%*%qq
            qq2<-matrix(0,L[spp], L[spp])
            qq2[which(vecP>=par[21]),]<-qq
            write.csv(qq2, file=paste("transmatM_", spp,time, ".csv", sep=""), row.names=F)
            d3prime[,time]<-Pnext
            
            Pnext[which(prez[,spp]==Psource)]=1
            Pfull[,time]<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) 
            Pfull_time[,time]<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) #threshold at 'discoverable' value
            reset_yrs<-YEAR/5
              
            }
            if (time>1)
            {

              dddd<-which(!(Pfull_time[1:length(which(Pfull_time[,time-1]!=0)),time-1]%in%Pfull_time[1:length(which(Pfull_time[,time]!=0)),time]))
              ffff<-which(prez[1:length(which(prez[,spp]!=0)), spp]%in%Pfull_time[dddd,time-1])
              if (time>floor(total_time))
              {ffff<-ffff[-which(ffff%in%c(mgmt[[time]]-1799,mgmt[[time]]-2*1799,mgmt[[time]]-3*1799,mgmt[[time]]))]}
        
              Pnext[ffff]<-par[21]
              Pfull[,time]<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) 
              Pfull_time[,time]<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) #threshold at 'discoverable' value
            }
            if (time==floor(YEAR/5))
              {
                dddd<-which(prez[1:length(which(prez[,spp]!=0)),spp]%in%prez2[1:length(which(prez2[,spp]!=0)),spp])
                cccc<-which(!(prez[1:length(which(prez[,spp]!=0)),spp]%in%prez2[1:length(which(prez2[,spp]!=0)),spp]))
                eeee<-which(Pnext[dddd]<par[21])
                Pnext[dddd[eeee]]<-par[21] #testing out setting to true presences
                Pnext[cccc]<-0 #remove this when not doing sdp
                Pfull_time[,time]<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) 
            }
            
            Pnext[which(prez[,spp]==Psource)]=1
            Pnext[which(Pnext>=par[21])]=Pnext[which(Pnext>=par[21])]*r0
            d4prime[,time]<-Pnext
            c_5[which(d4prime[,time]>=par[21]),time]=1
            c_6[which(d4prime[,time]<1),time]=1
            c_7[which(d4prime[,time]>=1),time]=1
            Pnext[which(Pnext>=1)]<-1
  
            vecP=Pnext
#            vecP[which(prez[,spp]==Psource)]=1

              vecP_time[,time]<-vecP
           }
            
      

          write.csv(dprime[,4:9], file=paste("dprime", ".csv", sep=""), row.names=F)
          write.csv(d2prime[,4:9], file=paste("d2prime", ".csv", sep=""), row.names=F)
          write.csv(d3prime[,4:9], file=paste("d3prime", ".csv", sep=""), row.names=F)
          write.csv(d4prime[,4:9], file=paste("d4prime", ".csv", sep=""), row.names=F)
          write.csv(vecP_time[,3:9], file=paste("vecptime",".csv", sep=""), row.names=F)
          write.csv(c_4[,4:9], file=paste("c_4",".csv", sep=""), row.names=F)
          write.csv(c_5[,4:9], file=paste("c_5",".csv", sep=""), row.names=F)
          write.csv(c_6[,4:9], file=paste("c_6",".csv", sep=""), row.names=F)
          write.csv(c_7[,4:9], file=paste("c_7",".csv", sep=""), row.names=F)
          
          M_big<-matrix(0,1799*5,5)
          M_big[1:1799,]<-1
          for (i in 1:5)
          {
            M_big[mgmt[[i+4]]+1799,i]<-1
            M_big[which(1:1799%in%c(mgmt[[i+4]],mgmt[[i+4]]-1799,mgmt[[i+4]]-(2*1799),mgmt[[i+4]]-(3*1799))),i]<-0
          }
          write.csv(M_big, file="management_test.csv", row.names=F)
    write.csv(d_out[,4:9], file="d_out.csv", row.names=F)
    