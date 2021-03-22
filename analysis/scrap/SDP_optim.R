  ## rewrite to get a suite of site vs. spread-focussed heuristic solutions

  rm(list=ls())
  library(rgenoud)
  n_spp=66
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
           
            r0<-par[22]
           budget_scen<-data.frame(site_bud=seq(0,1, length.out=11), spread_bud=seq(1,0,length.out=11))
           qz<-c(0.3,0.6,0.9)
           bios<-c(0.1,0.3,0.5)
           qerad<-0.99
           B=963943
           setwd("~/Desktop/OneDrive - McGill University/Grad/scripts/")
           
           # for (q_in in qz)
           # {
           #   for (q_out in qz)
           #   {
           #     for (qbio in bios)
           #     {
           #   
           # 
           #  for (scen in 1:11)
           #  {
           q_in=q_out=0.3
           qbio=0.5

           opt_manage<-function(par){
            par<-round(par)
           mgmt_itme<-matrix(par, nrow=1799*5, ncol=5)
           mgmt<-list()
           for (time in 6:10)
           {
             mgmt[[time]]<-which(mgmt_itme[1800:nrow(mgmt_itme),time-5]>0)
           }
           mgmt[[11]]<-0
              vecP<-rep(0,L[spp])
              for (rrr in 1:length(Psource))
              {vecP[which(prez[,spp]==Psource[rrr])]=1}
              
            Pfull<-matrix(0, 3372, total_time+7)
            Pfull_good<-matrix(0, 3372, total_time+7)
            Pfull_time<-Pfull
            vecP_time=d2prime=d3prime=d4prime=dprime=d_out=matrix(0,L[spp], total_time+7)
          

            #frac_spread=0.5
            #frac_site=0.5
            # frac_spread=budget_scen$spread_bud[scen]
            #  frac_site=budget_scen$site_bud[scen]
            c_4=c_5=c_6=c_7=matrix(0,L[spp], ncol(vecP_time))
            #mgmt<-list()
            for (time in 1:(total_time+7))
            {

              vecP[which(prez[,spp]==Psource)]=1
              dprime[,time]<-vecP
              d2prime[which(vecP<par[21]),time]<-0
              d2prime[,time]<-vecP
              d_out[,time]<-d2prime[,time]
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
              
              if (time<=floor(total_time)+1)
              {
              Pnext=(vecP[which(vecP>=par[21])])%*%(qq)
              qq2<-matrix(0,L[spp], L[spp])
              qq2[which(vecP>=par[21]),]<-qq
   #           write.csv(qq2, file=paste("transmatM_", spp,time, "site.csv", sep=""), row.names=F)
              Pnext[which(prez[,spp]==Psource)]=1
              Pnext[which(Pnext<0)]<-0
              Pfull[,time]<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) 
              Pfull_time[,time]<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) #threshold at 'discoverable' value
              d3prime[,time]<-Pnext
              
              }
              if (time>floor(total_time)+1)
              {
              vecP[which(prez[,spp]==Psource)]=1
              qq3<-as.matrix(read.csv(paste("transmatM_", spp,time,0,1, 0.3,0.3,0.1, ".csv", sep="_")))
              vecP[(subset(mgmt[[time]], mgmt[[time]]>(3*1799) & mgmt[[time]]<=(4*1799))-(3*1799))]<-(1-qerad)*vecP[(subset(mgmt[[time]], mgmt[[time]]>(3*1799) & mgmt[[time]]<=(4*1799))-(3*1799))]
              vecP[(subset(mgmt[[time]], mgmt[[time]]>2*1799 & mgmt[[time]]<=(3*1799))-(2*1799))]<-(1-qbio)*vecP[(subset(mgmt[[time]], mgmt[[time]]>(2*1799) & mgmt[[time]]<=(3*1799))-(2*1799))]
              
              dprime[,time]<-vecP
              d2prime[,time]<-vecP
              d2prime[which(vecP<par[21]),time]<-0
              c_4[which(vecP>=par[21]),time]<-1
                d_out[,time]<-d2prime[,time]
               d_out[(subset(mgmt[[time]], mgmt[[time]]>(1799) & mgmt[[time]]<=(2*1799))-1799),time]<-d2prime[(subset(mgmt[[time]], mgmt[[time]]>(1799) & mgmt[[time]]<=(2*1799))-1799),time]*(1-q_out)
               Pnext[which(1:1799%in%(mgmt[[time]]))]<-(((1-q_in)*((d_out[,time]%*%qq3)-d_out[,time]*diag(qq3)))+vecP*diag(qq3))[which(1:1799%in%(mgmt[[time]]))]
            
          Pnext[which(1:1799%in%(mgmt[[time]])==F)]<-((d_out[,time]%*%qq3)-(d_out[,time])*diag(qq3)+vecP*diag(qq3))[which(1:1799%in%(mgmt[[time]])==F)]
          
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
                if (time>(floor(total_time)+1))
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
              vecP_time[,time]<-vecP
            }
            cost_vec<-c(rep(0, 1799), rep((646863/309),1799*2), rep(141519, 1799), rep(31000*2500, 1799))
            cost_nonerad<-colSums(mgmt_itme[1:(4*1799),]*cost_vec[1:(4*1799)])
            cost_erad<-colSums(mgmt_itme[(4*1799+1):(5*1799),]*cost_vec[(4*1799+1):(5*1799)]*vecP_time[,6:10])
            if(any(cost_nonerad+cost_erad>963943)){
              return (1e+12)
            }
            return(sum(sweep(as.matrix(vecP_time)[,6:10],MARGIN=1,as.vector(V_i[prez[,1]]+1),"*")))
           }
mgmt_itme<-read.csv('./RoT/management_test_0_1_0.3_0.3_0.1_site.csv')
mgmt_itme<-round(mgmt_itme)
m<-genoud(starting.values =unlist(mgmt_itme), fn=opt_manage,data.type.int = TRUE,max=F, nvars=8995*5, Domains=cbind(rep(0,8995*5), rep(1,8995*5)), boundary.enforcement=2, solution.tolerance=1000, max.generations = 1, pop.size=1)
opt_manage(unlist(mgmt_itme))
