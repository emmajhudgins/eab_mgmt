rm(list=ls()) 
## reference existing github data here instead of local files
#Read in Data
setwd('../../UStreedamage/data/')
data<-read.csv('countydatanorm_march.csv', stringsAsFactors = FALSE) # spatial data
data2<-read.csv('spdat_clean_gdk.csv', stringsAsFactors = FALSE) # species data, see Hudgins et al. corrigendum for information on why ALB (spp=3) cannot be accurately fit
host.density2<-read.csv('hostden_clean_gdk.csv') # tree host density by pest species
prez<-read.csv('prez_clean_gdk.csv') # invasible host range (from FIA)
prez2<-read.csv('prez2_clean_gdk.csv') # pest presences for all species
prez2[,1]<-readRDS('../output/presences_time_eab.rds')[[1]][,5] # pest presences in 2020
L<-rep(0,64) # size of each pest's host range
for (sppp in 1:64)
{
  L[sppp]<-length(which(prez[,sppp]!=0))
}
currpopden<-as.matrix(read.csv("currpopden_5.csv", stringsAsFactors = FALSE))
currpopden2<-as.matrix(read.csv("future_scaled_pop2.csv"))
twenty35<-rowMeans(cbind(currpopden[,47], currpopden2[,4]))
twenty45<-rowMeans(cbind(currpopden2[,4], currpopden2[,5]))
twenty55<-rowMeans(cbind(currpopden2[,5], currpopden2[,6]))
currpopden<-cbind(currpopden, twenty35, currpopden2[,4], twenty45, currpopden2[,5], twenty55, currpopden2[,6])
V_i<-read.csv('../../eab_mgmt/data/streettrees_grid.csv')[,20]
Tr1<-function(x)
{
  sqrt((data$X_coord-data$X_coord[x])^2+(data$Y_coord-data$Y_coord[x])^2)
}
sources<-as.list(read.csv('Psources_notypos.csv')[,1])
dists<-sapply(1:3372, Tr1)
T1<-exp(-dists/50000)
rm(dists)
YEAR<-28 #assume initial invasion in in 1992
spp=1  
par<-rep(0,23)
par[c(1,21,22,4,18,20,8)]<-as.numeric(c(2.675,c(0.000538410692229749, 0.299034706404549, -0.525670755351726, 15.6132848183217,-0.163552592351765, 0.323831382884772)))
par[22]<-abs(par[22])+1


#Pest Parameters
total_time<-YEAR/5

constpD=rep(0,64)
constpD=matrix(rep(constpD),3372,64, byrow=TRUE)
constpD2<-matrix(rep(par[9]*data[,18]+par[10]*data[,16]+par[16]*data[,19]+par[17]*data[,20]+par[18]*data[,21]),3372,64)+par[19]*host.density2
constpD<-as.numeric(constpD)+constpD2
constpD3<-matrix(rep(par[4]*data[,19]+par[12]*data[,20]+par[3]*data[,21]+par[13]*data[,18]+par[15]*data[,16]),3372,64)+par[5]*host.density2

#Pest Parameters
Psource=sources[[spp]]
Discovery<-2020-YEAR            

T2<-T1[prez[1:L[spp],spp],prez[1:L[spp],spp]]
qq3_list<-list()
for (time in 6:11)
{
  qq3_list[[time]]<-as.matrix(read.csv(paste("../../eab_mgmt/data/transmatM_", spp,time,0,1,0.3, 0.3,0.1, ".csv", sep="_")))
}
adj_list<-read.csv("../../eab_mgmt/output/adj_list.csv")
r0<-par[22]
budget_scen<-data.frame(site_bud=seq(0,1, length.out=11), spread_bud=seq(1,0,length.out=11))
qz<-c(0.3,0.6,0.9)
bios<-c(0.1,0.3,0.5)
B=963943
for (q_out in qz)
{
  for (qbio in bios)
  {
    
    q_in=q_out
    for (scen in 1:11)
    {
      vecP<-rep(0,L[spp])
      for (rrr in 1:length(Psource))
      {vecP[which(prez[,spp]==Psource[rrr])]=1}
      
      bc_pp_out<-bc_pp_in<-pp_bio<-matrix(0,1799,total_time+6)
      Pfull<-matrix(0, 3372, total_time+6)
      Pfull_good<-matrix(0, 3372, total_time+6)
      Pfull_time<-Pfull
      vecP_time=d2prime=d3prime=d4prime=dprime=d_out=matrix(0,L[spp], total_time+6)
      
      
      frac_spread=budget_scen$spread_bud[scen]
      frac_site=budget_scen$site_bud[scen]
      c_4=c_5=c_6=c_7=c_8=matrix(0,L[spp], ncol(vecP_time))
      mgmt<-list()
      for (time in 1:(total_time+6))
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
        
        if (time<=floor(total_time))
        {
          Pnext=(vecP[which(vecP>=par[21])])%*%(qq)
          qq2<-matrix(0,L[spp], L[spp])
          qq2[which(vecP>=par[21]),]<-qq
          #write.csv(qq2, file=paste("transmatM_", spp,time, ".csv", sep=""), row.names=F)
          Pnext[which(prez[,spp]==Psource)]=1
          Pnext[which(Pnext<0)]<-0
          Pfull[,time]<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) 
          Pfull_time[,time]<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) #threshold at 'discoverable' value
          d3prime[,time]<-Pnext
          
        }
        if (time>floor(total_time))
        {
          vecP[which(prez[,spp]==Psource)]=1
          qq3<-matrix(0,L[spp], L[spp]) # full transition matrix
          qq3[which(vecP>=par[21]),]<-qq
          qq3<-qq3_list[[time]]
          # write.csv(qq3, file=paste("transmatM_", spp,time,frac_site, frac_spread, q_in,q_out,qbio, ".csv", sep="_"), row.names=F)
          qq2<-qq3 #transition matrix with 0 on diagonal
          diag(qq2)<-0
          pp<-sweep(qq2,2,vecP,'*')
          Pnext<-rep(0,L[spp])
          
          shortpath_in<-apply(pp,2,which.max)
          shortpath_out<-apply(pp,1,which.max)
          for (i in 1:1799)
          {
            bc_pp_out[i,time]<-sum(length(which(shortpath_out==i))*(V_i[prez[which(shortpath_out==i),1]]+1)*(1-vecP[which(shortpath_out==i)]))
            bc_pp_in[i,time]<-sum(length(which(shortpath_in==i))*(V_i[prez[i,1]]+1)*(1-vecP[i]))
            pp_bio[i,time]<-(vecP[i]*(V_i[prez[i,1]]+1)*(1-qbio))/141519
          }
          
          mgmt[[time]]<-vector()
          ce_site<-rank(pp_bio[,time],ties.method="random")
          cost_site<-c(rep(50000,1799))
          if (length(mgmt[[time-4]][which(mgmt[[time-4]]>2*L[spp])])>0){
            for(xx in 1:length(mgmt[[time-4]][which(mgmt[[time-4]]>2*L[spp])]))
            {
              c_8[unlist(adj_list[mgmt[[time-4]][which(mgmt[[time-4]]>2*L[spp])][xx]-(2*L[spp]),]),time]<-1
            }}
          c_8[mgmt[[time-2]][which(mgmt[[time-2]]>2*L[spp])]-(2*L[spp]),time]<-1
          
          ce_site[which(pp_bio[,time]==0)]<-NA
          ce_site[which(cost_site>B*frac_site)]<-NA
          ce_site[which(c_8[,time]==1)]<-NA
          tt=min(ce_site,na.rm=T)
          cost=cost2=0
          while(is.infinite(tt)==F&cost2<=B*frac_site& sum(ce_site[which(ce_site>=tt)], na.rm=T)!=0 & tt<1799)
          {
            if (which(ce_site==tt)%in%c(mgmt[[time]], mgmt[[time]]-1799,mgmt[[time]]-(2*1799),mgmt[[time]]-(3*1799),mgmt[[time]]+1799,mgmt[[time]]+(2*1799),mgmt[[time]]+(3*1799))==F)
            {mgmt[[time]]<-c(mgmt[[time]],which(ce_site==tt)+2*1799)
            cost=cost+cost_site[which(ce_site==tt)]}
            tt=min(ce_site[which(ce_site>tt)],na.rm=T)
            if (length(tt)>0){
              while ((cost_site[which(ce_site==tt)]>(B-cost2) ))
              {
                tt=ce_site[which(ce_site>tt)][order(ce_site[which(ce_site>tt)])][2]
                if (is.na(tt))
                {
                  break
                }
                if (tt>1799)
                {
                  break
                }}
            }
            if (is.na(tt))
            {
              break
            }
            if (tt<1799)
            {
              cost2=cost+cost_site[which(ce_site==tt)]
            }
            if (tt>1799)
            {
              break
            }
          }
          ce_spread<-rank(c(bc_pp_in[,time],bc_pp_out[,time]), ties.method="random")
          ce_spread[c(which(bc_pp_in[,time]==0), which(bc_pp_out[,time]==0)+1799)]<-NA
          tt=min(ce_spread,na.rm=T)
          cost3=cost4=0
          while(is.infinite(tt)==F& cost4<=B*frac_spread & sum(ce_spread[which(ce_spread>=tt)], na.rm=T)!=0)
          {
            if(which(ce_spread==tt)%in%c(mgmt[[time]], mgmt[[time]]-1799,mgmt[[time]]-(2*1799),mgmt[[time]]-(3*1799),mgmt[[time]]+1799,mgmt[[time]]+(2*1799),mgmt[[time]]+(3*1799))==F)
            {mgmt[[time]]<-c(mgmt[[time]],which(ce_spread==tt))
            cost3=cost3+(646863/309)}
            tt=min(ce_spread[which(ce_spread>tt)],na.rm=T)
            cost4=cost3+(646863/309)
          }
          vecP[which(c_8[,time]==1)]<-(1-qbio)*vecP[which(c_8[,time]==1)]
          
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
          if (time>(floor(total_time)))
            if (any(ffff%in%c(mgmt[[time]]-1799,mgmt[[time]]-2*1799,mgmt[[time]]-3*1799,mgmt[[time]])) )
            {
              {ffff<-ffff[-which(ffff%in%c(mgmt[[time]]-1799,mgmt[[time]]-2*1799,mgmt[[time]]-3*1799,mgmt[[time]]))]}}
          
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
        c_5[which(d3prime[,time]>=par[21]),time]=1
        c_6[which(d4prime[,time]<1),time]=1
        c_7[which(d4prime[,time]>=1),time]=1
        Pnext[which(Pnext>=1)]<-1
        
        vecP=Pnext
        vecP_time[,time]<-vecP
      }
      
      write.csv(vecP_time[,5:11], file=paste("../../eab_mgmt/output/vecptime",frac_site,frac_spread,q_in,qbio,"second.csv", sep="_"), row.names=F)
      
      # M_big<-matrix(0,1799*5,5)
      # M_big[1:1799,]<-1
      # for (i in 1:5)
      # {
      #   M_big[mgmt[[i+6]]+1799,i]<-1
      #   M_big[which(1:1799%in%c(mgmt[[i+6]],mgmt[[i+6]]-1799,mgmt[[i+6]]-(2*1799))),i]<-0
      # }
      # write.csv(M_big, file=paste("../../eab_mgmt/output/management_test",frac_site,frac_spread,q_in,qbio,".csv", sep="_"), row.names=F)
      #   
    }
  }
}
