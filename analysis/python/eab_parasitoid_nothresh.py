import sys
import subprocess
from gurobipy import *
import numpy
import pandas
import itertools
import random
import scipy
import re


#%% Spread model data

spp = 1 #which column of pest data is EAB?

B = 1650 # yearly budget in thousands $, equivalent to number of biocontrol sites allocated 2015-2020
HostVol = pandas.io.parsers.read_csv('../../data/streettrees_grid.csv') #modelled street tree count in each cell (ash trees in column 20) from Hudgins et al. 2021
prez=pandas.io.parsers.read_csv('../../data/prez_clean_gdk.csv') #current distribution of all pests (inc. EAB in column 1)
prez=prez.iloc[:,0] #subset to EAB
HostVol = HostVol/1e+06 # convert to millions of trees
L = 1799

#previous biocontrol history 
bio_hist = pandas.io.parsers.read_csv("../../data/biocontrol_history.csv")
bio_2010 = prez.index[prez.isin(bio_hist["V3"][bio_hist["V2"] == 2010])].tolist()
bio_2015 = prez.index[prez.isin(bio_hist["V3"][bio_hist["V2"] == 2015])].tolist()
bio_2020 = prez.index[prez.isin(bio_hist["V3"][bio_hist["V2"] == 2020])].tolist()


adj_mat= pandas.io.parsers.read_csv('../../output/adj_list.csv') #neighbour matrix for use in parasitoid dispersal
adj_mat.astype(int)
adj_list= list([] for x in range(L))
for i in range(L):
    adj_row=adj_mat.iloc[i,:]
    adj_list[i].append(adj_row.where(adj_row!=0).dropna())


source = 613
time = range(1,7+1) #35 years
time_sub = range(1,5+1) #25 years

r = 1.299034706404549 # pest growth rate from published model
phi =0.538410692229749 #pest detectability/dispersal density threshold from published model

for rr in range(0,3): #iterate over efficiency scenarios
    for qq in range(0,3): 
        #efficiencies 
        effs_quar=[0.3,0.6,0.9] #efficiency scenarios for quarantines
        effs_bio=[0.1,0.3,0.5] #for biocontrol
        sites2 = list(range(1,L+1))
        sites2.remove(613) # maintain source location at maximum propagule pressure
        sites = range(1,L+1)

        eff_nothing = 0
        n_actions = 4 #includes no action
        all_comb = range(1, L * n_actions+1)
        range_none = range(1, L+1)
        range_quar_in = range(L+1, 2*L+1)
        range_quar_out= range(2*L+1, 3*L+1)
        range_bio = range(3*L+1, 4*L+1)

        #Costs
        cost_quar_in=list(itertools.repeat(15.421, L)) # aphis 2020 biocontrol budget divided by eab quarantined grid cells
        cost_quar_out = list(itertools.repeat(15.421, L))
        cost_bio =list(itertools.repeat(50, L))
        cost_nothing = list(itertools.repeat(0, L))
        cost_vec = pandas.Series(cost_nothing+ cost_quar_in+ cost_quar_out+ cost_bio)


        #set efficiencies
        eff_quar_in = effs_quar[rr]
        eff_quar_out=effs_quar[rr]
        eff_bio = effs_bio[qq]
        eff_vec = pandas.Series(list(itertools.repeat(eff_nothing, L))+ list(itertools.repeat(eff_quar_in, L))+list(itertools.repeat(eff_quar_out, L))+list(itertools.repeat(eff_bio, L)))

        #%%Create optimization model
        #Budget
        # Add variables
        m = Model('spread')
        d = m.addVars(sites,range(1, len(time)+2), vtype=GRB.CONTINUOUS, name="d",lb=0, ub=5000)  # density in each cell
        M = m.addVars(all_comb, time_sub, vtype=GRB.BINARY, name="M")# decision variables
        B_each = m.addVars(time_sub,vtype=GRB.CONTINUOUS, name="B", lb=0,ub=B) #budget each year

        dprime = m.addVars(sites,time, vtype=GRB.CONTINUOUS, name="dprime", lb=0,ub=1000) #density variables before and after different model equations/constraints
        d2prime = m.addVars(sites,time, vtype=GRB.CONTINUOUS, name="d2prime", lb=0,ub=1000)
        d3prime = m.addVars(sites,time, vtype=GRB.CONTINUOUS, name="d3prime", lb=0,ub=3000/phi) #3000/phi is a reasonable upper bound to limit size of model matrix
        d4prime = m.addVars(sites,time, vtype=GRB.CONTINUOUS, name="d4prime",lb=0,ub=3000/phi)
        d_out=m.addVars(sites, time, vtype=GRB.CONTINUOUS, name="d_out", lb=0,ub=1000)


        c_4 = m.addVars(sites,time, vtype=GRB.BINARY, name="c_4") #different indicator constraints
        c_5 = m.addVars(sites,time, vtype=GRB.BINARY, name="c_5")
        c_6 = m.addVars(sites,time, vtype=GRB.BINARY, name="c_6")
        c_7 = m.addVars(sites,time, vtype=GRB.BINARY, name="c_7")
        c_8 = m.addVars(sites,range(1,7), vtype=GRB.BINARY, name="c_8")

        #Set objective - minimize exposed street ash volume over all timesteps (exposure calculated by multiplication of relative propagule pressure and # street trees in each cell, small number is added to avoid zeroes)
        m.setObjective(quicksum(d[loc, year]*(HostVol.iloc[prez.iloc[loc-1]-1,19]+0.000001) for loc in sites for year in range(1+1,len(time)+2)), GRB.MINIMIZE)
         

        #%% Add constraints 

        # Only allow one type of management per cell
        for ii in sites:
            j = (ii, L+ii, (2*L)+ii, (3*L)+ii)
            m.addConstrs(((quicksum(M[loc,year]for loc in j) == 1) for year in time_sub), name = "One Action") #limit management actions to one per site per timestep

        # Budget constraint

        m.addConstrs((B_each[year]==quicksum(M[loc,year]*cost_vec[loc-1] for loc in range(1,4*L+1)) for year in time_sub), name = "Budget")
        m.addConstrs((B_each[year] <=B for year in time_sub), name = "Total Budget")

        #initial density taken from R simulation
        vecptime = pandas.io.parsers.read_csv("../../output/vecptime_0_1_{0}_{1}_.csv".format(eff_quar_in, eff_bio))*1000 #density ix x1000 to reduce orders of magnitude in model matrix
        m.addConstrs(((d[loc,1]==vecptime.iloc[loc-1,0]) for loc in sites2), name="initial_den")

        #apply biocontrol to reduce density prior to dispersal

        #biological control has 50% efficacy after 1 timestep and 100% efficacy after 2
        m.addConstrs(((dprime[ii,year] >= d[ii,year]-((1000/phi) *c_8[ii,year-2])-((1000/phi) *c_8[ii,year-1])) for ii in sites for year in range(3,8)), name = "nobiocontrol)")
        m.addConstrs(((dprime[ii,year] >= (1-eff_vec[3*L+ii-1])*d[ii,year]-((1000/phi) *(1-c_8[ii,year-2]))) for ii in sites for year in range(3,8)), name = "biocontrol")
        m.addConstrs(((dprime[ii,year] >= d[ii,year]-((1000/phi) *c_8[ii,year-1])-((1000/phi) *c_8[ii,year-2])) for ii in sites for year in range(3,8)), name = "nobiocontrol)")
        m.addConstrs(((dprime[ii,year] >= (1-(eff_vec[3*L+ii-1]*0.5))*d[ii,year]-((1000/phi)*(1-c_8[ii,year-1])+((1000/phi)*(c_8[ii,year-2])))) for ii in sites for year in range(3,8)), name = "biocontrol") # 50% impact after 1 timestep, 100% after 2
        m.addConstrs(((dprime[ii,year] >= d[ii,year]-((1000/phi) *c_8[ii,year-1])) for ii in sites for year in range(2,3)), name = "nobiocontrol)")
        m.addConstrs(((dprime[ii,year] >= (1-(eff_vec[3*L+ii-1]*0.5))*d[ii,year]-((1000/phi)*(1-c_8[ii,year-1]))) for ii in sites for year in range(2,3)), name = "biocontrol") # 50% impact after 1 timestep, 100% after 2

        # add impact of biocontrol (c_8) in focal site and in adjacent sites to focal site after a 2-timestep time lag
        m.addConstrs(((c_8[ii,year]>=M[3*L+ii,year]) for ii in sites for year in range(2,5)), name = "biocontrol2")
        m.addConstrs(((c_8[ii,year]>=M[3*L+jj,year-2]) for jj in adj_list[ii-1][0].astype(int) for year in range(4,5)), name = "parasitoidDisp")

        # take into account historical biocontrol prior to today (because of time lags)
        start= list([x+1 for x in bio_2020])
        for j in range(0, len(bio_2015)):
            start.extend(x for x in [adj_list[q][0] for q in bio_2015][j] if x-1 not in bio_2020)
        nostart= list(x for x in sites if x not in start)
        second = list()
        for j in range(0, len(bio_2020)):
            second.extend(x for x in [adj_list[q][0] for q in bio_2020][j] if x-1 not in bio_2015)
        no_2= list(x for x in sites if x not in second)
#            third = list(bio_2010)
#            no_3= list(x for x in sites if x not in third)

            #c_8 is ==0 if sites not in biocontrol history or turned on by present-day management
        m.addConstrs(((c_8[ii,year]<=M[3*L+ii,year]) for ii in nostart for year in range(1,2)), name = "biocontrol2")
   
        m.addConstrs(((c_8[ii,year]<=M[3*L+ii,year]) for ii in no_2 for year in range(2,3)), name = "biocontrol2")
           # m.addConstrs(((c_8[ii,year]<=M[3*L+ii,year]) for ii in no_3 for year in range(3,4)), name = "biocontrol2")
        m.addConstrs(((c_8[ii,year]<=M[3*L+ii,year]+quicksum(M[3*L+jj,year-2] for jj in adj_list[ii-1][0].astype(int))) for ii in sites for year in range(3,6)), name = "biocontrol2")
        #m.addConstrs(((c_8[ii,year]<=M[3*L+ii,year]) for ii in sites for year in range(3,6)), name = "biocontrol2")
        m.addConstrs(((c_8[ii,year]==0) for ii in sites for year in range(6,7)), name = "biocontrol2")

#            ## sites not adjacent to biological control applied in 2020 (used below)
#            not_adj_bio_2020 = list(range(1,L+1))
#            for x in pandas.Series([adj_list[q] for q in bio_2020][0])[0].astype(int):
#                not_adj_bio_2020.remove(x)
#            m.addConstrs(((c_8[ii,year]<=quicksum(M[3*L+ii,year]+M[3*L+adj_list[ii-1][0][jj],year-2] for jj in range(len(adj_list[ii-1][0])))) for ii in not_adj_bio_2020 for year in range(4,5)), name="orstatement")

        m.addConstrs(((c_8[ii,1]==1) for ii in start), name = "biohist2")
        m.addConstrs(((c_8[ii,2]==1) for ii in second), name = "biohist2")
        #m.addConstrs(((c_8[ii,1]==0) for ii in nostart), name = "biohist2")

        #Lower density threshold - set d2prime to 0 if dprime<phi (c_4==0)
        m.addConstrs((c_4[loc,year]>=(dprime[loc,year]-phi)/(3000/phi) for loc in sites for year in time), name="LT") #if density is below detection threshold, set to zero
        m.addConstrs(((d2prime[loc,year]>=dprime[loc,year]-((3000/phi)*(1-c_4[loc,year]))) for loc in sites for year in time), name="LT3") #otherwise, keep it at d
        m.addConstrs(((d2prime[loc,year]<=dprime[loc,year]+((3000/phi)*(1-c_4[loc,year]))) for loc in sites for year in time), name="LT3") #otherwise, keep it at d
        m.addConstrs(((d2prime[loc,year]<=((3000/phi)*(c_4[loc,year]))) for loc in sites for year in time), name="LT3") #otherwise, keep it at d

        #Propagule pressure allowed to emigrate

        # propagule pressure impacted by outward quarantine if M[2:L+loc,year]==1 (reduced by 1-eff_quar_out)
        m.addConstrs(((d_out[loc,year] >= d2prime[loc,year]-((3000/phi)*M[2*L+loc,year])) for loc in sites for year in time_sub), name="noquar_2out") #if quarantine out action not chosen, dispersal out is equivalent to total possible dispersal out
        m.addConstrs(((d_out[loc,year] <= d2prime[loc,year]+((3000/phi)*M[2*L+loc,year])) for loc in sites for year in time_sub), name="noquar_2out") #if quarantine out action not chosen, dispersal out is equivalent to total possible dispersal out
        m.addConstrs(((d_out[loc,year] == d2prime[loc,year]) for loc in sites for year in range(6,8)), name="noquar_2out") #if quarantine out action not chosen, dispersal out is equivalent to total possible dispersal out


        #apply quarantine out to calculate new density that can disperse out
        m.addConstrs((d_out[ii,year] >= (1-eff_vec[2*L+ii-1])*d2prime[ii,year]-(3000/phi*(1-M[2*L+ii,year])) for ii in sites for year in time_sub), name = "quar_out") 
        m.addConstrs((d_out[ii,year] >= d2prime[ii,year]-((3000/phi)*(M[2*L+ii,year])) for ii in sites for year in time_sub), name = "quar_out") 
        m.addConstrs((d_out[ii,year] <= d2prime[ii,year] for ii in sites for year in time_sub), name = "quar_out") 
        m.addConstrs((d_out[ii,year] == d2prime[ii,year] for ii in sites for year in range(6,8)), name = "quar_out") 

        #quarantine in calculations require dispersal matrixes
        #dispersal transition matrices from publication (accounts for human population growth), rounding to reduce error when comparing to R output

        Tij = pandas.io.parsers.read_csv("../../data/transmatM__1_6_0_1_0.3_0.3_0.1_.csv")
        Tij= numpy.around(Tij, 6)
        Tij2=pandas.io.parsers.read_csv("../../data/transmatM__1_7_0_1_0.3_0.3_0.1_.csv")
        Tij2= numpy.around(Tij2, 6)
        Tij3 = pandas.io.parsers.read_csv("../../data/transmatM__1_8_0_1_0.3_0.3_0.1_.csv")
        Tij3= numpy.around(Tij3, 6)
        Tij4=pandas.io.parsers.read_csv("../../data/transmatM__1_9_0_1_0.3_0.3_0.1_.csv")
        Tij4= numpy.around(Tij4, 6)
        Tij5 = pandas.io.parsers.read_csv("../../data/transmatM__1_10_0_1_0.3_0.3_0.1_.csv")
        Tij5= numpy.around(Tij5, 6)
        Tij6=pandas.io.parsers.read_csv("../../data/transmatM__1_11_0_1_0.3_0.3_0.1_.csv")
        Tij6= numpy.around(Tij6, 6)
        Tij7=pandas.io.parsers.read_csv("../../data/transmatM__1_12_0_1_0.3_0.3_0.1_.csv")
        Tij7= numpy.around(Tij7, 6)

        #shorten list of possible immigration/emigration sites to reduce problem size
        tij = numpy.stack([Tij,Tij2,Tij3,Tij4,Tij5,Tij6, Tij7])
        full_out = [list([] for x in range(L)),list([] for x in range(L)),list([] for x in range(L)),list([] for x in range(L)),list([] for x in range(L)),list([] for x in range(L)),list([] for x in range(L)),list([] for x in range(L))]
        for i in range(0,7):
            for j in range(L):
                full_out[i][j].append(numpy.array(range(L))[tij[i,range(L),j]!=0]) # which sites are sources in influx to j

        #calculate immigration/emigration depending on quarantine_in management action status
        m.addConstrs(((d3prime[ii,year] >= tij[year-1,ii-1,ii-1]*d2prime[ii,year]+(1-eff_vec[L+ii-1])*quicksum(tij[year-1,j,ii-1]*d_out[j+1,year] for j in numpy.ndarray.tolist(full_out[year-1][ii-1][0]))-tij[year-1,ii-1,ii-1]*d_out[ii,year]-((3000/phi)*(1-M[L+ii,year]))) for ii in sites for year in time_sub), name = "quar_in")
        m.addConstrs(((d3prime[ii,year] >= tij[year-1,ii-1,ii-1]*d2prime[ii,year]+quicksum(tij[year-1,j,ii-1]*d_out[j+1,year] for j in numpy.ndarray.tolist(full_out[year-1][ii-1][0]))-tij[year-1,ii-1,ii-1]*d_out[ii,year]-((3000/phi)*M[L+ii,year]))for ii in sites for year in time_sub), name="noquar_in")
        m.addConstrs(((d3prime[ii,year] <= tij[year-1,ii-1,ii-1]*d2prime[ii,year]+(1-eff_vec[L+ii-1])*quicksum(tij[year-1,j,ii-1]*d_out[j+1,year] for j in numpy.ndarray.tolist(full_out[year-1][ii-1][0]))-tij[year-1,ii-1,ii-1]*d_out[ii,year]+((3000/phi)*(1-M[L+ii,year]))) for ii in sites for year in time_sub), name = "quar_in")
        m.addConstrs(((d3prime[ii,year] <= tij[year-1,ii-1,ii-1]*d2prime[ii,year]+quicksum(tij[year-1,j,ii-1]*d_out[j+1,year] for j in numpy.ndarray.tolist(full_out[year-1][ii-1][0]))-tij[year-1,ii-1,ii-1]*d_out[ii,year]+((3000/phi)*M[L+ii,year]))for ii in sites for year in time_sub), name="noquar_in")
        m.addConstrs(((d3prime[ii,year] == tij[year-1,ii-1,ii-1]*d2prime[ii,year]+quicksum(tij[year-1,j,ii-1]*d_out[j+1,year] for j in numpy.ndarray.tolist(full_out[year-1][ii-1][0]))-tij[year-1,ii-1,ii-1]*d_out[ii,year]) for ii in sites for year in range(6,8)), name="noquar_in")

        #Growth (only if c_5==1)
        m.addConstrs(((d4prime[loc, year] >= d3prime[loc, year]*r-((3000/phi)*(1-c_5[loc,year]))) for loc in sites2 for year in time), name="growth") #allow populations to grow if theyre over phi
        m.addConstrs(((d4prime[loc, year] <= d3prime[loc, year]*r+((3000/phi)*(1-c_5[loc,year]))) for loc in sites2 for year in time), name="growth")
        m.addConstrs(((d4prime[loc, year] >= d3prime[loc, year]) for loc in sites2 for year in time), name="growth")
        #c_5==1 if d3prime>phi
        m.addConstrs(((c_5[loc,year]>=(d3prime[loc,year]-phi)/(3000/phi)) for loc in sites2 for year in time), name="growth")#allow populations to grow if theyre over phi

        ## Cap d4 prime after growth by 1 - convert to starting density at next timestep
        # c_6==1 if d3prime<1, c_7==1 if d3prime>1
        #d at next timestep is equal to d4prime or 1, whichever is lower
        m.addConstrs((((c_6[loc,year]+c_7[loc,year] == 1)) for loc in sites2 for year in time), name="max_d1") #density is either above or below maximum density (1000)
        m.addConstrs(((d[loc, year+1] >= d4prime[loc, year]-((3000/phi)*(1-c_6[loc,year]))) for loc in sites2 for year in time), name="max_d2") #if it's below, don't set to 1000
        m.addConstrs(((d[loc, year+1] <= d4prime[loc, year]) for loc in sites2 for year in time), name="max_d2") #if it's below, don't set to 1000
        m.addConstrs(((d[loc, year+1] >= 1000-((3000/phi)*(1-c_7[loc,year]))) for loc in sites2 for year in time), name="max_d3") #if it is above 1000, reduce back to 1000
        m.addConstrs(((d[loc, year+1] <= 1000+((3000/phi)*(1-c_7[loc,year]))) for loc in sites2 for year in time), name="max_d3") #if it is above 1000, reduce back to 1000
        m.addConstrs(((c_6[loc,year]<=1+((1000-d4prime[loc,year])/(3000/phi))) for loc in sites2 for year in time), name="max_d4")
        m.addConstrs(((c_7[loc,year]>=(d4prime[loc,year]-1000)/(3000/phi)) for loc in sites2 for year in time), name="max_d5")
        m.addConstrs(((c_6[loc,year]>=-c_7[loc,year]+((d4prime[loc,year])/(3000/phi))) for loc in sites2 for year in time), name="max_d6")

        #maintain source at maximum propagule pressure
        m.addConstrs(((d[source,year]==1000) for year in range(1, len(time)+2)), name="source_den") #pest density has to be at maximum in source cell
        m.addConstrs(((d4prime[source,year]==1000*r) for year in range(1, len(time)+1)), name="source_den")
        #no impact of biocontrol in first year beyond what was already included in input data
        m.addConstrs(((dprime[ii,year] == d[ii,year]) for ii in sites for year in range(1,2)), name = "nobiocontrol)" ) #when no biocontrol, don't reduce density

        #keep sites at threshold if management has not taken place
        m.addConstrs(((d4prime[loc, year] >= (c_4[loc,year]*phi*r)-((3000/phi)*(M[L+loc, year]+M[2*L+loc,year]))) for loc in sites for year in range(1,2)) , name="growth")#if no management taken, populations can't go extinct (maintained at phi*r)
        m.addConstrs(((d4prime[loc, year] >= (c_4[loc,year]*phi*r)-((3000/phi)*(M[L+loc, year]+M[2*L+loc,year]+c_8[loc,year-1]))) for loc in sites for year in range(2,3)) , name="growth")#if no management taken, populations can't go extinct (maintained at phi*r)
        m.addConstrs(((d4prime[loc, year] >= (c_4[loc,year]*phi*r)-((3000/phi)*(M[L+loc, year]+M[2*L+loc,year]+c_8[loc,year-1]+c_8[loc,year-2]))) for loc in sites for year in range(3,6)) , name="growth")#if no management taken, populations can't go extinct (maintained at phi*r)
        m.addConstrs(((d4prime[loc, year] >= (c_4[loc,year]*phi*r-c_8[loc,year-1])) for loc in sites for year in range(6,7)) , name="growth")#if no management taken, populations can't go extinct (maintained at phi*r)
        m.addConstrs(((d4prime[loc, year] >= (c_4[loc,year]*phi*r)) for loc in sites for year in range(7,8)) , name="growth")#if no management taken, populations can't go extinct (maintained at phi*r)

        #minimum pest density for biocontrol to be feasible 
        #m.addConstrs(((M[3*L+ii,year]==0) for ii in [value for value in range(1,1800) if numpy.array(vecptime)[value-1,1] < 27.916] for year in range(1,6)), name = "biocontrol5") # biocontrol only above a minimum density (above average density of initially invaded cells)
        m.addConstrs(((M[2*L+loc,year]<=(1+((d[loc,year]-phi)/(3000/phi)))) for loc in sites for year in range(1,6)), name = "biocontrol5") # biocontrol only above a minimum density (above average density of initially invaded cells)
        m.addConstrs(((M[1*L+loc,year]<=(1+((999-d[loc,year])/(3000/phi)))) for loc in sites for year in range(1,6)), name = "biocontrol5") # biocontrol only above a minimum density (above average density of initially invaded cells)

        mgmt=pandas.io.parsers.read_csv("../../output/management_test_1_0_{0}_{1}_.csv".format(eff_quar_in, eff_bio))#starting condition

        #       implement starting condition
        for sites3 in range(1, L*n_actions+1):
            for year in range(1,5+1):
                    M[sites3, year].start=mgmt.iloc[sites3-1,year-1]  

        m.setParam('LogFile', 'eab_parasitoid_{0}_{1}.log'.format(rr,qq)) #unique logfile
        m.setParam('MIPGap',0.01)
        m.setParam('MIPFocus',1)
        m.setParam('RINS',100)
        m.setParam('TimeLimit',28800)
        m.setParam('NodeFileStart', 10000)
        m.setParam('Threads', 4)
        m.setParam('Heuristics',1)
        m.update()

         #%%Solve & Print
        m.optimize()
        M2 = dict(m.getAttr('X', M)) #save management matrix
        M3 = pandas.Series(M2)
        M4 = M3.unstack()
        M4.to_csv('../../output/M_{0}_{1}_nothresh.csv'.format(eff_quar_out,eff_bio), index=False,header=None)
        # m.write('../../output/model_{0}_{1}.mps'.format(eff_quar_out,eff_bio)) #create .mps file for gurobi_cl
       # m.write('../../output/modelstart_{0}_{1}.sol'.format(eff_quar_out,eff_bio)) # use as initial solution in gurobi_cl
        M2 = dict(m.getAttr('X', d)) #save management matrix
        M3 = pandas.Series(M2)
        M4 = M3.unstack()
        M4.to_csv('../../output/vecptime_{0}_{1}_nothresh.csv'.format(eff_quar_out,eff_bio), index=False,header=None)
