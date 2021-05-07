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

#EAB
spp = 1

B = 963.943 # APHIS 2020 budget related to EAB biocontrol + tribal response (no survey)
HostVol = pandas.io.parsers.read_csv('../../data/streettrees_grid.csv')
prez=pandas.io.parsers.read_csv('../../data/prez_clean_gdk.csv')
prez=prez.iloc[:,0]
HostVol = HostVol/1e+06
L = 1799

adj_mat= pandas.io.parsers.read_csv('../../output/adj_list.csv')
adj_mat.astype(int)
adj_list= list([] for x in range(L))
for i in range(L):
    adj_row=adj_mat.iloc[i,:]
    adj_list[i].append(adj_row.where(adj_row!=0).dropna())


source = 613
time = range(1,5+1) #25 years
r = 1.299034706404549 # from published model
phi =0.538410692229749


#%% Preprocessing

#efficiencies

effs_quar=[0.3,0.6,0.9]
effs_bio=[0.1,0.3,0.5]
sites2 = list(range(1,L+1))
sites2.remove(613)
sites = range(1,L+1)
eff_nothing = 0
n_actions = 4
all_comb = range(1, L * n_actions+1)
range_none = range(1, L+1)
range_quar_in = range(L+1, 2*L+1)
range_quar_out= range(2*L+1, 3*L+1)
range_bio = range(3*L+1, 4*L+1)

#Costs
cost_quar_in=list(itertools.repeat(646.863/309, L)) # aphis 2020 biocontrol budget divided by eab quarantined grid cells
cost_quar_out = list(itertools.repeat(646.863/309, L))
cost_bio =list(itertools.repeat(50, L))
cost_nothing = list(itertools.repeat(0, L))
cost_vec = pandas.Series(cost_nothing+ cost_quar_in+ cost_quar_out+ cost_bio)

for rr in range(1,3):
    for qq in range(2,3):

        #%%Create optimization model
        #Budget
        # Add variables
        m = Model('spread')
        d = m.addVars(sites,range(1, len(time)+2), vtype=GRB.CONTINUOUS, name="d",lb=0, ub=5000)  # density in each cell
        M = m.addVars(all_comb, time, vtype=GRB.BINARY, name="M")# decision variables
        B_each = m.addVars(time,vtype=GRB.CONTINUOUS, name="B", lb=0,ub=B)

        dprime = m.addVars(sites,time, vtype=GRB.CONTINUOUS, name="dprime", lb=0,ub=1000)
        d2prime = m.addVars(sites,time, vtype=GRB.CONTINUOUS, name="d2prime", lb=0,ub=1000)
        d3prime = m.addVars(sites,time, vtype=GRB.CONTINUOUS, name="d3prime", lb=0,ub=3000/phi)
        d4prime = m.addVars(sites,time, vtype=GRB.CONTINUOUS, name="d4prime",lb=0,ub=3000/phi)
        d_out=m.addVars(sites, time, vtype=GRB.CONTINUOUS, name="d_out", lb=0,ub=1000)


        c_4 = m.addVars(sites,time, vtype=GRB.BINARY, name="c_4")
        c_5 = m.addVars(sites,time, vtype=GRB.BINARY, name="c_5")
        c_6 = m.addVars(sites,time, vtype=GRB.BINARY, name="c_6")
        c_7 = m.addVars(sites,time, vtype=GRB.BINARY, name="c_7")
        c_8 = m.addVars(sites,range(1,4), vtype=GRB.BINARY, name="c_8")

        #Set objective
        m.setObjective(quicksum(d[loc, year]*(HostVol.iloc[prez.iloc[loc-1]-1,19]+0.000001) for loc in sites for year in range(1+1,len(time)+2)), GRB.MINIMIZE)
         

        #%% Add constraints 

        # Only allow one type of management per cell
        for ii in sites:
            j = (ii, L+ii, (2*L)+ii, (3*L)+ii)
            m.addConstrs(((quicksum(M[loc,year]for loc in j) == 1) for year in time), name = "One Action") 

         # Budget constraint
        m.addConstrs((B_each[year]==quicksum(M[loc,year]*cost_vec[loc-1] for loc in range(1,4*L+1)) for year in time), name = "Budget")
        m.addConstrs((B_each[year] <=B for year in time), name = "Total Budget")


        #%% Not actual constraints but model equations

        #Management
        for loc in sites:
            for year in time:
                #Lower threshold
                m.addConstr(c_4[loc,year]>=(dprime[loc,year]-phi)/(1000/phi), name="LT")
                m.addConstr((c_4[loc,year]==0)>>(d2prime[loc,year]==0), name="LT2")
                m.addConstr((c_4[loc,year]==1)>>(d2prime[loc,year]==dprime[loc,year]), name="LT3")

                #Dispersal 
                #indicators for when quarantining in, out, or both, maxmimum density=no quarantines
                m.addConstr((M[2*L+loc,year]==0)>>(d_out[loc,year] == d2prime[loc,year]), name="noquar_2out")

        for loc in sites2:
            for year in time:
                #Growth
                m.addConstr((c_5[loc,year]==1)>>(d4prime[loc, year] == d3prime[loc, year]*r), name="growth")
                m.addConstr((d4prime[loc, year] >= d3prime[loc, year]), name="growth")
                m.addConstr((M[loc, year]==1)>>(d4prime[loc, year] >= c_4[loc,year]*phi*r) , name="growth")
                m.addConstr((c_5[loc,year]>=(d3prime[loc,year]-phi)/(3000/phi)), name="growth")
                ## Cap d4 prime after growth by 1 - convert to starting density at next timestep
                m.addConstr(((c_6[loc,year]+c_7[loc,year] == 1)), name="max_d1")
                m.addConstr((c_6[loc,year]==1)>>(d[loc, year+1] == d4prime[loc, year]), name="max_d2")
                m.addConstr((c_7[loc,year]==1)>>(d[loc, year+1] == 1000), name="max_d3")
                m.addConstr((c_6[loc,year]<=1+((1000-d4prime[loc,year])/(3000/phi))), name="max_d4")
                m.addConstr((c_7[loc,year]>=(d4prime[loc,year]-1000)/(3000/phi)), name="max_d5")
                m.addConstr((c_6[loc,year]>=-c_7[loc,year]+((d4prime[loc,year])/(3000/phi))), name="max_d6")


        m.addConstrs(((d[source,year]==1000) for year in range(1, len(time)+2)), name="source_den")
        m.addConstrs(((d4prime[source,year]==1000*r) for year in range(1, len(time)+1)), name="source_den")
        m.addConstrs(((dprime[ii,year] == d[ii,year]) for ii in sites for year in range(1,4)), name = "nobiocontrol)" )


        m.setParam('MIPGap', 0.1)
        m.setParam('Presolve', 2)
        m.setParam('LogFile', 'eab_parasitoid_{0}_{1}.log'.format(rr,qq))
        m.setParam('Threads', 2)
        #m.setParam('Heuristics', 1)
        m.setParam('SolutionLimit',1)
        m.setParam('NodeFileStart',1000)
        m.update()

        eff_quar_in = effs_quar[rr]
        eff_quar_out=effs_quar[rr]
        eff_bio = effs_bio[qq]
        eff_vec = pandas.Series(list(itertools.repeat(eff_nothing, L))+ list(itertools.repeat(eff_quar_in, L))+list(itertools.repeat(eff_quar_out, L))+list(itertools.repeat(eff_bio, L)))
        for ii in sites:
            for year in time:
                m.addConstr((M[2*L+ii,year]==1)>>(d_out[ii,year] == (1-eff_vec[2*L+ii-1])*d2prime[ii,year]), name = "quar_out") 
        for ii in sites:
            for year in range(3,5):
                m.addConstr((c_8[ii,year-2]==1)>>(dprime[ii,year] == (1-eff_vec[3*L+ii-1])*d[ii,year]), name = "biocontrol")
                m.addConstr((M[3*L+ii,year-2]==1)>>(c_8[ii,year-2]==1), name = "biocontrol2")
                m.addConstr(c_8[ii,year-2]<=M[3*L+ii,year-2], name="orstatement")
                m.addConstr((c_8[ii,year-2]==0)>>(dprime[ii,year] == d[ii,year]), name = "nobiocontrol)" )

        for ii in sites:
            for year in range(5,6):
                m.addConstr((c_8[ii,year-2]<=quicksum(M[3*L+ii,year-2]+M[3*L+adj_list[ii-1][0][jj],year-4] for jj in range(len(adj_list[ii-1][0])))), name="orstatement")
                m.addConstr((c_8[ii,year-2]==0)>>(dprime[ii,year] == d[ii,year]), name = "nobiocontrol)" )
                m.addConstr((c_8[ii,year-2]==1)>>(dprime[ii,year] == (1-eff_vec[3*L+ii-1])*d[ii,year]), name = "biocontrol")
                m.addConstr((M[3*L+ii,year-2]==1)>>(c_8[ii,year-2]==1), name = "biocontrol2")
                for jj in adj_list[ii-1][0].astype(int):
                    m.addConstr((M[3*L+jj,year-4]==1)>>(c_8[ii,year-2]==1), name = "parasitoidDisp") 
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
        vecptime = pandas.io.parsers.read_csv("../../output/vecptime_0_1_0.3_0.3_.csv")*1000
        mgmt=pandas.io.parsers.read_csv("../../output/management_test_0_1_0.3_0.3_.csv")


        tij = numpy.stack([Tij,Tij2,Tij3,Tij4,Tij5,Tij6])
        full_out = [list([] for x in range(L)),list([] for x in range(L)),list([] for x in range(L)),list([] for x in range(L)),list([] for x in range(L)),list([] for x in range(L))]
        for i in range(0,6):
            for j in range(L):
                full_out[i][j].append(numpy.array(range(1,L+1))[tij[i,range(L),j]!=0]) # which sites are sources in influx to j

        for ii in sites:
            for year in time:
                m.addConstr((M[L+ii,year]==1)>>(d3prime[ii,year] == tij[year-1,ii-1,ii-1]*d2prime[ii,year]+(1-eff_vec[L+ii-1])*(quicksum(tij[year-1,j-1,ii-1]*d_out[j,year] for j in numpy.ndarray.tolist(full_out[year-1][ii-1][0]))-tij[year-1,ii-1,ii-1]*d_out[ii,year])), name = "quar_in")
                m.addConstr((M[L+ii,year]==0)>>(d3prime[ii,year] == tij[year-1,ii-1,ii-1]*d2prime[ii,year]+quicksum(tij[year-1,j-1,ii-1]*d_out[j,year] for j in numpy.ndarray.tolist(full_out[year-1][ii-1][0]))-tij[year-1,ii-1,ii-1]*d_out[ii,year]), name="noquar_in")
        m.addConstrs(((d[loc,1]==vecptime.iloc[loc-1,0]) for loc in sites2), name="initial_den")


        for sites3 in range(1, L*n_actions+1):
            for year in range(1,5+1):
                M[sites3, year].start=mgmt.iloc[sites3-1,year-1]    
        # for year in range(1,1+1):
        #     B_each[year].start=sum(mgmt.iloc[range(0,4*1799),year-1]*cost_vec)        
        #for sites3 in range(1, L+1):
        #    for year in range(1,5+1):
        #        d_out[sites3,year].start=douts.iloc[sites3-1,year-1]
        #        d2prime[sites3,year].start = d2primes.iloc[sites3-1,year-1]
        #        d3prime[sites3,year].start = d3primes.iloc[sites3-1,year-1]
        #        d4prime[sites3,year].start = d4primes.iloc[sites3-1,year-1]
        #        c_4[sites3,year].start = c_4s.iloc[sites3-1,year-1]
        #        c_5[sites3,year].start = c_5s.iloc[sites3-1,year-1]
        #        c_6[sites3,year].start = c_6s.iloc[sites3-1,year-1]
        #        c_7[sites3,year].start = c_7s.iloc[sites3-1,year-1]
        #        d[sites3,year].start = vecptime.iloc[sites3-1,year-1]  
        #    for year in range(6,6+1):
        #        d[sites3,year].start = vecptime.iloc[sites3-1,year-1]
        m.update()

         #%%Solve & Print
        #m.optimize()
        # M2 = dict(m.getAttr('X', M))
        # M3 = pandas.Series(M2)
        # M4 = M3.unstack()
        # M4.to_csv('../../output/M3_{0}_{1}_{2}.csv'.format(eff_quar_in,eff_quar_out,eff_bio), index=False,header=False)
        m.write('../../output/model_{0}_{1}_{2}.mps'.format(eff_quar_in,eff_quar_out,eff_bio))