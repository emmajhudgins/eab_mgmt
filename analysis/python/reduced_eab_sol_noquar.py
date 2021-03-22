supercomputer = False
import sys
import commands
from gurobipy import *
import numpy
import pandas
import itertools
import random
import scipy
import re

if supercomputer == True:
    args = sys.argv
    nm = args[1]
    thread= args[2]
    control_path=commands.getoutput("cat sc_full_path")
    i=int(commands.getoutput("echo `ssh snowi '~/get_rep 1 %s %s %s'`"%(nm,thread,control_path)))

    
if supercomputer == False:
    nm='offline'
    thread=0
    args=sys.argv
    i=1



#%% Spread model data

#EAB
spp = 1

B = 963.943 # APHIS 2020 budget related to EAB biocontrol + tribal response (no survey)
HostVol = pandas.io.parsers.read_csv('streettrees_grid.csv')
prez=pandas.io.parsers.read_csv('prez_clean_gdk.csv')
prez=prez.iloc[:,0]
HostVol = HostVol/1e+06


source = 613


L = 1799

time = range(1,5+1) #25 years
r = 1.299034706404549 # from published model
phi =0.538410692229749


#%% Preprocessing

#efficiencies

eff_erad = 0.99 #maybe test sensitivity but don't present
erad_thresh=1 #0.1% of cell covered is max can eradicate (2.5km2)
#maybe differentiate rapid response from larger scale eradication
effs_quar=[0.3,0.6,0.9]
effs_bio=[0.1,0.3,0.5]
sites2 = range(1,L+1)
sites2.remove(613)
sites = range(1,L+1)
eff_nothing = 0
n_actions = 5
all_comb = range(1, L * n_actions+1)
range_none = range(1, L+1)
range_quar_in = range(L+1, 2*L+1)
range_quar_out= range(2*L+1, 3*L+1)
range_bio = range(3*L+1, 4*L+1)
range_erad = range(4*L+1, 5*L+1)

#Costs
# alternative cost (2269665/2500)
cost_erad = (31*2.500)#sharov & liebhold,
cost_quar_in=list(itertools.repeat(646.863/309, L)) # aphis 2020 biocontrol budget divided by eab quarantined grid cells
cost_quar_out = list(itertools.repeat(646.863/309, L))
cost_bio =list(itertools.repeat(141.519, L))
cost_nothing = list(itertools.repeat(0, L))
cost_vec = pandas.Series(cost_nothing+ cost_quar_in+ cost_quar_out+ cost_bio)
#cost_vec = cost_vec/1000


#%%Create optimization model
#Budget
# Add variables
m = Model('spread')
d = m.addVars(sites,range(1, len(time)+2), vtype=GRB.CONTINUOUS, name="d")  # density in each cell
M = m.addVars(all_comb, time, vtype=GRB.BINARY, name="M")# decision variables

dprime = m.addVars(sites,time, vtype=GRB.CONTINUOUS, name="dprime") 
d2prime = m.addVars(sites,time, vtype=GRB.CONTINUOUS, name="d2prime")
d3prime = m.addVars(sites,time, vtype=GRB.CONTINUOUS, name="d3prime")
d4prime = m.addVars(sites,time, vtype=GRB.CONTINUOUS, name="d4prime")
d_out=m.addVars(sites, time, vtype=GRB.CONTINUOUS, name="d_out")
B_nonerad = m.addVars(time,vtype=GRB.CONTINUOUS, name="B_nonerad")  
B_erad = m.addVars(sites,time,vtype=GRB.CONTINUOUS, name="B_erad")  


c_4 = m.addVars(sites,time, vtype=GRB.BINARY, name="c_4")
c_5 = m.addVars(sites,time, vtype=GRB.BINARY, name="c_5")
c_6 = m.addVars(sites,time, vtype=GRB.BINARY, name="c_6")
c_7 = m.addVars(sites,time, vtype=GRB.BINARY, name="c_7")

   

#%% Add constraints 

# Only allow one type of management per cell
for ii in sites:
    j = (ii, L+ii, (2*L)+ii, (3*L)+ii, (4*L)+ii)
    m.addConstrs(((quicksum(M[loc,year]for loc in j) == 1) for year in time), name = "One Action") 

m.addConstrs((((d[loc,year]-phi)/1000+(1-M[4*L+loc,year])>=0) for loc in sites2 for year in time), name = "no erad when abs") 

m.addConstrs((((d[loc,year]-erad_thresh)/1000+M[4*L+loc,year]<=1) for loc in sites2 for year in time), name = "no erad above thresh") 

   # Budget constraint
m.addConstrs((B_nonerad[year]==quicksum(M[loc,year]*cost_vec[loc-1] for loc in range(1,4*L+1)) for year in time), name = "Budget")
m.addConstrs((B_erad[loc,year]>=0 for loc in sites for year in time), name = "Budget")
m.addConstrs((B_erad[loc,year]>=((d[loc,year])*cost_erad)-cost_erad*1000*(1-M[4*L+loc,year]) for loc in sites for year in time), name = "Budget")
m.addConstrs((B_erad[loc,year]<=cost_erad*1000*(M[4*L+loc,year]) for loc in sites for year in time), name = "Budget")
m.addConstrs((B_nonerad[year]+quicksum(B_erad[loc,year]for loc in sites) <=B for year in time), name = "Budget")


#%% Not actual constraints but model equations

#m.addConstrs(((d[loc,1]==d0.iloc[loc-1]) for loc in sites), name="initial d")

#Management
for j in (2,3):
    m.addConstrs(((dprime[loc,year] >= d[loc,year]-1000*(1-M[(j-1)*L+loc,year])) for loc in sites for year in time), name = "Effect of 2Mgmt" )
    m.addConstrs(((dprime[loc,year] <= d[loc,year]) for loc in sites for year in time), name = "Effect of 2Mgmt" )



#Lower threshold
m.addConstrs(((d2prime[loc,year]>=0) for loc in sites for year in time), name="LT")
m.addConstrs(((c_4[loc,year]>=(dprime[loc,year]-phi)/1000) for loc in sites for year in time), name="LT")
m.addConstrs(((c_4[loc,year]<=1+((dprime[loc,year]-phi)/1000)) for loc in sites for year in time), name="LT")

m.addConstrs(((d2prime[loc,year]<=1000*(c_4[loc,year])) for loc in sites for year in time), name="LT")
m.addConstrs(((d2prime[loc,year]>=dprime[loc,year]-1000*(1-c_4[loc,year])) for loc in sites for year in time), name="LT")

#Dispersal 
#indicators for when quarantining in, out, or both, maxmimum density=no quarantines
m.addConstrs((d_out[ii,year] >= d2prime[ii,year]-1000*(M[2*L+ii,year]) for ii in sites for year in time), name="noquar_2out")
m.addConstrs((d_out[ii,year] >= 0 for ii in sites for year in time), name="quar_2out")
m.addConstrs((d_out[ii,year] <= d2prime[ii,year] for ii in sites for year in time), name = "quar_2out" )


m.addConstrs((d3prime[ii,year] >= 0 for ii in sites for year in time), name="quar_2in")

#Growth

m.addConstrs(((d4prime[loc, year] >= d3prime[loc, year]) for loc in sites2 for year in time), name="growth")
m.addConstrs(((d4prime[loc, year] >=(d3prime[loc, year]*r)-(1000*(1-c_5[loc,year]))) for loc in sites2 for year in time), name="growth")
m.addConstrs(((d4prime[loc, year] <=(phi*r)+3*(1000/phi)*(c_5[loc,year])) for loc in sites2 for year in time), name="growth")
m.addConstrs((((c_5[loc,year] >=(d3prime[loc,year]-phi)/(1000/phi))) for loc in sites for year in time), name="growth")
m.addConstrs((((c_5[loc,year] <=(1+((d3prime[loc,year]-phi)/(1000/phi))))) for loc in sites for year in time), name="growth")

m.addConstrs(((d4prime[loc, year] >= (M[loc, year]*phi*r-3*(1000/phi)*(1-c_4[loc,year]))) for loc in sites2 for year in range(1,6)),name="noext")
#m.addConstrs(((d4prime[loc, 1] >= (M[loc, 1]*phi*r-5*(1-c_4[loc,1]))) for loc in range(1, L+1)),name="noext")
#m.addConstrs(((c_5[loc,1] >=d[loc,1]-phi) for loc in sites), name="noext")


## Cap d4 prime after growth by 1 - convert to starting density at next timestep

m.addConstrs(((c_6[loc,year]+c_7[loc,year] == 1) for loc in sites for year in time), name="max_d1")
m.addConstrs(((d[loc, year+1] >= d4prime[loc, year]-100*(1000*(1-c_6[loc,year]))) for loc in sites for year in time), name="max_d2")
m.addConstrs(((d[loc, year+1] >= 1000*c_7[loc,year]) for loc in sites for year in time), name="max_d3")
m.addConstrs(((c_6[loc,year] >= (1000-d4prime[loc,year])/(1000/phi)) for loc in sites for year in time), name="max_d4")
m.addConstrs((c_7[loc,year] >= ((d4prime[loc,year]-1000)/(1000/phi)) for loc in sites for year in time), name="max_d5")
   

m.addConstrs(((d[source,year]==1000) for year in range(1, len(time)+1)), name="source_den")
m.addConstrs(((d4prime[source,year]==1000*r) for year in range(1, len(time)+1)), name="source_den")

m.addConstrs((d[loc, year]<=(1000) for loc in sites for year in range(1, len(time)+2)),name="noext")

m.setObjective(quicksum(d[loc, year]*(HostVol.iloc[prez.iloc[loc-1]-1,19]+0.000001) for loc in sites for year in range(1+1,len(time)+2)), GRB.MINIMIZE)
# m.setParam('MIPGap', 0.1)
# m.setParam('Presolve', 2)
# m.setParam('MIPFocus', 3)
# m.setParam('Heuristics', 1)
# m.setParam('Presparsify', 1)
# m.setParam('Cuts', 3)
# m.setParam('Aggregate',1)
# m.setParam('SubMIPNodes',1000)
# m.setParam('StartNodeLimit', 1e+09)
m.setParam('OptimalityTol', 1e-05)
m.setParam('FeasibilityTol', 1e-05)
m.setParam('IntFeasTol', 1e-05)
# m.setParam('ImproveStartTime', 3600)
#    m.setParam('MinRelNodes', -1)
#    m.setParam('ZeroObjNodes', -1)
#    m.setParam('PumpPasses', -1)
m.setParam('TimeLimit', 3600)
  
m.update()
            
           
eff_quar_in = 0.9
eff_quar_out=0.9
eff_bio = 0.3
eff_vec = pandas.Series(list(itertools.repeat(eff_nothing, L))+ list(itertools.repeat(eff_quar_in, L))+list(itertools.repeat(eff_quar_out, L))+list(itertools.repeat(eff_bio, L))+list(itertools.repeat(eff_erad, L)))
m.addConstrs((d_out[ii,year] >= ((1-eff_vec[2*L+ii-1])*d2prime[ii,year]-1000*(1-M[2*L+ii,year])) for ii in sites for year in time), name = "quar_out" )
for j in (1,4,5):
    m.addConstrs(((dprime[loc,year] >= ((1-eff_vec[(j-1)*L+loc-1])*d[loc,year]-1000*(1-M[(j-1)*L+loc,year]))) for loc in sites for year in time), name = "Effect of Mgmt" )

Tij = pandas.io.parsers.read_csv("transmatM__1_6_0_1_0.3_0.3_0.1_.csv")
Tij= numpy.around(Tij, 6)
Tij2=pandas.io.parsers.read_csv("transmatM__1_7_0_1_0.3_0.3_0.1_.csv")
Tij2= numpy.around(Tij2, 6)
Tij3 = pandas.io.parsers.read_csv("transmatM__1_8_0_1_0.3_0.3_0.1_.csv")
Tij3= numpy.around(Tij3, 6)
Tij4=pandas.io.parsers.read_csv("transmatM__1_9_0_1_0.3_0.3_0.1_.csv")
Tij4= numpy.around(Tij4, 6)
Tij5 = pandas.io.parsers.read_csv("transmatM__1_10_0_1_0.3_0.3_0.1_.csv")
Tij5= numpy.around(Tij5, 6)
Tij6=pandas.io.parsers.read_csv("transmatM__1_11_0_1_0.3_0.3_0.1_.csv")
Tij6= numpy.around(Tij6, 6)
c_4s = pandas.io.parsers.read_csv("c_4_0_1_0.3_0.3_0.1_.csv")
c_5s = pandas.io.parsers.read_csv("c_5_0_1_0.3_0.3_0.1_.csv")
c_6s = pandas.io.parsers.read_csv("c_6)0)1)0.3)0.3)0.1).csv")
c_7s = pandas.io.parsers.read_csv("c_7)0)1)0.3)0.3)0.1).csv")
d2primes = pandas.io.parsers.read_csv("d2prime_0_1_0.3_0.3_0.1_.csv")*1000
d3primes = pandas.io.parsers.read_csv("d3prime_0_1_0.3_0.3_0.1_.csv")*1000
d4primes = pandas.io.parsers.read_csv("d4prime_0_1_0.3_0.3_0.1_.csv")*1000
vecptime = pandas.io.parsers.read_csv("vecptime_0_1_0.3_0.3_0.1_.csv")*1000
dprimes = pandas.io.parsers.read_csv("dprime_0_1_0.3_0.3_0.1_.csv")*1000
douts = pandas.io.parsers.read_csv("d_out_0_1_0.3_0.3_0.1_.csv")*1000
mgmt=pandas.io.parsers.read_csv("management_test_0_1_0.3_0.3_0.1_.csv")


tij = numpy.stack([Tij,Tij2,Tij3,Tij4,Tij5,Tij6])
full_out = [list([] for x in xrange(L)),list([] for x in xrange(L)),list([] for x in xrange(L)),list([] for x in xrange(L)),list([] for x in xrange(L)),list([] for x in xrange(L))]
for i in range(0,6):
    for j in range(L):
        full_out[i][j].append(numpy.array(range(1,L+1))[tij[i,range(L),j]!=0]) # which sites are sources in influx to j
m.addConstrs((d3prime[ii,year] >= tij[year-1,ii-1,ii-1]*d2prime[ii,year]+(1-eff_vec[L+ii-1])*(quicksum(tij[year-1,j-1,ii-1]*d_out[j,year] for j in numpy.ndarray.tolist(full_out[year-1][ii-1][0]))-tij[year-1,ii-1,ii-1]*d_out[ii,year])-(1000/phi)*(1-M[L+ii,year]) for ii in sites for year in time), name = "quar_in")
m.addConstrs(((d3prime[ii,year] >= tij[year-1,ii-1,ii-1]*d2prime[ii,year]+quicksum(tij[year-1,j-1,ii-1]*d_out[j,year] for j in numpy.ndarray.tolist(full_out[year-1][ii-1][0]))-tij[year-1,ii-1,ii-1]*d_out[ii,year]-(1000/phi)*(1-M[ii,year]-M[4*L+ii,year]-M[3*L+ii,year]-M[2*L+ii,year])) for ii in sites for year in time), name="noquar_in")
m.addConstrs(((d3prime[ii,year] <= tij[year-1,ii-1,ii-1]*d2prime[ii,year]+quicksum(tij[year-1,j-1,ii-1]*d_out[j,year] for j in numpy.ndarray.tolist(full_out[year-1][ii-1][0]))-tij[year-1,ii-1,ii-1]*d_out[ii,year]) for ii in sites for year in time), name="noquar_in")
m.addConstrs(((d[loc,1]==vecptime.iloc[loc-1,0]) for loc in sites2), name="initial_den")
       # m.write('model_{0}_{1}.mps'.format(qin,qbio))
# for sites3 in range(1, L*n_actions+1):
#     for year in range(1,5+1):
#         M[sites3, year].start=mgmt.iloc[sites3-1,year-1]    
# for year in range(1,5+1):
#     B_nonerad[year].start=sum(mgmt.iloc[range(0,4*1799),year-1]*cost_vec)        
# #for sites3 in range(1, L+1):
# #    for year in range(1,5+1):
# #        B_erad[sites3, year].start=mgmt.iloc[4*1799+sites3-1,year-1]*vecptime.iloc[sites3-1,year-1]*cost_erad
# #        d_out[sites3,year].start=douts.iloc[sites3-1,year-1]
# #        d2prime[sites3,year].start = d2primes.iloc[sites3-1,year-1]
# #        d3prime[sites3,year].start = d3primes.iloc[sites3-1,year-1]
# #        d4prime[sites3,year].start = d4primes.iloc[sites3-1,year-1]
# #        c_4[sites3,year].start = c_4s.iloc[sites3-1,year-1]
# #        c_5[sites3,year].start = c_5s.iloc[sites3-1,year-1]
# #        c_6[sites3,year].start = c_6s.iloc[sites3-1,year-1]
# #        c_7[sites3,year].start = c_7s.iloc[sites3-1,year-1]
# #        d[sites3,year].start = vecptime.iloc[sites3-1,year-1]  
# #        dprime[sites3,year].start=dprimes.iloc[sites3-1,year-1]
# #    for year in range(6,6+1):
# #        d[sites3,year].start = vecptime.iloc[sites3-1,year-1]
# m.update()
# m.setParam('NodeFileStart',500)
# #m.setParam('Method',3)
# #m.setParam('VarBranch',3)
# #m.setParam('Cuts',0)
                        
#             #%%Solve & Print
# m.optimize()
m.read('model_{0}_noquar.sol'.format(eff_bio))
m.setParam('SolutionLimit',1)
m.optimize()
M2 = dict(m.getAttr('X', M))
M3 = pandas.Series(M2)
M4 = M3.unstack()
M4.to_csv('M3_{0}_noquar.csv'.format(eff_bio), index=False,header=False)
# m.write('model_{0}_{1}_{2}.mps'.format(eff_quar_in,eff_quar_out,eff_bio))
#            c = m.getConstrs()
#            IIS=m.getAttr('IISConstrs')
#            whichin=pandas.Series([i==1 for i in IIS])
#            c2=numpy.where(whichin==True)[0]
#            [c[i] for i in numpy.array(c2)]
#input()
