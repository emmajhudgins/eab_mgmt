<img src="https://github.com/emmajhudgins/eab_mgmt/blob/main/results_preview.png" alt="results preview" width="400" align="right"/>

_repo to accompany_
# Optimal Emerald Ash Borer (*Agrilus planipennis*) control across the United States
Emma J. Hudgins, Jeffrey O. Hanson, Christian MacQuarrie, Denys Yemshanov, Eve McDonald-Madden, Chris Baker, Matthew Holden, Iadine Chadès, & Joseph R. Bennett

Substantial conceptual contributions by Martin Péron, past PhD student of Iadine Chadès.

An optimal control approach to EAB management for urban tree persistence in explicit space and time.

**The larger postdoc project context:**

The goal is to find rules of thumb for the optimal management of North American invasive forest pests, and to create an interactive open source tool to display these management optimizations.

**Previous work:**

I previously created a cost forecast for the damages caused by 57 US invasive insect pest species to street trees in the next 30 years. To do this, a model I built for countrywide urban trees, spread model predictions from another part of my thesis, and a model of host tree death given pest exposure. I found that damages will be dominated by emerald ash borer induced ash tree mortality, and that the Midwest would see the greatest costs.

**Motivation for EAB pilot study**	

Given that the overwhelming majority of forecasted damages were due to EAB, combined with recent news that the EAB quarantine was being lifted, I thought it would make a good initial test for this approach using US data where I already had models built. The eventual goal is to build similar models for Canada and apply the framework to several pests of interest to ascertain if there are similarities in best strategies across species.


**Methods:**

we reformulated the pest dispersal model as an optimal control framework. Quarantines could limit propagule pressure entering or leaving cells, and biocontrol and eradication could decrease population sizes in each cell.
Contrary to past frameworks, this approach accounts for the interactivity of management implications across timesteps. 
We group potential strategies into site and spread focussed, either tackling flow of propagules or local population sizes. 
We tested our optimization model against rules of thumb that used different proportions of site and spread-focussed measures
We chose to minimize future exposed street ash (gain in propagule pressure*# urban trees)

**Results:** 

We found that the best management strategy always included a combination of site-focussed (eradication and biological control) and spread-focussed (quarantine) management measures. Optimal strategies vastly outperformed common rules-of-thumb such as targeting the sites causing the greatest number of new establishments or receiving the highest propagule pressure for control actions. 

These findings support a multipronged EAB management approach, where regional context is considered, and quarantine regions are refined rather than abandoned. Quarantines targeting inflow of propagules are high priority for preventing new jumps to areas like Seattle, while quarantines focussing on outflow dominate the core of the invasion. Quarantines  get replaced with biocontrol in more efficient scenarios, and Eradications are rare and are only recommended when biocontrol is inefficient.

**File organization structure**

The analysis section is split into R and Python (3) components. The python scripts require the gurobipy library, which interfaces with GUROBI 9.1.0 software (free for academic use, www.gurobi.com). All R components were build in R 4.1.0.

R Scripts
1. 1-biocontrol_history.R - this script uses parasitoid release and recapture data from APHIS' MapBioControl database (see www.mapbiocontrol.org for data inquiries) to examine recapture patterns across space and time for the four parasitoid species for EAB, which are then used to indicate which sites had a history of biocontrol release at a level that could impact EAB density in line with biocontrol efficiency scenarios (10,000 released and 100 recovered parasitoids).
2. SDP_calcRoT_spreadfirst_noerad.R - calculates rule of thumb scenarios across different budget allocations to quarantines vs. biologican control (in 20% increments) and different levels of efficiency of each action. Saves management action locations over time and resulting EAB propagule pressure for plotting and comparison with optimizations. Spread-focused management (quarantines) is allocated first in this script.
3. SDP_calcRoT_spreadsecon_noerad.R - calculates rule of thumb scenarios across different budget allocations to quarantines vs. biologican control (in 20% increments) and different levels of efficiency of each action. Saves management action locations over time and resulting EAB propagule pressure for plotting and comparison with optimizations. Site-focused management (biological control) is allocated first in this script.
4. evaluate_SDP.R - this script checks the value of the objective function for both optimized runs and rule-of-thumb scenarios and examines the proportion of different actions taken 
5. evaluate_SDP_budget.R - this script checks the value of the objective function for both optimized runs, budget-allocation, and sensitivity analysis and plots the comparison between objective values.
6. pp_compare.R - this script performs the ANOVA and Tukey tests between the propagule pressure in cells allocated to different action types across optimization scenarios.
7. SDP_plot.R - this script plots the spatial pattern in propagule pressure, exposed ash street trees, and management actions for both optimizations and rule of thumb scenarios, as well as more comparisons between optimality and rules of thumb.


Python Scripts
1. eab_parasitoid.py - basic optimization across efficiency scenarios
2. eab_parasitoid_bud.py - fixed percentage budget allocation optimizations (20% increments)
3. eab_parasitoid_fastdisp_fastgrowth.py - sensitivity analysis to faster parasitoid dispersal and growth
4. eab_parasitoid_fastdisp.py - sensitivity analysis to faster parasitoid dispersal
5. eab_parasitoid_nothresh.py - sensitivity analysis with removal of requirement to impose biological control in high-EAB density sites.
