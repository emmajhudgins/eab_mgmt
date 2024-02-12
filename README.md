<img src="https://github.com/emmajhudgins/eab_mgmt/blob/main/results_preview.png" alt="results preview" width="400" align="right"/>

_repo to accompany_
# Spread management priorities to limit emerald ash borer (*Agrilus planipennis*) impacts to United States street trees

An optimal control approach to EAB management for urban tree persistence in explicit space and time.


**Methods:**

We reformulated an existing pest dispersal model as an optimal control framework. Quarantines could limit propagule pressure entering or leaving cells, and biocontrol and eradication could decrease population sizes in each cell.

Contrary to past frameworks, this approach accounts for the interactivity of management implications across timesteps. 

We group potential strategies into site and spread focussed, either tackling flow of propagules or local population sizes. 

We tested our optimization model against rules of thumb that used different proportions of site and spread-focussed measures.

We chose to minimize future exposed street ash (EAB propagule pressure* number of urban trees)

**Results:** 

We found that the best management strategy always included a combination of site-focussed (eradication and biological control) and spread-focussed (quarantine) management measures. Optimal strategies vastly outperformed common rules-of-thumb such as targeting the sites causing the greatest number of new establishments or receiving the highest propagule pressure for control actions. 

These findings support a multipronged EAB management approach, where regional context is considered, and quarantine regions are refined rather than abandoned. Quarantines targeting inflow of propagules are high priority for preventing new jumps to areas like Seattle, while quarantines focussing on outflow dominate the core of the invasion. 

**File organization structure**

The analysis section is split into R and Python (3) components. The python scripts require the gurobipy library, which interfaces with GUROBI software (free for academic use, www.gurobi.com; built in v9.1.0). All R components were built in R 4.1.0.

### R Scripts
0. 0-forecasted_presences_eab.R - this script updates a model of EAB spread to allow establishment anywhere there is urban ash, and to fit urban ash density-specific spread terms(https://doi.org/10.1111/1365-2664.14141), producing forecasts of EAB spread to 2050 in the absence of management and a dispersal kernel for use in optimizations.
1. 1-biocontrol_history.R - this script uses private parasitoid release and recapture data from APHIS' MapBioControl database accessed with an MOU agreement (see www.mapbiocontrol.org for data inquiries) to examine recapture patterns across space and time for the four parasitoid species for EAB, which are then used to indicate which sites had a history of biocontrol release at a level that could impact EAB density in line with biocontrol efficiency scenarios (10,000 released and 100 recovered parasitoids).
2. SDP_calcRoT_spreadfirst_noerad_street.R - calculates rule of thumb scenarios across different budget allocations to quarantines vs. biologican control (in 20% increments) and different levels of efficiency of each action. Saves management action locations over time and resulting EAB propagule pressure for plotting and comparison with optimizations. Spread-focused management (quarantines) is allocated first in this script.
3. SDP_calcRoT_spreadsecon_noerad_street.R - calculates rule of thumb scenarios across different budget allocations to quarantines vs. biologican control (in 20% increments) and different levels of efficiency of each action. Saves management action locations over time and resulting EAB propagule pressure for plotting and comparison with optimizations. Site-focused management (biological control) is allocated first in this script.
4. evaluate_SDP.R - this script checks the value of the objective function for both optimized runs and rule-of-thumb scenarios and examines the proportion of different actions taken 
5. evaluate_SDP_budget.R - this script checks the value of the objective function for both optimized runs, budget-allocation, and sensitivity analysis and plots the comparison between objective values.
6. pp_compare.R - this script performs the ANOVA and Tukey tests between the propagule pressure in cells allocated to different action types across optimization scenarios.
7. SDP_plot.R - this script plots the spatial pattern in propagule pressure, exposed ash street trees, and management actions for both optimizations and rule of thumb scenarios, as well as more comparisons between optimality and rules of thumb.
8. cost_estimate.R - this script calculates a rough estimate of the maximal cost savings from this approach relative to a biocontrol-only, rule-of-thumb approach using previously published ash size class distributions and cost data from Hudgins et al. (2022) and Aukema et al. (2011).
9. SDP_calcRoT_spreadfirst_noerad_street_noaction.R - this script calculates ash exposure in the absence of any management action by setting the management budget to $0.

### Python Scripts
1. eab_parasitoid_street.py - basic optimization across efficiency scenarios
2. eab_parasitoid_street_bud.py - fixed percentage budget allocation optimizations (20% increments)
3. eab_parasitoid_street_fastdisp_fastgrowth.py - sensitivity analysis to faster parasitoid dispersal and growth (produces \*_fdfg\* files)
4. eab_parasitoid_street_fastdisp.py - sensitivity analysis to faster parasitoid dispersal
5. eab_parasitoid_street_nothresh.py - sensitivity analysis with removal of requirement to impose biological control in high-EAB density sites.
6. eab_parasitoid_street_extend.py - sensitivity analysis to extending the solving time by another 8 hours, loading the solve from script 1 as an initial solution.
7. eab_parasitoid_street_upper.py - sensitivity analysis where the upper EAB density in sites selected for management was capped at 0.75
8. eab_parasitoid_street_multact.py - sensitivity analysis where >1 management action types were allowed to be selected in a site
9. eab_parasitoid_street_doublecost.py - sensitivity analysis where the cost of biocontrol was doubled
10. eab_parasitoid_street_halfcost.py - sensitivity analysis where the cost of biocontrol was halved

