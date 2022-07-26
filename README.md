_repo to accompany_
# Optimal Emerald Ash Borer (Agrilus planipennis) control across the United States
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

** File organization structure **

The analysis section is split into R and Python (3) components. The python scripts require the gurobipy library, which interfaces with GUROBI 9.1.0 software (free for academic use, [!www.gurobi.com])