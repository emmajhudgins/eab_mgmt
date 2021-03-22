_repo to accompany_
# Optimal Emerald Ash Borer (Agrilus planipennis) control across the United States
Emma J. Hudgins, Jeffrey O. Hanson, Christian MacQuarrie, Denys Yemshanov, Richard Schuster, Eve McDonald-Madden, Chris Baker, Matthew Holden, Iadine Chadès, & Joseph R. Bennett

An optimal control approach to EAB management for urban tree persistence in explicit space and time.

**The larger postdoc project context:**

The goal is to find rules of thumb for the optimal management of Canadian invasive forest pests, and to create an interactive open source tool to display these management optimizations.

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

**Current Work** 
- Better estimates of biocontrol cost incorporated
- Incorporating a model for the dispersal of biocontrol agents themselves. Chris MacQuarrie has data that could potentially be used to estimate a very basic parasitoid dispersal kernel.
- Eradication is not a viable EAB strategy due to the high detection threshold of EAB, so I will be reformulating eradication as an insecticide treatment protocol
- incorporating a time lag such that biocontrol can only be feasible after a certain number of years of EAB infestation (potentially estimating from US county detection and MapBioControl)
- create some budget allocation scenarios rather than allowing any % of the budget to go to biocontrol, because this would likely have more uptake by government, and could show diminishing returns of increased biocontrol investment (but keep in mind different jusidictions for each management option)
