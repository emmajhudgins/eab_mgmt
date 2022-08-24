#R script to perform a rough estimate of the cost savings associated with the largest improvement between optimal and rule-of-thumb scenarios.

#trees
pred_trees_huge<-readRDS('./data/pred_trees_huge_public.RDS') #from Hudgins et al. 2022 (J Appl Ecol)
ashes<-subset(pred_trees_huge,genus=="fraxinus")
#costs from Aukema et al. 2011 corrected to 2019 USD
small_st_treat<-450*(117.244/103.157)
med_st_treat<-600*(117.244/103.157)
large_st_treat<-1200*(117.244/103.157)

#weighted average cost by tree size across all communities modelled in Hudgins et al. 2022
(sum(ashes$pred_small)/sum(ashes$pred_large+ashes$pred_small+ashes$pred_med)*small_st_treat+sum(ashes$pred_med)/sum(ashes$pred_large+ashes$pred_small+ashes$pred_med)*med_st_treat+sum(ashes$pred_large)/sum(ashes$pred_large+ashes$pred_small+ashes$pred_med)*large_st_treat)*(1239440-366020)
