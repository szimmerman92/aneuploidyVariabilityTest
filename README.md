# aneuploidyVariabilityTest
Create a test to see if chromosome 18 is duplicated or removed more often than other chromosomes.

## To do

1. Change resampling method to *bootstrap* 
	- take N samples *with replacement*, where N is the the number of cells.
	- 1000 bootstraps
	
2. Compare bootstrap distribution between each pair of chromosomes
	- Wilcoxon or t-test to assess for similarity in variability between pairs of chromosomes
	- Possibly KS test as well
	- Generate heatmap with p-values (23x23 matrix)
	- 

3. Linear regression - identify features that predict chromosome variability
	- Do this by chromosome and by bin
	- Example features include # of genes, size of chromosome/segment, # of TFs in segment, gene density, open chromatin regions
	- Y is the observed frequency, Xs are the observed features

4. Visualization of p-values - use all pairwise comparisons. Plot the number of significantly different chromosomes on the x axis.

5. Network - assess similarity between pairwise comparisons.

6. Use association rules to determine multiple regions that are co-deleted or co-amplified
	- See "mining of massive datasets," Chapter 6.1
