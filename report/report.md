
# European Winter Temperatures and the Quasi-Biennial Oscillation

## Introduction

- mechanisms that influence European winter temperature variability
- past advances in applying complex networks to similar issues

## Methods

### Data Preparation

- removal of linear background (e.g. climate change)
- normalize to mean zero and variance one
- downsampling to daily or weekly resolution
- (denoising)
- keep only winter temperatures between 01.12 and 01.03
- (band pass filtering?) show fourier spectrum

### Construction of the Complex Network

- pearson correlation coefficient
- (mutual information)
- RMD from joint recurrence plot
- correlation matrix

- for each pair of time series (X, Y)
- construct S surrogates for Y
- calculate the synchronization measure between X and all surrogates to get a distribution
- calculate the synchronization measure between X and the original Y
- construct a 95% confidence interval and reject the original synchronization measure if it lies outside
- create a significance matrix by calculating the distance to the mean of the test distribution

- goal is to ge ta p=5% link density network (Donges)
- for this take the upper (1 - p) quantile of the significance distribution (of the XY pairs)
- connect all pairs of X and Y that fall in this quantile with links
- this implies that the upper (1 - p) quantile of the significance distribution only contains significant values
- WHAT TO DO WHEN THIS IS NOT THE CASE OR WHEN THERE ARE MORE SIGNIFICANT LINKS?
- TELECONNECTIONS SHOULD BE INCLUDED? HOW TO DO THAT?
- The resulting network has thus only significant connections and a link density of 5%

### Correction for Background Effects

- surrogate networks should share the same link density and link distance distribution
- calculate a distance matrix between all XY pairs using great circle distance of the lat/lon 
  coordinates (faster than geodesic)
- calculate the link distance distribution by looking the distance up in the distance matrix 
- for each possible link calculate the possibility p_link(d) of connection as the ratio of seen links with this distance 
  over the number of possible links with this distance
- connect the links randomly with probability p_link(d)

- construct NS network surrogates for which you calculate a metric (centrality, betweenness)
- remove the mean of this distribution from the original network's metric and estimate the error using the variance

### Evolving Network Analysis

- break down the whole time series into windows (say 20, depends on calculation time)
- generate networks for each window and calculate network measures 
- correct again for boundary effects


