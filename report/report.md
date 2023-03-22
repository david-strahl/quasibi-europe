
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

To build a complex network from the gridded reanalysis data, I considered each grid-cell as a node within the network.
For a grid of N?? by M?? cells, this gives a total of NxM nodes. A link within the network represent an effect of 
one node on another. Therefore, no self loops are allowed within the network and links are non-directional. This 
reduces the number of total links to
N*(N - 1)/2

To determine whether a pair of nodes will be connected in the network I look at the dependence between time series 
of the corresponding grid cells. This dependence may either be quantified using the linear Pearson correlation or 
non-linear measures like the mutual information of the recurrence based measure of dependence. 



To calculate the recurrence based measure of dependence of two time series x and y one firstly 
builds recurrence plots by thresholding the distance of visited in phase space:
R_ij^^ = Heaviside(|xi - xj| - eps)
where |.| is an appropriate distance metric (here absolute difference) and eps is a distance threshold. The mean 
over i is called the column recurrence rate while the mean over i and j is called the (total) recurrence rate. The 
threshold may also be selected by 


The constructed network should of course only include correlations that are statistically significant. Both (Donges) 
and (Goswami) described a simple procedure than can be used to test whether the correlation between two grid-cells is 
statistically significant. It builds on the comparison against a test distribution derived from surrogate time series.
Assume you have two time series x and y. Firstly, calculate the correlation measure between them. Secondly, generate 
an appropriate number of surrogates for y. Here, I used twin surrogates as described by (Goswani). Thirdly, calculate 
the correlation measure between x and all surrogates of y build a test distribution. Finally, construct a 95% 
confidence interval around the mean of the test distribution and reject all values that fall within this interval as 
statistically insignificant. The difference to the mean in units of standard deviation can be used to rank statistical 
significance of all links, where of course values that are insignificant are disregarded.

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


