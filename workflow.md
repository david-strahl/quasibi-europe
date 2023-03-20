# Workflow

## Data Cleanup

OPTIONAL:
- (Denoising using frequency filter below 12 month period and above 34 periods.)
- Remove linear trend if still existent.

ALWAYS:
- Normalize to zero mean and variance.
- Downsample to daily and weekly resolution by averaging.
- Subset only the winter part (21 Dez - 21 March).

Keep all combinations as data, e.g.:
- ``denoised_detrended_daily_winter.npz``
- ``weekly_winter.npz``

## Construction of Stationary Correlation Network

- Joint-Recurrence and MI or Pearson based network (for comparison) of whole timeseries.
  - What embedding dimension for the state space?
- Store the whole adjecency matrix.
- Threshold the adjecency matrix to build the network.
  - Only maintain statistically significant connections (relative to shuffled or surrogate timeseries)
  - Teleconnections should be included.
  - For simplicity you could go with the edge density 0.005 that is given by Donges and derive the threshold from it.
    All extensions to the method can then be made later.
- Calculate possible network metrics (degree centrality, AWC, local and global clustering coefficients, closeness 
  centrality, betweenness centrality, average pathlength) and store them to facilitate later analysis.

Questions that need to be answered here:
- How to do the Joint-Recurrence? Multiply the recurrence plots.
- Can you get PyUnicorn to work so that your coding workload is reduced? DONE

## Construction of Non-Stationary Correlation Network

- Split the time series into a reasonable amount of windows (DECIDE ON THIS SOON).
- Do the same stationary analysis for each window to find variations in the metrics.
- Calculate Hamming distance between windows. 

Questions that need to be answered here:
- Should windows overlap?
- How many windows are resonable (depends on runtime of stationary analysis per time step)

## Questions to Answer in the Report

- Can one construct a complex network of European winter temperatures based on join recurrences?
- Does the resulting network uncover features missed with the Pearson correlation coefficient?
- Can you generate an evolving network for European winter temperatures? 
- Do network metrics show significant changes in values?
- How can this guide further analysis? 

## Needed Plots
- Example of processed time series.
- Stationary JR vs. Pearson (vs. MI).
- Time series of interesting network measures.



## Takeaways from the Papers

### Donges
- Peasrson and MI networks show similar structures on the local and mesoscale and differ only slightly on the global 
  scale.
- Comparison should be done at equal edge densities and not thresholds for edge significance.
- Significance test against shuffled or surrogate timeseries.
- Spearman ranked rho for correlation between network metrics.