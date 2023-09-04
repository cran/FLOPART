# FLOPART: Functional Labeled Optimal Partitioning

[https://github.com/tdhock/FLOPART/workflows/R-CMD-check/badge.svg](https://github.com/tdhock/FLOPART/actions)

Functional Labeled Optimal Partitioning (FLOPART) is an optimal peak
detection algorithm with label constraints. The dynamic programming
algorithm computes a segmentation and corresponding set of peaks which
minimizes the penalized Poisson loss (for non-negative integer count
data), subject to the following constraints:

* SAME as previous PeakSeg packages: there are two states (up/peak and
  down/background). Changes from up/peak state must be non-increasing
  and go to down/background state. Changes from down/background state
  must be non-decreasing and go to up/peak state.
* LESS constraints relative to previous PeakSeg packages: there are no
  constraints for the model to start and end in the down/background
  state.
* MORE constraints relative to previous PeakSeg packages: there are
  additional constraints defined by labels (see below).

# Installation

```r
if(require("remotes"))install.packages("remotes")
remotes::install_github("tdhock/FLOPART")
```

# Usage

The main driver function is FLOPART, which takes three arguments
* `coverage` is a data table with columns chromStart, chromEnd, count.
* `label` is a data table with columns chromStart, chromEnd, annotation.
* `penalty` is a non-negative numeric value, larger for fewer changepoints/peaks.

* chromStart and chromEnd columns are integer positions on a chromosome. 
* count is an integer data value to segment.
* annotation is a character string, one of noPeaks, peakStart, peakEnd.

The constraints for the different label types are as follows: 

* noPeaks: the model must stay in the down/background state in this label.
* peakStart: at the start of the label the model must be in the
  down/background state, there must be exactly one non-decreasing
  change in this label, and the model must be in the up/peak state at
  the end of the label.
* peakEnd: at the start of the label the model must be in the up/peak
  state, there must be exactly one non-increasing change in this
  label, and the model must be in the down/background state at the end
  of the label.
  
For a simple real data example, 

```r
> library(data.table) # for print method.
> data("Mono27ac.simple", package="FLOPART")
> Mono27ac.simple
$coverage
      chrom chromStart chromEnd count
   1: chr11     145000   146765     0
   2: chr11     146765   146807     1
   3: chr11     146807   175254     0
   4: chr11     175254   175296     1
   5: chr11     175296   175738     0
  ---                                
2975: chr11     326980   326981     5
2976: chr11     326981   326983     6
2977: chr11     326983   326985     5
2978: chr11     326985   326992     4
2979: chr11     326992   327000     3

$label
   chrom chromStart chromEnd annotation
1: chr11     180000   200000    noPeaks
2: chr11     208000   220000    peakEnd
3: chr11     300000   308250  peakStart
4: chr11     308260   320000    peakEnd
```

To run FLOPART on these data,

```r
fit <- with(Mono27ac.simple, FLOPART::FLOPART(coverage, label, penalty=1400))
```

The `fit` is a list of results; we can use `head` to look at the first
few items in each of these results.

```r
> lapply(fit, head)
$coverage_dt
   chromStart chromEnd count weight
1:     145000   146765     0   1765
2:     146765   146807     1     42
3:     146807   175254     0  28447
4:     175254   175296     1     42
5:     175296   175738     0    442
6:     175738   175780     1     42

$label_dt
   chromStart chromEnd annotation type firstRow lastRow
1:     180000   200000    noPeaks    0       14     118
2:     208000   220000    peakEnd   -1      725    1322
3:     300000   308250  peakStart    1     2723    2770
4:     308260   320000    peakEnd   -1     2773    2871

$cost_mat
           [,1]       [,2]
[1,] 0.00000000 0.00000000
[2,] 0.11067717 0.11067717
[3,] 0.01052251 0.01052251
[4,] 0.01909784 0.01909784
[5,] 0.01886280 0.01886280
[6,] 0.02660139 0.02660139

$intervals_mat
     [,1] [,2]
[1,]    1    1
[2,]    2    1
[3,]    3    3
[4,]    3    5
[5,]    4    4
[6,]    4    5

$segments_dt
          mean firstRow lastRow state chromStart chromEnd     status
1:  0.06556501        1     197     0     145000   206725 background
2: 13.17743878      198    1132  2986     206725   209216       peak
3:  0.22431609     1133    1423     0     209216   236120 background
4:  7.38781362     1424    1792  2986     236120   237515       peak
5:  0.19459495     1793    2079     0     237515   267598 background
6:  3.68970814     2080    2552  2986     267598   270853       peak
```

The main output of the algorithm is `segments_dt` (shown above) which
is a data table of optimal segment means, given the penalty and
constraints. Each row of this data table describes a segment in the
optimal model:

* mean is the mean value (in units of the count column from the coverage data)
* chromStart/chromEnd give the start/end positions of each segment, in
  units of chromStart/chromEnd from the coverage
  data. firstRow/lastRow are similar but in units of data points (rows
  in coverage_dt).
* state/status show if the segment is in peak or background.

Other outputs which could be of interest:

* `cost_mat` is the optimal cost computed by dynamic programming at
  each data point (row) in each state (column).
* `intervals_mat` is the number of intervals (function pieces) stored
  in the piecewise cost function at each data point (row) in each
  state (column). This is useful for empirical analysis of time
  complexity, since the computation time is linear in the number of
  intervals.
* `coverage_dt` is a modified version of the initial coverage data
  set. For example there will be additional rows if there are labels
  with start/end values that are not present in the initial coverage.
* `label_dt` has additional columns firstRow/lastRow which correspond
  to the start/end positions of the labels, in units of the rows in
  `coverage_dt`.

# Related work

* [PeakSegOptimal](https://github.com/tdhock/PeakSegOptimal) and
  [PeakSegDisk](https://github.com/tdhock/PeakSegDisk) implement
  solvers for the Poisson model without label constraints, but with an
  additional constraint that the model must start and end in the
  down/background state.
* [LOPART](https://github.com/tdhock/LOPART) implements a solver for
  the normal model with 0/1 label constraints, but no inequality
  constraints (each change can be in any direction, up or down).
* [Repo with figures from related paper in progress](https://github.com/tdhock/LabeledFPOP-paper).
