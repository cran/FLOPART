% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FLOPART.R
\name{FLOPART}
\alias{FLOPART}
\title{Functional Labeled Optimal Partitioning}
\usage{
FLOPART(coverage, label, penalty)
}
\arguments{
\item{coverage}{data frame of coverage}

\item{label}{data frame of labels}

\item{penalty}{non-negative penalty constant}
}
\value{
list with named elements: coverage_dt is a data table with
  columns chromStart, chromEnd, count, weight; label_dt is a data
  table with columns chromStart, chromEnd, annotation, type,
  firstRow, lastRow; cost_mat is a Nx2 numeric matrix of optimal
  penalized Poisson loss values up to each data point and in each
  state; intervals_mat is a Nx2 integer matrix of counts of
  intervals used to store the optimal cost function, useful for
  analyzing time/space complexity; segments_dt is a data table
  with columns chromStart, chromEnd, status, mean.
}
\description{
Main function for computing optimal segmentation model with
Poisson loss, up-down constraints, and label constraints.
}
\examples{
library(data.table)
data("Mono27ac.simple", package="FLOPART")
Mono27ac.simple
label.pen <- 1400
fit <- with(Mono27ac.simple, FLOPART::FLOPART(coverage, label, label.pen))
lapply(fit, head)

## Plot data and model.
ann.colors <- c(
  noPeaks="orange",
  peakStart="#efafaf",
  peakEnd="#ff4c4c")
model.color <- "blue"
(peaks.dt <- fit[["segments_dt"]][status=="peak"][, peak.y := -2][])
if(require("ggplot2")){
  ggplot()+
    ggtitle("Model with label constraints (FLOPART)")+
    scale_fill_manual("label", values=ann.colors)+
    geom_rect(aes(
      xmin=chromStart, xmax=chromEnd,
      ymin=-Inf, ymax=Inf,
      fill=annotation),
      alpha=0.5,
      color="grey",
      data=Mono27ac.simple[["label"]])+
    geom_step(aes(
      chromStart, count),
      data=Mono27ac.simple[["coverage"]],
      color="grey50")+
    geom_step(aes(
      chromStart, mean),
      data=fit[["segments_dt"]],
      color=model.color)+
    geom_segment(aes(
      chromStart, peak.y,
      xend=chromEnd, yend=peak.y),
      color=model.color,
      size=1,
      data=peaks.dt)+
    geom_point(aes(
      chromEnd, peak.y),
      color=model.color,
      shape=21,
      fill="white",
      data=peaks.dt)+
    theme_bw()+
    theme(panel.spacing=grid::unit(0, "lines"))
}

## To analyze computational complexity, plot number of intervals
## stored in cost function, versus data point, for each cost status.
imat <- fit[["intervals_mat"]]
interval.dt <- data.table(
  intervals=as.integer(imat),
  status=c("peak", "background")[as.integer(col(imat))],
  data.i=as.integer(row(imat)))
if(require("ggplot2")){
  ggplot()+
    scale_fill_manual("label", values=ann.colors)+
    geom_rect(aes(
      xmin=firstRow-0.5, xmax=lastRow+0.5,
      ymin=-Inf, ymax=Inf,
      fill=annotation),
      alpha=0.5,
      color="grey",
      data=fit[["label_dt"]])+
    geom_line(aes(
      data.i, intervals, color=status),
      size=1,
      data=interval.dt)
}
}
\author{
Toby Dylan Hocking
}
