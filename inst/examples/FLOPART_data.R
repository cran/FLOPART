library(data.table)
d <- function(chromStart, chromEnd, count){
  data.table(chromStart, chromEnd, count)
}
(cov.dt <- rbind(
  d(0, 10, 53),
  d(10, 20, 124)))
l <- function(chromStart, chromEnd, annotation){
  data.table(chromStart, chromEnd, annotation)
}
lab.dt <- rbind(
  l(2, 7, "noPeaks"),
  l(8, 15, "peakStart"))
FLOPART::FLOPART_data(cov.dt)
FLOPART::FLOPART_data(cov.dt, lab.dt)
data("Mono27ac", package="FLOPART")
sapply(Mono27ac, dim)
converted <- with(Mono27ac, FLOPART::FLOPART_data(coverage, labels))
sapply(converted, dim)
