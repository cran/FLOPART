label_colors <- c(
  noPeaks="#f6f4bf",
  peakStart="#ffafaf",
  peakEnd="#ff4c4c")

##' FLOPART needs at most one label per coverage data row, which may
##' not be the case for arbitrary coverage/labels.
##' @title Convert data for input to FLOPART
##' @param coverage data frame of coverage with columns chromStart,
##'   chromEnd, count
##' @param label data frame of labels with with columns chromStart,
##'   chromEnd, annotation
##' @return named list: coverage_dt is data table representing a
##'   run-length encoding of the input coverage data, with additional
##'   rows if there are label chromStart/chromEnd values not in the
##'   set of coverage positions; label_dt is a data table with one row
##'   per label, and additional columns firstRow/lastRow which refer
##'   to row numbers of coverage_dt, 0-based for passing to C++ code.
##' @author Toby Dylan Hocking
##' @example inst/examples/FLOPART_data.R
FLOPART_data <- function(coverage, label){
  i.chromStart <- i.chromEnd <- count <- chromStart <- chromEnd <-
    annotation <- labelStart <- labelEnd <- type <- firstRow <-
      lastRow <- . <- run.i <- NULL
  if(missing(label)){
    label <- data.table(
      chromStart=integer(), chromEnd=integer(), annotation=character())
  }
  if(!is.data.frame(coverage)){
    stop("coverage must be a data frame")
  }
  if(any(max(coverage[["chromEnd"]]) < label[["chromEnd"]])){
    stop("label ends must be on or before last coverage")
  }
  if(any(label[["chromStart"]] < min(coverage[["chromStart"]]))){
    stop("label starts must be on or after first coverage")
  }
  label_code <- get_label_code()
  cov_dt <- data.table(coverage)
  cov_dt[, run.i := cumsum(c(1, diff(count)!=0))]
  stop_count_data_missing <- function(){
    stop("count data missing, meaning that some chromStart are not equal to previous chromEnd, please fix by adding rows with count=0")
  }
  compressed <- cov_dt[, {
    if(any(chromStart[-1] != chromEnd[-.N])){
      stop_count_data_missing()
    }
    .(chromStart=chromStart[1], chromEnd=chromEnd[.N])
  }, by=.(run.i, count)]
  df.names <- list("label", "compressed")
  key.vec <- c("chromStart", "chromEnd")
  uniq.pos.list <- list()
  dt.list <- list()
  for(data.type in df.names){
    data.dt <- data.table(get(data.type), key=key.vec)
    for(k in key.vec){
      uniq.pos.list[[paste(data.type, k)]] <- data.dt[[k]]
    }
    dt.list[[data.type]] <- data.dt
  }
  uniq.pos <- unique(sort(unlist(uniq.pos.list)))
  uniq.dt <- data.table(
    chromStart=uniq.pos[-length(uniq.pos)],
    chromEnd=uniq.pos[-1])
  with.counts <- dt.list[["compressed"]][
    uniq.dt,
    .(chromStart=i.chromStart, chromEnd=i.chromEnd, count,
      weight = i.chromEnd - i.chromStart),
    on=.(chromStart < chromEnd, chromEnd > chromStart)]
  if(any(is.na(with.counts$count))){
    stop_count_data_missing()
  }
  label.index.dt <- with.counts[
    dt.list[["label"]],
    .(firstRow=.I[1]-1L,
      lastRow=.I[.N]-1L,
      labelStart=i.chromStart,
      labelEnd=i.chromEnd,
      type=label_code[paste(annotation)],
      annotation),
    by=.EACHI,
    on=.(chromEnd > chromStart, chromStart < chromEnd)
  ][, .(
    chromStart=labelStart,
    chromEnd=labelEnd,
    annotation, type,
    firstRow, lastRow
  )]
  list(coverage_dt=with.counts, label_dt=label.index.dt)
}

##' Main function for computing optimal segmentation model with
##' Poisson loss, up-down constraints, and label constraints.
##' @title Functional Labeled Optimal Partitioning
##' @param coverage data frame of coverage
##' @param label data frame of labels
##' @param penalty non-negative penalty constant
##' @return list with named elements: coverage_dt is a data table with
##'   columns chromStart, chromEnd, count, weight; label_dt is a data
##'   table with columns chromStart, chromEnd, annotation, type,
##'   firstRow, lastRow; cost_mat is a Nx2 numeric matrix of optimal
##'   penalized Poisson loss values up to each data point and in each
##'   state; intervals_mat is a Nx2 integer matrix of counts of
##'   intervals used to store the optimal cost function, useful for
##'   analyzing time/space complexity; segments_dt is a data table
##'   with columns chromStart, chromEnd, status, mean.
##' @author Toby Dylan Hocking
##' @example inst/examples/FLOPART.R
FLOPART <- function(coverage, label, penalty){
  status <- state <- . <- chromStart <- chromEnd <- NULL
  data.list <- FLOPART_data(coverage, label)
  result <- with(data.list, FLOPART_interface(
    coverage_dt[["count"]],
    coverage_dt[["weight"]],
    penalty,
    label_dt[["type"]],
    label_dt[["firstRow"]],
    label_dt[["lastRow"]]))
  segs <- as.data.table(result[["segments_df"]])
  if(nrow(segs)==0)warning("there is no feasible model given label constraints; fix by modifying labels")
  name.vec <- c(firstRow="chromStart", lastRow="chromEnd")
  for(xRow in names(name.vec)){
    chromX <- name.vec[[xRow]]
    row.vec <- segs[[xRow]]
    pos.vec <- data.list[["coverage_dt"]][[chromX]]
    set(segs, j=chromX, value=pos.vec[row.vec])
  }
  segs[, status := ifelse(state==0, "background", "peak")]
  result[["segments_df"]] <- NULL
  result[["segments_dt"]] <- segs[, .(
    chromStart, chromEnd, status, mean)]
  c(data.list, result)
}
