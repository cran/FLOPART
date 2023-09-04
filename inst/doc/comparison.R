## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(data.table)
data("Mono27ac.simple", package="FLOPART")
Mono27ac.simple

## -----------------------------------------------------------------------------
ann.colors <- c(
  noPeaks="orange",
  peakStart="#efafaf",
  peakEnd="#ff4c4c")
if(require("ggplot2")){
  ggplot()+
    theme_bw()+
    scale_fill_manual("label", values=ann.colors)+
    geom_rect(aes(
      xmin=chromStart-0.5, xmax=chromEnd+0.5,
      ymin=-Inf, ymax=Inf,
      fill=annotation),
      alpha=0.5,
      color="grey",
      data=Mono27ac.simple[["label"]])+
    geom_step(aes(
      chromStart, count),
      data=Mono27ac.simple[["coverage"]],
      color="grey50")
}

## -----------------------------------------------------------------------------
label.pen <- 1400
fit <- with(Mono27ac.simple, FLOPART::FLOPART(coverage, label, label.pen))
lapply(fit, head)
str(fit)

## -----------------------------------------------------------------------------
penalty.vec <- c(
  "7"=1400,
  "6"=1450,
  "5"=1460)
unlabeled.segs.dt.list <- list()
for(penalty in penalty.vec){
  ufit <- FLOPART::FLOPART(Mono27ac.simple[["coverage"]], penalty=penalty)
  pen.segs <- ufit[["segments_dt"]]
  pen.peaks <- pen.segs[status=="peak"]
  n.peaks <- nrow(pen.peaks)
  unlabeled.segs.dt.list[[paste(penalty)]] <- data.table(
    penalty, n.peaks, pen.segs)
}
(unlabeled.segs.dt <- do.call(rbind, unlabeled.segs.dt.list))

## -----------------------------------------------------------------------------
model.color <- "blue"
if(require("ggplot2")){
  ggplot()+
    ggtitle("Models without label constraints")+
    geom_step(aes(
      chromStart, count),
      data=Mono27ac.simple[["coverage"]],
      color="grey50")+
    geom_step(aes(
      chromStart, mean),
      data=unlabeled.segs.dt,
      color=model.color)+
    theme_bw()+
    theme(panel.spacing=grid::unit(0, "lines"))+
    facet_grid(penalty + n.peaks ~ ., labeller=label_both)
}

## -----------------------------------------------------------------------------
cons.segs <- fit[["segments_dt"]]
model.list <- c(unlabeled.segs.dt.list, list(FLOPART=data.table(
  penalty=label.pen, n.peaks=sum(cons.segs[["status"]]=="peak"), cons.segs)))
peaks.dt.list <- list()
err.dt.list <- list()
segs.dt.list <- list()
for(model in names(model.list)){
  model.segs <- model.list[[model]]
  model.peaks <- model.segs[status=="peak"]
  model.err <- if(requireNamespace("PeakError")){
    perr <- PeakError::PeakErrorChrom(model.peaks, Mono27ac.simple[["label"]])
    data.table(perr)[, .(chromStart, chromEnd, status)]
  }else{
    data.table(chromStart=integer(), chromEnd=integer(), status=character())
  }
  segs.dt.list[[model]] <- data.table(model, model.segs)
  err.dt.list[[model]] <- data.table(model, model.err)
  peaks.dt.list[[model]] <- data.table(model, model.peaks)
}
(segs.dt <- do.call(rbind, segs.dt.list))
(err.dt <- do.call(rbind, err.dt.list))
(peaks.dt <- do.call(rbind, peaks.dt.list))

## -----------------------------------------------------------------------------
peak.y <- -2
if(require("ggplot2")){
  ggplot()+
    ggtitle("Models with label constraints (FLOPART) and without (penalty values)")+
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
      data=segs.dt,
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
    theme(panel.spacing=grid::unit(0, "lines"))+
    facet_grid(model ~ ., labeller=label_both)+
    scale_linetype_manual(
      "error type",
      values=c(
        correct=0,
        "false negative"=3,
        "false positive"=1))+
    geom_rect(aes(
      xmin=chromStart, xmax=chromEnd,
      ymin=-Inf, ymax=Inf,
      linetype=status),
      fill=NA,
      color="black",
      size=1,
      data=err.dt)
}

