# Script to tune treeClust (and RPART) hyperparameters using the REACTOME test set (area under the precision-recall curve, AUPRC, as the final read-out).
# 2018 P. Grabowski, Rappsilber lab

library(data.table)
library(treeClust)
library(reshape2)
library(ROCR)
library(MESS)
set.seed(42)
setwd("./")

auc2 = function (x, y, from = min(x), to = max(x), type = c("linear","spline"), absolutearea = FALSE, ...){
  type = match.arg(type)
  if (length(x) != length(y))
    stop("x and y must have the same length")
  if (length(unique(x)) < 2)
    return(NA)
  if (type == "linear") {
    if (absolutearea)
      y = y - min(y)
    values = approx(x, y, xout = sort(unique(c(from, to,
                                                x[x > from & x < to]))), ...)
    res = 0.5 * sum(diff(values$x) * (values$y[-1] + values$y[-length(values$y)]))
    if (absolutearea)
      res = res - min(y) * (max(x) - min(x))
  }
  else {
    if (absolutearea)
      myfunction = function(z) {
        abs(splinefun(x, y, method = "natural")(z))
      }
    else myfunction = splinefun(x, y, method = "natural")
    res = integrate(myfunction, lower = from, upper = to, subdivisions = 500L, stop.on.error = FALSE, rel.tol = 0.00005 + .Machine$double.eps^0.25)$value
  }
  res
}

## A function to calculate AUPRC using the bench (REACTOME) data
calc_AUPRC = function(tc_dists){
  tc_sims = 1 - as.matrix(tc_dists)
  tc_sims = data.table(melt(tc_sims))
  colnames(tc_sims) = c("Protein_1","Protein_2","treeClust_sim")
  tc_sims = tc_sims[!is.na(treeClust_sim),]
  setkey(tc_sims, Protein_1, Protein_2)

  ## Calculate AUPRC
  merge_bench = merge(bench, tc_sims)
  merge_bench = merge_bench[,-grep("Protein",colnames(merge_bench)),with=F]
  labels = merge_bench$Class
  merge_bench$Class = NULL
  pred = prediction( predictions = merge_bench$treeClust_sim,
                     labels = labels,
                     label.ordering = c("FP", "TP"))
  perf = performance(pred, measure = "prec", x.measure = "rec")
  precision = perf@y.values[[1]][500:length(perf@y.values[[1]])]
  recall = perf@x.values[[1]][500:length(perf@x.values[[1]])]
  auprc = auc2(recall,precision, type="spline")
  auprc
}

## Load saved benchmark data (down-sampled to 10% TPs)
bench = fread("Reactome_TP_FP_10perc_subset_for_GS.csv")
setkey(bench, Protein_1, Protein_2)

## Load median-normalized ProteomeHD data
prohd = fread("ProteomeHD_v1_1.csv")
prohd = prohd[!duplicated(prohd$Simplified_protein_ID),]
prohd = data.frame(prohd, row.names = 2)
prohd = prohd[,-c(1:3)]

## Select only proteins with at least 95 ratios
feature_count = apply(prohd, 1, function(x){ sum(!is.na(x))})
prohd_ratios_min95 = prohd[feature_count >= 95,]

## Run the grid search
for( serule_i in seq(1,2,by=0.2)){
  for( cp_i in seq(0.005,0.15,by=0.01)){
    print("Running hyperparam search...")
    tc_dist_temp = treeClust.dist( prohd_ratios_min95,
                                   d.num = 2,
                                   verbose = TRUE,
                                   rcontrol = rpart.control(cp = cp_i),
                                   control = treeClust.control(parallelnodes = 2, serule = serule_i) )
    curr_auc = calc_AUPRC(tc_dist_temp)
    print("AUC:")
    print(curr_auc)
  }
}
