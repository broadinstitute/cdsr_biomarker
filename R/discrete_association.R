require(tidyverse)
require(magrittr)
require(stats)
require(data.table)

#' Performs t-tests on a matrix of binary features and a vector y.
#'
#' @param X n x m matrix of binary features (missing values are allowed).
#' @param y length n vector of numerical values (missing values are allowed).
#'
#' @return A table with the following columns:
#' feature: the column names of X.
#' effect_size: an estimate of the difference between groups.
#' t.stat: the t-statistic associated with the test on group differences.
#' p: the p-value of the t-statistic.
#' q: the multiple hypothesis corrected p-value.
#' 
#' @export
#'
discrete_test <- function(X, y) {
  y <- y[is.finite(y)]  # only finite values
  X <- X[, apply(X, 2, function(x) !any(is.na(x)))]
  
  out_table <- tibble::tibble()  # tracks output
  feats <- colnames(X)  # features to run tests on
  overlap <- dplyr::intersect(rownames(X), names(y))
  
  if (length(overlap) < 10) {
    stop("Not enough overlapping data")
  }
  
  y <- y[overlap]
  X <- X[overlap,]
  
  # loop through each feature (e.g. lineage)
  for (feat in feats) {
    has_feat <- rownames(X[X[,feat] == 1,])
    group <- y[has_feat]
    others <- dplyr::setdiff(y, group)
    
    # need more than one data point to generate a group mean and test
    if (length(group) >= 5 & length(others) >= 5) {
      
      # get t test results (two-sample unpaired)
      t <- stats::t.test(group, others)
      
      # compile results and add to final table
      result <- tibble::tibble(feature = feat,
                               effect_size = t$estimate[1] - t$estimate[2],
                               t.stat = t$statistic["t"], p = t$p.value)
      out_table %<>%
        dplyr::bind_rows(result)
    }
  }
  
  if (nrow(out_table) < 1) {
    stop("All group sizes too small")
  }
  
  # calculate q values (BH procedure)
  out_table$q <- stats::p.adjust(c(out_table$p,
                                   rep(1, length(feats) - nrow(out_table))),
                                 method = "BH")[1:nrow(out_table)]
  
  return(out_table)
}
