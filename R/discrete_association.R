#' Performs t-tests on a matrix of binary features and a vector y.
#'
#' @param X n x m matrix of binary features (missing values are allowed).
#' @param y length n vector of numerical values (missing values are allowed).
#'
#' @return A table with the following columns:
#' \describe{
#' feature: the column names of X. \cr
#' effect_size: an estimate of the difference between groups. \cr
#' t.stat: the t-statistic associated with the test on group differences. \cr
#' p.value: the p-value of the t-statistic. \cr
#' q.value: the multiple hypothesis corrected p-value. \cr
#' }
#' 
#' @export discrete_test
discrete_test <- function(X, y) {
  
  y <- y[is.finite(y)]  # only finite values
  X <- X[, apply(X, 2, function(x) !any(is.na(x)))]
  
  overlap <- dplyr::intersect(rownames(X), names(y))
  
  if (length(overlap) < 10) {
    stop("Not enough overlapping data")
  }
  
  y <- y[overlap]
  X <- X[overlap,]
  
  out_table <- run_lm_stats_limma(X, y, target_type = "feature")
  out_table %<>%
    dplyr::select(feature, EffectSize, t_stat, p.value, q.value, ) %>%
    dplyr::rename(effect_size = EffectSize)
  
  return(out_table)
}
