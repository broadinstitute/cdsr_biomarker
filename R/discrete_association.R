#' Estimate linear-model stats for a matrix of data using limma with empirical Bayes moderated t-stats for p-values
#'
#' @param mat: Nxp data matrix with N cell lines and p genes
#' @param vec: N vector of independent variables. Can be two-group labels as factors, bools, or can be numeric
#' @param covars: Optional Nxk matrix of covariates
#' @param weights: Optional N vector of precision weights for each data point
#' @param target_type: Name of the column variable in the data (default 'Gene')
#' @param limma_trend: Whether to fit an intensity trend with the empirical Bayes variance model
#'
#' @return: data frame of stats
#' @export
#'
#' @examples
#' CRISPR = load.from.taiga(data.name='avana-2-0-1-d98f',
#' data.version=1,
#' data.file='ceres_gene_effects',
#' transpose = T)
#' is_panc <- load.from.taiga(data.name = 'ccle-lines-lineages') %>% .[, 'pancreas']
#' ulines <- intersect(rownames(CRISPR), names(is_panc))
#' lim_res <- run_lm_stats_limma(CRISPR[ulines,], is_panc[ulines])
#' @export run_lm_stats_limma
run_lm_stats_limma <- function(mat, vec, covars = NULL, weights = NULL,
                               target_type = 'Gene', limma_trend = FALSE) {
  require(limma)
  require(magrittr)
  require(tibble)
  require(plyr)
  require(dplyr)
  
  udata <- which(!is.na(vec))
  if (!is.numeric(vec)) {
    pred <- factor(vec[udata])
    stopifnot(length(levels(pred)) == 2) #only two group comparisons implemented so far
    n_out <- colSums(!is.na(mat[udata[pred == levels(pred)[1]],,drop=F]))
    n_in <- colSums(!is.na(mat[udata[pred == levels(pred)[2]],,drop=F]))
    min_samples <- pmin(n_out, n_in) %>% set_names(colnames(mat))
  } else {
    pred <- vec[udata]
    min_samples <- colSums(!is.na(mat[udata,]))
  }
  #there must be more than one unique value of the independent variable
  if (length(unique(pred)) <= 1) {
    return(NULL)
  }
  #if using covariates add them as additional predictors to the model
  if (!is.null(covars)) {
    if (!is.data.frame(covars)) {
      covars <- data.frame(covars)
    }
    combined <- covars[udata,, drop = FALSE]
    combined[['pred']] <- pred
    form <- as.formula(paste('~', paste0(colnames(combined), collapse = ' + ')))
    design <- model.matrix(form, combined)
    design <- design[, colSums(design) != 0, drop = FALSE]
  } else {
    design <- model.matrix(~pred)
  }
  if (!is.null(weights)) {
    if (is.matrix(weights)) {
      weights <- t(weights[udata,])
    } else{
      weights <- weights[udata]
    }
  }
  fit <- limma::lmFit(t(mat[udata,]), design, weights = weights)
  fit <- limma::eBayes(fit, trend = limma_trend)
  targ_coef <- grep('pred', colnames(design), value = TRUE)
  results <- limma::topTable(fit, coef = targ_coef, number = Inf)
  
  if (colnames(results)[1] == 'ID') {
    colnames(results)[1] <- target_type
  } else {
    results %<>% rownames_to_column(var = target_type)
  }
  results$min_samples <- min_samples[results[[target_type]]]
  
  two_to_one_sided <- function(two_sided_p, stat, test_dir) {
    #helper function for converting two-sided p-values to one-sided p-values
    one_sided_p <- two_sided_p / 2
    if (test_dir == 'right') {
      one_sided_p[stat < 0] <- 1 - one_sided_p[stat < 0]
    } else {
      one_sided_p[stat > 0] <- 1 - one_sided_p[stat > 0]
    }
    return(one_sided_p)
  }
  results %<>% set_colnames(revalue(colnames(.), c('logFC' = 'EffectSize', 'AveExpr' = 'Avg', 't' = 't_stat', 'B' = 'log_odds',
                                                   'P.Value' = 'p.value', 'adj.P.Val' = 'q.value', 'min_samples' = 'min_samples'))) %>% na.omit()
  results %<>% dplyr::mutate(p.left = two_to_one_sided(p.value, EffectSize, 'left'),
                             p.right = two_to_one_sided(p.value, EffectSize, 'right'),
                             q.left = p.adjust(p.left, method = 'BH'),
                             q.right = p.adjust(p.right, method = 'BH'))
  return(results)
}

#' Performs t-tests on a matrix of binary features and a vector y.
#'
#' @param X n x m matrix of binary features (missing values are allowed).
#' @param y length n vector of numerical values (missing values are allowed).
#'
#' @return A table with the following columns:
#' feature: the column names of X.
#' effect_size: an estimate of the difference between groups.
#' t.stat: the t-statistic associated with the test on group differences.
#' p.value: the p-value of the t-statistic.
#' q.value: the multiple hypothesis corrected p-value.
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
