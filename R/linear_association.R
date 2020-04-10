require(tidyverse)
require(magrittr)
require(ashr)

#' Calculates the linear association between a matrix of features A and a vector y. 
#'
#' @param A n x m numerical matrix of features (missing values are allowed).
#' @param y length n vector of numerical values (missing values are not allowed)
#' @param W n x p numerical matrix of confounders (missing values are not allowed).
#' @param robust.se boolean flag to control if homoskedastic or heteroskedastic.
#' @param scale.A boolean to control if the columns of A to be scaled or not.
#' @param scale.y boolen to control if the variance of y should be scaled to 1.
#' @param is.dependent should y be considered as a dependent or independent variable
#' @param var.th float minimum variation of included feature
#' @param coef.var boolean to control if the coefficient of variation is used.
#'
#' @return A table with the following columns:
#' \describe{
#' feature:  column names of A \cr
#' p.val: p-values computed for the null-hypothesis beta = 0 \cr
#' q.val: q-values computed using Benjamini-Hotchberg method \cr
#' n: number of non-missing values for each column of A \cr
#' betahat: estimated linear coefficient controlled for W \cr
#' sebetahat: standard error for the estimate of beta \cr
#' sebetahat: standard error for the estimate of beta \cr
#' lfdr: local and global fdr values computed using adaptive shrinkage as \cr
#' qvalue: local and global fdr values computed using adaptive shrinkage as \cr
#' PosteriorMean: moderated effect size estimates again computed with adaptive shrinkage \cr
#' PosteriorSD: standard deviations of moderated effect size estimates \cr
#' z.score: z-score PosteriorMean/PosteriorSD}
#'
#' @export
#' 
lin_associations = function(A, y, W = NULL, robust.se = F, 
                            scale.A = T, scale.y = F, is.dependent = T, var.th = 0, coef.var = F){
  
  require(corpcor); require(tibble); require(dplyr); require(ashr)
  
  mean.na = function(x) mean(x, na.rm = T)
  
  N = length(y);
  y_ = scale(y, center = T, scale = scale.y)
  if (coef.var) {
    A_ = A[, apply(A,2,function(x) var(x,na.rm = T)/mean(x,na.rm = T)^2) > var.th]
    A_ = scale(A_, center = T, scale = scale.A)
  } else {
    A_ = A[, apply(A,2,function(x) var(x,na.rm = T)) > var.th]
    A_ = scale(A_, center = T, scale = scale.A)
  }
  N.A = apply(A_, 2, function(x) sum(is.finite(x)))
  
  
  if(is.null(W)) {
    W_ = matrix(rep(1,N))
  } else {
    W_ = cbind(W, rep(1,N))}
  
  
  P = diag(rep(1,N)) -  W_ %*% corpcor::pseudoinverse(t(W_) %*% W_) %*% t(W_)
  y_ = c(P %*% y_)
  
  A.NA = !is.finite(A_)
  A_ = A_; A_[A.NA] = 0
  A_ = P %*% A_; A_[A.NA] = NA

  
  if(is.dependent){
    K = apply(A_^2, 2, mean.na) 
    beta.hat =  apply(A_ * y_, 2, mean.na) / K; 
    eps.hat = y_ - t(t(A_)* beta.hat)
  } else {
    K = apply((A_ < Inf) * y_^2, 2, mean.na)
    beta.hat =  apply(A_ * y_, 2,  mean.na) / K; 
    eps.hat = A_ - as.matrix(y_) %*% t(beta.hat)}
  
  if(!robust.se){
    res = tibble::tibble(feature = colnames(A_),
                         beta.hat = beta.hat,
                         beta.se = sqrt(apply(eps.hat^2, 2, mean.na) / K / (N.A - 1 - dim(W_)[2])),
                         p.val = 2*pnorm(-abs(beta.hat / beta.se)),
                         q.val.BH = p.adjust(p.val, method = "BH"), 
                         n = N.A)
  } else{
    res = tibble::tibble(feature = colnames(A_),
                         beta.hat = beta.hat,
                         beta.se = sqrt(apply(A_^2 * eps.hat^2, 2, mean.na) / (N.A - 1 - dim(W_)[2]))/K , 
                         p.val = 2*pnorm(-abs(beta.hat / beta.se)),
                         q.val.BH = p.adjust(p.val, method = "BH"),
                         n = N.A)}
  
  res$beta.se[res$beta.se < .000001] <- .000001 # deals with bug in ash
  res = dplyr::bind_cols(res[,c(1,4,5,6)], 
                         ashr::ash(res$beta.hat, res$beta.se)$result[,c(1,2,7,8,9,10)])
  
  res %<>%dplyr::mutate(z.score = PosteriorMean/PosteriorSD)
  
  return(res)
}  

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

