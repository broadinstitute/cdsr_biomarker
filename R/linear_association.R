library(tidyverse)
library(magrittr)
library(useful)
library(taigr)

#' Calculates the linear association between a matrix or features A and a vector y. 
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
#' @return a data frame with the following columns
#'feat.A : column names of A
#'p.val : p-values computed for the null-hypothesis beta = 0
#'q.val : q-values computed using Benjamini-Hotchberg method
#'n : number of non-missing values for each column of A
#'betahat : estimated linear coefficient controlled for W
#'sebetahat ; standard error for the estimate of beta
#'lfdr and qvalue: local and global fdr values computed using adaptive shrinkage as 
#'described in doi:10.1093/biostatistics/kxw041 
#'PosteriorMean and PosteriorSD : moderated effect size estimates and corresponding standard deviations, 
#'again computed with adaptive shrinkage.
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
    A_ = A[, apply(A,2,function(x) var(x,na.rm = T)/mean(x,na.rm = T)) > var.th]
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
  
  res = dplyr::bind_cols(res[,c(1,4,5,6)], 
                         ashr::ash(res$beta.hat, res$beta.se)$result[,c(1,2,7,8,9,10)])
  return(res)
}  
