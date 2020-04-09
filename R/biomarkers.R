require(tidyverse)
require(taigr)
require(magrittr)
require(readr)
require(here)

source(here("R/discrete_association.R"))
source(here("R/linear_association.R"))
source(here("R/random_forest.R"))

#' Gets the biomarkers for each column in a matrix of responses using:
#' random_forest
#' lin_associations
#' discrete_test
#'
#' @param Y n x p numerical matrix of responses to perturbations,
#'   rownames must be Arxspan IDs (missing values are allowed).
#' @param p_cutoff maximum p-value in output tables (to keep them smaller), defaults to 0.1
#' @param out_path if specified, string path to folder to write results to as a .csv
#'
#' @return a list with components
#' \describe{
#'   \item{rf_table}{A table with the following columns:
#'   feature: the column names of X.
#'   RF.imp.mean: an estimate of the importance of that feature for model accuracy.
#'   RF.imp.sd: the standard deviation of the importance estimate.
#'   RF.imp.stability: the proportion of models that used this feature.
#'   rank: the rank of the feature in terms of importance for the model.
#'   MSE: the mean-squared error of the model.
#'   MSE.se: the standard error of the MSE.
#'   R2: the \eqn{R^2} of the model.
#'   PearsonScore: the Pearson correlation of predicted and observed responses.
#'   pert: the column names of Y
#'   }
#'   \item{lin_table}{A table with the following columns:
#'   feat.A : column names of A
#'   p.val : p-values computed for the null-hypothesis beta = 0
#'   q.val : q-values computed using Benjamini-Hotchberg method
#'   n : number of non-missing values for each column of A
#'   betahat : estimated linear coefficient controlled for W
#'   sebetahat ; standard error for the estimate of beta
#'   lfdr and qvalue: local and global fdr values computed using adaptive shrinkage as
#'   described in doi:10.1093/biostatistics/kxw041
#'   PosteriorMean and PosteriorSD : moderated effect size estimates and corresponding standard deviations,
#'   again computed with adaptive shrinkage.
#'   pert: the column names of Y
#'   }
#'   \item{disc_table}{A table with the following columns:
#'   feature: the column names of X.
#'   effect_size: an estimate of the difference between groups.
#'   t.stat: the t-statistic associated with the test on group differences.
#'   p.value: the p-value of the t-statistic.
#'   q.value: the multiple hypothesis corrected p-value.
#'   pert: the column names of Y
#'   }
#' }
#'
#' @export
#'
get_biomarkers <- function(Y, p_cutoff=0.1, out_path=NULL) {
  require(magrittr)
  
  # lists tracking which feature sets are associated with which functions
  rf_data <- c("all_features", "ccle_features")
  discrete_data <- c("lineage", "mutation")
  linear_data <- c("copy_number", "dependency_shRNA", "dependency_XPR",
                   "expression", "miRNA", "repurposing", "RPPA", "total_proteome")

  # output tables
  linear_table <- tibble::tibble()
  discrete_table <- tibble::tibble()
  random_forest_table <- tibble::tibble()

  # linear associations
  for(feat in linear_data) {

    # load feature set
    X <- taigr::load.from.taiga(data.name='biomarker-features-699b',
                                data.version=5, data.file=feat, quiet=T)

    # for each perturbation get results
    for(pert in colnames(Y)) {
      # select specific column and filter to finite values
      y <- Y[,pert]; names(y) <- rownames(Y)
      y <- y[is.finite(y)]

      # get overlapping data
      overlap <- dplyr::intersect(rownames(X), names(y))
      # calculate correlations
      res.lin <- lin_associations(X[overlap,], y[overlap])

      # if specified p-value cutoff, filter values above cutoff
      if(!is.null(p_cutoff)) {
        res.lin %<>% dplyr::filter(p.val <= p_cutoff)
      }

      # append to output tables
      linear_table %<>% dplyr::bind_rows(res.lin %>%
          dplyr::mutate(pert = pert, feature_type = feat))
    }
  }
  # repeat for discrete t-test
  for(feat in discrete_data) {

    X <- taigr::load.from.taiga(data.name='biomarker-features-699b',
                                data.version=5, data.file=feat, quiet=T)

    # only do test for colummns with more that 1 in group
    X <- X[,apply(X, 2, function(x) sum(x) > 1)]

    for(pert in colnames(Y)) {
      y <- Y[,pert]; names(y) <- rownames(Y)
      y <- y[is.finite(y)]

      overlap <- dplyr::intersect(rownames(X), names(y))
      res.disc <- discrete_test(X[overlap,], y[overlap])

      if(!is.null(p_cutoff)) {
        res.disc %<>% dplyr::filter(p.value <= p_cutoff)
      }

      discrete_table %<>% dplyr::bind_rows(res.disc %>%
         dplyr::mutate(pert = pert, feature_type = feat))
    }
  }
  # repeat for random forest
  for(feat in rf_data) {
    X <- taigr::load.from.taiga(data.name='biomarker-features-699b',
                                data.version=5, data.file=feat, quiet=T)
    for(pert in colnames(Y)) {
      y <- Y[,pert]; names(y) <- rownames(Y)
      y <- y[is.finite(y)]

      overlap <- dplyr::intersect(rownames(X), names(y))
      res.rf <- random_forest(X[overlap,], y[overlap])
      random_forest_table %<>% dplyr::bind_rows(res.rf$model_table %>%
           dplyr::mutate(pert = pert, feature_set = feat))
    }
  }

  # write .csv if path given
  if(!is.null(out_path)) {
    readr::write_csv(random_forest_table, paste(out_path, "rf_table.csv", sep = "/"))
    readr::write_csv(linear_table, paste(out_path, "lin_association_table.csv", sep = "/"))
    readr::write_csv(discrete_table, paste(out_path, "discrete_table.csv", sep = "/"))
  }

  # return
  return(list("rf_table" = random_forest_table,
              "lin_table" = linear_table,
              "disc_table" = discrete_table))
}
