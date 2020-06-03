require(tidyverse)
require(taigr)
require(magrittr)
require(readr)
require(cdsrmodels)

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
#'   \item{rf_table}{A table with the following columns: \cr
#'   feature: the column names of X. \cr
#'   RF.imp.mean: an estimate of the importance of that feature for model accuracy. \cr
#'   RF.imp.sd: the standard deviation of the importance estimate. \cr
#'   RF.imp.stability: the proportion of models that used this feature. \cr
#'   rank: the rank of the feature in terms of importance for the model. \cr
#'   MSE: the mean-squared error of the model. \cr
#'   MSE.se: the standard error of the MSE. \cr
#'   R2: the \eqn{R^2} of the model. \cr
#'   PearsonScore: the Pearson correlation of predicted and observed responses. \cr
#'   pert: the column names of Y \cr
#'   }
#'   \item{lin_table}{A table with the following columns (values calculated with adaptive shrinkage): \cr
#'   betahat: estimated linear coefficient controlled for W. \cr
#'   sebetahat: standard error for the estimate of beta. \cr
#'   NegativeProb: posterior probabilities that beta is negative. \cr
#'   PositiveProb: posterior probabilities that beta is positive \cr
#'   lfsr: local and global fsr values. \cr
#'   svalue: s-values.\cr
#'   lfdr: local and global fdr values. \cr
#'   qvalue: q-values for the effect size estimates. \cr
#'   PosteriorMean: moderated effect size estimates. \cr
#'   PosteriorSD: standard deviations of moderated effect size estimates. \cr
#'   dep.var: dependent variables (column names of Y). \cr
#'   ind.var: independent variables (column names of X). \cr
#'   }
#'   \item{disc_table}{A table with the following columns:
#'   feature: the column names of X. \cr
#'   effect_size: an estimate of the difference between groups. \cr
#'   t.stat: the t-statistic associated with the test on group differences. \cr
#'   p.value: the p-value of the t-statistic. \cr
#'   q.value: the multiple hypothesis corrected p-value. \cr
#'   pert: the column names of Y \cr
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
      res.lin <- cdsrmodels::lin_associations(X[overlap,], y[overlap])$res.table
      res.lin <- tibble::as_tibble(res.lin)

      # append to output tables
      linear_table %<>% dplyr::bind_rows(res.lin %>%
          dplyr::mutate(feature_type = feat))
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
      res.disc <- cdsrmodels::discrete_test(X[overlap,], y[overlap])
      
      # if cutoff specified filter p values above cutoff
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
      res.rf <- cdsrmodels::random_forest(X[overlap,], y[overlap])
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

#' Generates the multi response profile biomarker report
#'
#' @param out_path string path to folder. If Y is NULL this folder should contain biomarker
#'  output files rf_table.csv etc... If Y is given the biomarker output files will be written to this
#'  folder
#' @param title string title for the report
#' @param Y optional n x p numerical matrix of responses to perturbations,
#'   rownames must be Arxspan IDs (missing values are allowed).
#' @param meta_data optional a dataframe containing meta data to include in the report
#'
#' @export
#'
generate_multi_profile_biomarker_report <- function(out_path, title, Y = NULL, meta_data = NULL) {
  if(!is.null(Y)) {
    get_biomarkers(Y, out_path = out_path)
    Y %>% as_tibble(rownames = "arxspan_id") %>% write_csv(paste(out_path, "data.csv", sep = "/"))
  }
  if(!is.null(meta_data)) {
    meta_data %>% write_csv(paste(out_path, "meta_data.csv", sep = "/"))
  }
  rmarkdown::render(system.file("reports", "multi_profile_biomarker_report.Rmd",
                                package = "cdsrbiomarker"),
                    params = list(in_path = out_path, title = title),
                    output_dir = out_path,output_file = title)
}

#' Generates the single response profile biomarker report
#'
#' @param out_path string path to folder. If Y is NULL this folder should contain biomarker
#'  output files rf_table.csv etc... If Y is given the biomarker output files will be written to this
#'  folder
#' @param title string title for the report
#' @param Y optional n x p numerical matrix of responses to perturbations,
#'   rownames must be Arxspan IDs (missing values are allowed).
#' @param meta_data optional a dataframe containing meta data to include in the report
#'
#' @export
#'
generate_single_profile_biomarker_report <- function(out_path, title, Y = NULL, meta_data = NULL) {
  if(!is.null(Y)) {
    get_biomarkers(Y, out_path = out_path)
    Y %>% as_tibble(rownames = "arxspan_id") %>% write_csv(paste(out_path, "data.csv", sep = "/"))
  }
  if(!is.null(meta_data)) {
    meta_data %>% write_csv(paste(out_path, "meta_data.csv", sep = "/"))
  }
  rmarkdown::render(system.file("reports", "multi_single_biomarker_report.Rmd",
                                package = "cdsrbiomarker"),
                    params = list(in_path = out_path, title = title),
                    output_dir = out_path,output_file = title)
}
