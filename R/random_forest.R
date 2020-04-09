require(tidyverse)
require(magrittr)
require(ranger)
require(data.table)

#' Fits a random forest from a matrix or features X and a vector y.
#'
#' @param X n x m numerical matrix of features (missing values are allowed).
#' @param y length n vector of numerical values (missing values are allowed).
#' @param k integer number of cross validation cycles to perform.
#' @param n number of features to be considered in the model (after correlation filter).
#'
#' @return a list with components
#' \describe{
#'   \item{model_table}{A table with the following columns:
#'   feature: the column names of X.
#'   RF.imp.mean: an estimate of the importance of that feature for model accuracy.
#'   RF.imp.sd: the standard deviation of the importance estimate.
#'   RF.imp.stability: the proportion of models that used this feature.
#'   rank: the rank of the feature in terms of importance for the model.
#'   MSE: the mean-squared error of the model.
#'   MSE.se: the standard error of the MSE.
#'   R2: the \eqn{R^2} of the model.
#'   PearsonScore: the Pearson correlation of predicted and observed responses.
#'   }
#'   \item{predictions}{A vector of the responses predicted by the random forest.}
#' }
#'
#' @export
#'
random_forest <- function(X, y, k = 10, n = 250){
  y <- y[is.finite(y)]  # only finite values
  X <- X[, apply(X, 2, function(x) !any(is.na(x)))]
  cl <- sample(dplyr::intersect(rownames(X), names(y)))  # overlapping rows
  X <- scale(X[cl,])  # scales the columns of X and selects overlapping lines
  y = y[cl]  # selects overlapping lines from y
  colnames(X) %<>% make.names()  # ensure unique column names

  N = floor(length(cl)/k)  # size of each test set (dataset size / num cv cycles)
  yhat_rf <- y  # vector to store predicted values of y

  SS = tibble()  # empty table to store output

  # run cross validation k times
  for (kx in 1:k) {
    test <- cl[((kx - 1) * N + 1):(kx * N)]  # select test set (N points)
    train <- dplyr::setdiff(cl, test)  # everything else is training

    # select training and test data from X, assumes no NAs in X
    X_train <- X[train,]
    X_test <- X[test,]

    X_train <- X_train[,apply(X_train,2,var) > 0]  # remove columns with variance 0

    # select top n correlated features in X (this filters to "relevant" features)
    # allows for faster model fitting
    X_train <- X_train[,rank(-abs(cor(X_train,y[train]))) <= n]

    # fit a random forest model using the ranger package
    # uses impurity as varaible importance metric
    rf <- ranger::ranger(y ~ .,
                         data = as.data.frame(cbind(X_train, y = y[train])),
                         importance = "impurity")

    # add predictions for test set to prediction vector
    yhat_rf[test] <- predict(rf, data = as.data.frame(X_test[, colnames(X_train)]))$predictions

    # extract variable importance metrics from RF model
    ss <- tibble::tibble(feature = names(rf$variable.importance),
                         RF.imp = rf$variable.importance / sum(rf$variable.importance),
                         RF.mse = mean((yhat_rf[test] - y[test])^2),
                         fold = kx)

    # add results from this step of cross-validation to overall model
    SS %<>%
      dplyr::bind_rows(ss)
  }

  # separate into RF table and one large table of model level statistics
  # RF table will have mean variable importances for each feature

  model_table <- tibble::tibble()  # model level statistics (empty for now)

  # importances table in matrix form (feature by CV run with values as importances)
  RF.importances <- SS %>%
    dplyr::distinct(feature, RF.imp, fold) %>%
    reshape2::acast(feature ~ fold, value.var = "RF.imp")

  # RF table with average importance values across validation runs
  # mean and sd are of the importance values
  # stability measures how many models used a featuer (divided by number of CV folds)
  RF.table <- tibble::tibble(feature = rownames(RF.importances),
                             RF.imp.mean = RF.importances %>%
                               apply(1, function(x) mean(x, na.rm = T)),
                             RF.imp.sd = RF.importances %>%
                               apply(1, function(x) sd(x, na.rm = T)),
                             RF.imp.stability = RF.importances %>%
                               apply(1, function(x) mean(!is.na(x)))) %>%
    # only keep features used in > half the models
    dplyr::filter((RF.imp.stability > 0.5), feature != "(Intercept)") %>%
    dplyr::arrange(desc(RF.imp.mean)) %>%
    dplyr::mutate(rank = 1:n())

  # table of RF model level statistics
  # MSE = mean-squared error of predictions, MSE.se = standard deviation of MSE
  # R2 and Pearson are measures of model accuracy
  RF.table %<>%
    dplyr::mutate(MSE = mean(dplyr::distinct(SS, RF.mse, fold)$RF.mse,
                             na.rm =T),
                  MSE.se = sd(dplyr::distinct(SS, RF.mse, fold)$RF.mse,
                              na.rm = T),
                  R2 = 1 - MSE/var(y[1:(k*N)], na.rm = T),
                  PearsonScore = cor(y[1:(k*N)], yhat_rf[1:(k*N)],
                                     use = "pairwise.complete.obs"))

  # return importance table, model level table, and predictions
  return(list("model_table" = RF.table,
              "predictions" = yhat_rf))
}
