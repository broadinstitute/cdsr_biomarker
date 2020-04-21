require(tidyverse)
require(magrittr)
require(ranger)

#' Fits a random forest from a matrix or features X and a vector y.
#'
#' @param X n x m numerical matrix of features (missing values are allowed).
#' @param y length n vector of numerical values (missing values are allowed).
#' @param k integer number of cross validation cycles to perform.
#' @param n number of features to be considered in the model (after correlation filter).
#'
#' @return a list with components
#' \describe{
#'   \item{model_table}{A table with the following columns: \cr
#'   feature: the column names of X. \cr
#'   RF.imp.mean: an estimate of the importance of that feature for model accuracy. \cr 
#'   RF.imp.sd: the standard deviation of the importance estimate. \cr 
#'   RF.imp.stability: the proportion of models that used this feature. \cr
#'   rank: the rank of the feature in terms of importance for the model. \cr 
#'   MSE: the mean-squared error of the model. \cr 
#'   MSE.se: the standard error of the MSE. \cr
#'   R2: the \eqn{R^2} of the model. \cr
#'   PearsonScore: the Pearson correlation of predicted and observed responses. \cr 
#'   }
#'   \item{predictions}{A vector of the responses predicted by the random forest.}
#' }
#'
#' @export
#'
random_forest <- function(X, y, k = 10, n = 500){
  y <- y[is.finite(y)]  # only finite values
  X <- X[, apply(X, 2, function(x) !any(is.na(x)))]
  
  # make sure data aligns
  cl <- sample(dplyr::intersect(rownames(X), names(y)))  # overlapping rows
  X <- scale(X[cl,])  # scales the columns of X and selects overlapping lines
  y <- y[cl]  # selects overlapping lines from y
  
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
                         fold = kx)

    # add results from this step of cross-validation to overall model
    SS %<>%
      dplyr::bind_rows(ss)
  }

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

  # model level statistics
  # MSE = mean-squared error of predictions, MSE.se = standard deviation of MSE
  # R2 and Pearson are measures of model accuracy
  mse <- mean((yhat_rf - y)^2)
  mse.se <- sqrt(var((yhat_rf - y)^2))/sqrt(length(y))
  r2 <- 1 - (mse/var(y))
  ps <- cor(y, yhat_rf, use = "pairwise.complete.obs")
  RF.table %<>%
    dplyr::mutate(MSE = mse,
                  MSE.se = mse.se,
                  R2 = r2,
                  PearsonScore = ps)

  # return importance table, model level table, and predictions
  return(list("model_table" = RF.table,
              "predictions" = yhat_rf))
}
