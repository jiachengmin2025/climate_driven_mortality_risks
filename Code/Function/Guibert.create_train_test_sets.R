#' Guibert comparator: create expanding-window train/test inputs
#' Splits region-specific mortality data and rebuilds DLNM cross-basis matrices
#' for each training/testing fold.
#'
#' @param dat_list Named list of mortality-rate data.frames/matrices by region.
#' @param cb_list Named list of UTCI lag matrices aligned with \code{dat_list}.
#' @param train_size Integer length of the training window.
#' @param test_size Integer end index of the testing window.
#' @param lag Integer maximum lag passed to \code{dlnm::crossbasis()}.
#' @param var.df Degrees of freedom for the exposure dimension.
#' @param lag.df Degrees of freedom for the lag dimension.
#' @param degree Polynomial degree passed to \code{dlnm::equalknots()}.
#'
#' @return A nested list with \code{train_set} and \code{test_set}. Each element
#'   contains \code{dat_list} and \code{cb_list} entries for each region.
Guibert.create_train_test_sets <- function(dat_list, cb_list, train_size, test_size, lag, var.df, lag.df, degree) {
  
  # Initialize lists for train and test sets
  train_test_data <- list(train_set = list(), test_set = list())
  
  # Iterate over each city (excluding Common)
  for (city_name in names(dat_list)) {
    # Train and test split for dat_list
    city_train <- dat_list[[city_name]][1:train_size, ]
    cb_train <- cb_list[[city_name]][1:train_size, ]
    
    city_test <- dat_list[[city_name]][(train_size + 1):test_size, ]
    cb_test <- cb_list[[city_name]][(train_size + 1):test_size, ]
    
    # Smoothing cross-basis matrix (train)
    varknots_train <- equalknots(cb_train, fun = "ns", df = var.df, degree = degree)
    lagknots_train <- logknots(lag, lag.df)
    cb_train <- crossbasis(cb_train, lag = lag,
                           argvar = list(fun = "ns", knots = varknots_train, df = var.df),
                           arglag = list(knots = lagknots_train, df = lag.df))
    
    # Smoothing cross-basis matrix (test) *** use training set setting! ***
    varknots_test <- equalknots(cb_train, fun = "ns", df = var.df, degree = degree)
    lagknots_test <- logknots(lag, lag.df)
    cb_test <- crossbasis(cb_test, lag = lag,
                          argvar = list(fun = "ns", knots = varknots_test, df = var.df),
                          arglag = list(knots = lagknots_test, df = lag.df))
    
    # Store train and test sets in the corresponding lists
    train_test_data$train_set$dat_list[[city_name]] <- city_train
    train_test_data$train_set$cb_list[[city_name]] <- cb_train
    
    train_test_data$test_set$dat_list[[city_name]] <- city_test
    train_test_data$test_set$cb_list[[city_name]] <- cb_test
  }
  
  return(train_test_data)
}
