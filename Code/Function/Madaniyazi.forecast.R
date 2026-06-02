#' Madaniyazi comparator: forecast mortality rates
#' Uses a fitted DLNM--Madaniyazi GLM to forecast deaths, then converts
#' predictions to annualized mortality rates.
#'
#' @param fit Fitted model object returned by \code{Madaniyazi.fit()}.
#' @param exposure Numeric vector of forecast-period exposures.
#' @param utci_crossbasis Cross-basis matrix for UTCI over the forecast horizon.
#' @param weekly_trend Numeric vector giving the forecast-period linear weekly
#'   trend.
#' @param week_of_year_basis Matrix or data.frame containing forecast-period
#'   cyclic seasonal basis terms.
#'
#' @return Numeric vector of annualized mortality-rate forecasts,
#'   computed as predicted deaths divided by exposure and multiplied by 52.
Madaniyazi.forecast <- function(fit, exposure, utci_crossbasis,
                                weekly_trend, week_of_year_basis) {
  pred_death <- predict(
    fit,
    newdata = list(exposure = exposure,
                   utci_crossbasis = utci_crossbasis,
                   weekly_trend = weekly_trend,
                   week_of_year_basis = week_of_year_basis),
    type = "response",
    se.fit = TRUE
  )$fit

  pred_death / exposure * 52
}
