#' Madaniyazi comparator: fit one age-specific regional model
#' Fits the DLNM--Madaniyazi comparator using death counts, an exposure offset,
#' a UTCI cross-basis, a linear weekly trend, and cyclic week-of-year seasonality.
#'
#' @param death Numeric vector of weekly death counts for one region and age group.
#' @param exposure Numeric vector of weekly population exposures aligned with
#'   \code{death}.
#' @param utci_crossbasis Cross-basis matrix for UTCI over the training window.
#' @param weekly_trend Numeric vector giving the linear weekly time trend.
#' @param week_of_year_basis Matrix or data.frame containing cyclic seasonal
#'   basis terms for week of year.
#'
#' @return A fitted quasipoisson \code{glm} object with log link and
#'   \code{offset(log(exposure))}.
Madaniyazi.fit <- function(death, exposure, utci_crossbasis,
                           weekly_trend, week_of_year_basis) {
  madaniyazi_formula <- death ~ offset(log(exposure)) +
    utci_crossbasis +
    weekly_trend +
    week_of_year_basis

  glm(madaniyazi_formula,
      family = quasipoisson(link = "log"),
      data = list(death = death,
                  exposure = exposure,
                  utci_crossbasis = utci_crossbasis,
                  weekly_trend = weekly_trend,
                  week_of_year_basis = week_of_year_basis),
      na.action = na.exclude)
}
