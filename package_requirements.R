pkgs <- c(
  # Data import / export
  "readxl",     # read Excel files
  "writexl",    # write Excel files

  # Data manipulation / reshaping / imputation helpers
  "dplyr",      # data manipulation (filter/mutate/summarise/joins)
  "tidyr",      # pivoting + tidying
  "reshape2",   # melt/cast reshaping (legacy but still used)
  "zoo",        # time series infrastructure + na.* + rolling ops
  "ISOweek",    # ISO week/date conversion utilities

  # Modelling
  "dlnm",       # distributed lag non-linear models
  "splines",    # spline bases (ns, bs)
  "forecast",   # forecasting + ARIMA/ETS helpers
  "astsa",      # applied time series analysis utilities
  "MTS",        # multivariate time series models/tools
  "mgcv",       # GAMs/GAMMs
  "demography", # mortality/demographic models
  "Metrics",    # Evaluation metric

  # Diagnostics / statistical tests
  "seastests",  # seasonality tests
  "tseries",    # time series tests (ADF, KPSS, etc.)
  "nortest",    # normality tests
  "hwwntest",   # heteroskedasticity / related tests

  # Visualization
  "ggplot2",    # plotting
  "patchwork",  # compose ggplots
  "gridExtra",  # arrange grobs/plots
  "ggpubr",     # publication helpers for ggplot
  "scales"      # axis/scale formatting and transforms
)

missing <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
invisible(lapply(pkgs, library, character.only = TRUE))
