#' Fit the selected SARIMA model for a time-varying factor
#'
#' @param x Numeric time series of \code{K(t)} or \code{kappa(t, i)} values.
#'
#' @return Fitted \code{forecast::auto.arima} model.
fit_sarima_factor <- function(x) {
  forecast::auto.arima(
    ts(as.numeric(x), frequency = 52),
    seasonal = TRUE,
    max.d = 0,
    allowmean = FALSE
  )
}

#' Format SARIMA order for LaTeX
#'
#' @param fit Fitted \code{forecast::auto.arima} model.
#'
#' @return Character string containing a LaTeX-formatted SARIMA order.
format_sarima_order <- function(fit) {
  arma <- fit$arma
  paste0(
    "ARIMA$(", arma[1], ",", arma[6], ",", arma[2], ")",
    "\\times(", arma[3], ",", arma[7], ",", arma[4], ")_{", arma[5], "}$"
  )
}

#' Convert a SARIMA coefficient p-value to significance stars
#'
#' @param p_value Numeric p-value.
#'
#' @return Character string containing LaTeX superscript significance stars.
sarima_significance_stars <- function(p_value) {
  ifelse(p_value < 0.01, "^{***}",
         ifelse(p_value < 0.05, "^{**}",
                ifelse(p_value < 0.1, "^{*}", "")))
}

#' Format one SARIMA coefficient for LaTeX
#'
#' @param fit Fitted \code{forecast::auto.arima} model.
#' @param term Coefficient name, e.g. \code{"ar1"} or \code{"sma1"}.
#'
#' @return Character string containing the formatted coefficient estimate,
#'   significance stars, or \code{"$-$"} when the term is absent.
format_sarima_coef <- function(fit, term) {
  estimates <- coef(fit)
  if (!term %in% names(estimates)) return("$-$")

  var_coef <- fit$var.coef
  if (is.null(var_coef) || !term %in% rownames(var_coef)) {
    return(paste0("$", formatC(estimates[[term]], format = "f", digits = 3), "$"))
  }

  std_error <- sqrt(diag(var_coef))[term]
  p_value <- 2 * pnorm(abs(estimates[[term]] / std_error), lower.tail = FALSE)
  paste0("$", formatC(estimates[[term]], format = "f", digits = 3),
         sarima_significance_stars(p_value), "$")
}

#' Build the LaTeX table for selected SARIMA orders
#'
#' @param order_table Data.frame with \code{time_series}, \code{DLNM_LC}, and
#'   \code{DLNM_LL} columns.
#'
#' @return Character string containing a complete LaTeX table.
make_sarima_order_latex <- function(order_table) {
  lines <- c(
    "\\begin{table}[H]",
    "    \\centering",
    "    \\begin{tabular}{ccc}",
    "    \\toprule",
    "    \\small\\textbf{Time series model} & \\small\\textbf{DLNM--LC} & \\small\\textbf{DLNM--LL} \\\\",
    "    \\midrule"
  )

  for (i in seq_len(nrow(order_table))) {
    lines <- c(lines, paste0(
      "    ", order_table$time_series[[i]], " & ",
      order_table$DLNM_LC[[i]], " & ",
      order_table$DLNM_LL[[i]], "\\\\"
    ))
  }

  paste(c(
    lines,
    "    \\bottomrule",
    "    \\end{tabular}",
    "    \\caption{The selected optimal SARIMA model for time-varying factors.}",
    "    \\label{tab:optimal_sarima}",
    "\\end{table}"
  ), collapse = "\n")
}

#' Build the LaTeX table for SARIMA coefficient estimates
#'
#' @param sarima_models Nested list of fitted SARIMA models for DLNM--LC and
#'   DLNM--LL time-varying factors.
#'
#' @return Character string containing a complete LaTeX table.
make_sarima_coef_latex <- function(sarima_models) {
  coef_terms <- c("ar1", "ar2", "ma1", "ma2", "sar1", "sma1")
  header <- paste0(
    " & $\\hat{\\phi}_1$ & $\\hat{\\phi}_2$ & ",
    "$\\hat{\\theta}_1$ & $\\hat{\\theta}_2$ & ",
    "$\\hat{\\Phi}_1$ & $\\hat{\\Theta}_1$\\\\"
  )
#' SARIMA table internal: build one coefficient row
#'
#' @param label LaTeX row label.
#' @param fit Fitted SARIMA model.
#'
#' @return Character string containing one LaTeX table row.
  make_row <- function(label, fit) {
    paste0(
      "    ", label, " & ",
      paste(vapply(coef_terms, function(term) format_sarima_coef(fit, term),
                   character(1)), collapse = " & "),
      "\\\\"
    )
  }

  paste(c(
    "\\begin{table}[H]",
    "  \\centering",
    "  \\begin{tabular}{ccccccc}",
    "    \\toprule",
    paste0("    \\small\\textbf{DLNM--LC}", header),
    "    \\midrule",
    make_row("$\\kappa(t, \\text{Athens})$", sarima_models$LC$Athens),
    make_row("$\\kappa(t, \\text{Lisbon})$", sarima_models$LC$Lisbon),
    make_row("$\\kappa(t, \\text{Rome})$", sarima_models$LC$Rome),
    "    \\midrule",
    paste0("    \\small\\textbf{DLNM--LL}", header),
    "    \\midrule",
    make_row("$K(t)$", sarima_models$LL$Common),
    make_row("$\\kappa(t, \\text{Athens})$", sarima_models$LL$Athens),
    make_row("$\\kappa(t, \\text{Lisbon})$", sarima_models$LL$Lisbon),
    make_row("$\\kappa(t, \\text{Rome})$", sarima_models$LL$Rome),
    "    \\bottomrule",
    "  \\end{tabular}",
    "  \\vspace{6pt}\\\\",
    "  \\footnotesize{Level of significance: $0.01$:***,\\quad $0.05$:**,\\quad $0.1$:*}",
    "  \\caption{The coefficients of selected optimal SARIMA models on time-varying factors.}",
    "  \\label{tab:optimal_sarima_coef}",
    "\\end{table}"
  ), collapse = "\n")
}

#' Build SARIMA order and coefficient tables from fitted DLNM models
#'
#' @param b Fitted multi-region \code{DLNM_LL()} object.
#' @param b1 Fitted Athens \code{DLNM_LC()} object.
#' @param b2 Fitted Lisbon \code{DLNM_LC()} object.
#' @param b3 Fitted Rome \code{DLNM_LC()} object.
#' @param output_dir Directory where SARIMA tables and model objects are saved.
#'
#' @return A list containing fitted SARIMA models, the order table, and LaTeX
#'   strings for both SARIMA tables.
build_sarima_tables_from_fits <- function(b, b1, b2, b3, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  sarima_models <- list(
    LC = list(
      Athens = fit_sarima_factor(b1$final_LC$k_t),
      Lisbon = fit_sarima_factor(b2$final_LC$k_t),
      Rome = fit_sarima_factor(b3$final_LC$k_t)
    ),
    LL = list(
      Common = fit_sarima_factor(b$final_LL$Kt),
      Athens = fit_sarima_factor(b$final_LL$kt$Athens),
      Lisbon = fit_sarima_factor(b$final_LL$kt$Lisbon),
      Rome = fit_sarima_factor(b$final_LL$kt$Rome)
    )
  )

  sarima_order_table <- data.frame(
    time_series = c(
      "$K(t)$",
      "$\\kappa(t, \\text{Athens})$",
      "$\\kappa(t, \\text{Lisbon})$",
      "$\\kappa(t, \\text{Rome})$"
    ),
    DLNM_LC = c(
      "--",
      format_sarima_order(sarima_models$LC$Athens),
      format_sarima_order(sarima_models$LC$Lisbon),
      format_sarima_order(sarima_models$LC$Rome)
    ),
    DLNM_LL = c(
      format_sarima_order(sarima_models$LL$Common),
      format_sarima_order(sarima_models$LL$Athens),
      format_sarima_order(sarima_models$LL$Lisbon),
      format_sarima_order(sarima_models$LL$Rome)
    ),
    stringsAsFactors = FALSE
  )

  sarima_order_latex <- make_sarima_order_latex(sarima_order_table)
  sarima_coef_latex <- make_sarima_coef_latex(sarima_models)

  saveRDS(sarima_models, file.path(output_dir, "optimal_sarima_models.rds"))
  saveRDS(sarima_order_table, file.path(output_dir, "optimal_sarima_orders.rds"))
  writeLines(sarima_order_latex, file.path(output_dir, "optimal_sarima_table.tex"))
  writeLines(sarima_coef_latex, file.path(output_dir, "optimal_sarima_coef_table.tex"))

  list(models = sarima_models, order_table = sarima_order_table,
       order_latex = sarima_order_latex, coef_latex = sarima_coef_latex)
}
