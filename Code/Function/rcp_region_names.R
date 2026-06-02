# Region and age-group names used by the RCP projection workflow.
# Internal names preserve legacy data labels; output names use Athens/Lisbon/Rome.
rcp_internal_region_names <- c("Attiki", "Lisbon", "Roma")
rcp_region_names <- c("Athens", "Lisbon", "Rome")
rcp_display_names <- c(Attiki = "Athens", Lisbon = "Lisbon", Roma = "Rome")
rcp_output_region_names <- rcp_region_names
rcp_age_names <- c("Y20_64", "Y65_74", "Y75_84", "Y_GE85")

#' Rename list entries using a region-name map
#'
#' @param x Object with optional \code{names()}, usually a region-indexed list.
#' @param name_map Named character vector mapping current names to new names.
#'
#' @return \code{x} with matching names replaced by \code{name_map}.
rename_region_names <- function(x, name_map) {
  current_names <- names(x)
  if (!is.null(current_names)) {
    names(x) <- ifelse(current_names %in% names(name_map),
                       unname(name_map[current_names]),
                       current_names)
  }
  x
}

#' Convert display region names to legacy internal names
#'
#' @param x Object with optional region names.
#'
#' @return \code{x} with \code{Athens} renamed to \code{Attiki} and
#'   \code{Rome} renamed to \code{Roma}.
to_internal_region_names <- function(x) {
  rename_region_names(x, c(Athens = "Attiki", Rome = "Roma"))
}

#' Convert legacy internal names to display names
#'
#' @param x Object with optional region names.
#'
#' @return \code{x} with legacy names converted to \code{Athens},
#'   \code{Lisbon}, and \code{Rome}.
to_output_region_names <- function(x) {
  rename_region_names(x, rcp_display_names)
}

normalise_region_names <- to_output_region_names
