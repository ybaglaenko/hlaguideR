#' Restrict guides to an editing window
#'
#' Keeps rows whose `position` column falls inside the inclusive window. Rows
#' with a missing position are dropped only when `keep_na` is `FALSE`; the
#' default is to keep them, so that guides lacking a scored position are never
#' removed silently.
#'
#' @param guides A data frame carrying a position column.
#' @param min,max Inclusive window bounds, counted from the protospacer's 5'
#'   end.
#' @param position_col Name of the column holding the position.
#' @param keep_na Whether to retain rows with a missing position.
#' @return `guides`, filtered.
#' @examples
#' g <- splice_guides(example_record())
#' nrow(g)
#' nrow(filter_editing_window(g, 3, 10))
#' @export
filter_editing_window <- function(guides, min = 3L, max = 10L,
                                  position_col = "position", keep_na = TRUE) {
  if (!is.data.frame(guides)) stop("`guides` must be a data frame.", call. = FALSE)
  if (!position_col %in% names(guides))
    stop("Column '", position_col, "' not found in `guides`.", call. = FALSE)
  p <- guides[[position_col]]
  inside <- !is.na(p) & p >= min & p <= max
  out <- guides[inside | (keep_na & is.na(p)), , drop = FALSE]
  rownames(out) <- NULL
  out
}
