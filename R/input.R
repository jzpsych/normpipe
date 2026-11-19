# ==============================================================================
# input.R -- checks of the samples handed to run_norming()
#
# A sample is list(meta, items): one row per respondent in `meta` (with `id`)
# and the item responses in `items`. The rows of `items` must belong to the
# rows of `meta`. When `items` has an `id` column, the rows are matched by it;
# otherwise both must have the same number of rows, and the order is taken as
# given.
# ==============================================================================

#' Check a sample `list(meta, items)` and bring `items` into the row order of
#' `meta`
#'
#' @param x `list(meta, items)` or `NULL`.
#' @param what Name of the argument, for messages.
#' @param unique_ids Must the ids of `meta` be unique (reference and
#'   calibration samples); donor samples may repeat ids across samples.
#' @return `x`, with `items` reordered by `items$id` (which is then dropped)
#'   when it has an `id` column; otherwise unchanged.
#' @noRd
align_sample <- function(x, what, unique_ids = TRUE) {
  if (is.null(x)) return(NULL)
  if (!is.list(x) || is.null(x$items) || !is.data.frame(x$items))
    stop("`", what, "` must be a list with a data frame `items` (and `meta`).", call. = FALSE)
  items <- x$items; meta <- x$meta
  has_id <- "id" %in% names(items)
  if (is.null(meta)) {
    if (has_id) x$items <- items[setdiff(names(items), "id")]
    return(x)
  }
  if (!is.data.frame(meta)) stop("`", what, "$meta` must be a data frame.", call. = FALSE)
  if (!"id" %in% names(meta)) stop("`", what, "$meta` needs a column `id`.", call. = FALSE)
  ids <- as.character(meta$id)
  if (anyNA(ids)) stop("`", what, "$meta$id` has missing values.", call. = FALSE)
  if (unique_ids && anyDuplicated(ids))
    stop("`", what, "$meta$id` is not unique (e.g. '", ids[anyDuplicated(ids)], "').", call. = FALSE)
  if (has_id) {
    iid <- as.character(items$id)
    if (anyDuplicated(ids))
      stop("`", what, "`: `items$id` can only be matched when `meta$id` is unique.", call. = FALSE)
    if (anyNA(iid) || anyDuplicated(iid) || length(iid) != length(ids) || !setequal(iid, ids))
      stop("`", what, "`: the ids in `items$id` are not the ids of `meta$id` (",
           sum(!ids %in% iid), " id(s) of meta without item data, ", sum(!iid %in% ids),
           " id(s) of the item data not in meta).", call. = FALSE)
    rows <- match(ids, iid)
    x$items <- items[rows, setdiff(names(items), "id"), drop = FALSE]
    return(x)
  }
  if (nrow(items) != nrow(meta))
    stop("`", what, "`: `items` has ", nrow(items), " rows, `meta` has ", nrow(meta), ". ",
         "The rows of `items` must be the respondents of `meta` in the same order; ",
         "an `id` column in `items` lets normpipe match them.", call. = FALSE)
  x
}
