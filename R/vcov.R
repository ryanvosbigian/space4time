
#' returns variance-covariance matrix
#'
#' @param object A `s4t_cjs` object
#' @param ... Other arguments. Not used (for generic consistency)
#' @return The variance-covariance matrix for the CJS portion of the model
#' @export
vcov.s4t_cjs <-  function(object, ...) {
  solve(object$res$hessian)
}
