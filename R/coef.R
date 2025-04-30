


#' Extract fitted model coefficients
#'
#' @description
#' `coef` extracts fitted model coefficients from
#'    fitted model objects. `coefficients.s4t_cjs` is an alias for it
#'
#' @param object a fitted model object from `s4t_cjs_ml()`
#' @param ... Other arguments. Not used (needed for generic consistency)
#' @return A named vector of coefficients
#' @export
#' @method coef s4t_cjs
#' @order 1
#'
#' @examples
#' # not run
#'
coef.s4t_cjs <- function(object, ...) {
  parnames <- object$estimated_parameters[,"parameter"]
  ests <- as.vector(object$estimated_parameters[,"estimate"])

  names(ests) <- parnames
  return(ests)
}

#' @rdname coef.s4t_cjs
#' @method coefficients s4t_cjs
#' @order 2
#' @export
coefficients.s4t_cjs <- coef.s4t_cjs
