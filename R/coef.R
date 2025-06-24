


#' Extract fitted model coefficients
#'
#' @description
#' `coef` extracts fitted model coefficients from
#'    fitted model objects. `coefficients.s4t_cjs_ml` is an alias for it
#'
#' @param object a fitted model object from `s4t_cjs_ml()`
#' @param ... Other arguments. Not used (needed for generic consistency)
#' @return A named vector of coefficients
#'
#' @importFrom stats coef
#' @export
#'
#' @method coef s4t_cjs_ml
#'
#' @examples
#' # not run
#'
coef.s4t_cjs_ml <- function(object, ...) {
  parnames <- object$estimated_parameters[,"parameter"]
  ests <- as.vector(object$estimated_parameters[,"estimate"])

  names(ests) <- parnames
  return(ests)
}

#' @rdname coef.s4t_cjs_ml
#' @method coefficients s4t_cjs_ml
#'
#' @importFrom stats coefficients
#'
#' @export
#'
coefficients.s4t_cjs_ml <- coef.s4t_cjs_ml


#' Extract fitted model coefficients
#'
#' @description
#' `coef` extracts fitted model coefficients from
#'    fitted model objects. `coefficients.s4t_ageclass_model` is an alias for it
#'
#' @param object a fitted model object from `s4t_ageclass_model()`
#' @param ... Other arguments. Not used (needed for generic consistency)
#' @return A named vector of coefficients
#'
#' @export
#'
#' @method coef s4t_ageclass_model
#'
#' @examples
#' # not run
#'
coef.s4t_ageclass_model <- function(object, ...) {
  parnames <- object$estimated_parameters[,"parameter"]
  ests <- as.vector(object$estimated_parameters[,"estimate"])

  names(ests) <- parnames
  return(ests)
}

#' @rdname coef.s4t_ageclass_model
#' @method coefficients s4t_ageclass_model
#'
#' @export
#'
coefficients.s4t_ageclass_model <- coef.s4t_ageclass_model

