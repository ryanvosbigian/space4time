
#' summary of s4t_cjs object
#'
#' @param object A s4t_cjs object
#' @param ... Not used
#' @returns estimated parameters table and AIC
#' @export
summary.s4t_cjs <- function(object, ...) {
  print(object$estimated_parameters,digits = 4)
  print(data.frame(AIC = round(object$AIC,2)))
  return(list(estimated_parameters = object$estimated_parameters),
         AIC = data.frame(AIC = round(object$AIC,2)))
}


#' summary of s4t_cjs_rstan object
#'
#' @description
#' Calls `summary,stanfit-method`.
#'
#'
#' @param object A s4t_cjs_rstan object.
#' @param pars a character vector of parameter names.
#' @param probs a numerical vector of quantiles of interest.
#' @param ... Additional arguments to pass to `summary,stanfit-method`.
#'
#' @returns Returns a named list with ...
#' @export
summary.s4t_cjs_rstan <- function(object,
                                  pars,
                                  probs = c(0.025, 0.25, 0.50, 0.75, 0.975),
                                  ...) {
  rstan::summary(object$res,pars = pars, probs = probs, ...)
}

#' summary of s4t_ageclass_cjs object
#'
#' @param object A s4t_ageclass_cjs object
#' @param ... Not used
#' @returns estimated parameters table and AIC
#' @export
summary.s4t_ageclass_model <- function(object, ...) {
  print(object$estimated_parameters,digits = 4)
  # data.frame(AIC = object$AIC)
  print(data.frame(AIC = round(object$AIC,2)))

}


#' ANOVA table for mark-recapture model fits
#'
#' @description
#' Compute analysis of variance table.
#'
#' @param object object of class `s4t_cjs`
#' @param ... object of class `s4t_cjs`
#' @returns An object of class `anova` inheriting form class `data.frame` ???
#'
#' @export
anova.s4t_cjs <- function(object, ...) {
  mCall <- match.call(expand.dots = TRUE)
  dots <- list(...)
  .sapply <- function(L, FUN, ...) unlist(lapply(L, FUN, ...))
  modp <- (as.logical(vapply(dots, is, NA, "s4t_cjs")))

  if (any(modp)) {
    mods <- c(list(object),dots[modp])
    # nobs.vec <- vapply(mods, nobs, 1L) #

    stop("blah")
  }
  stop("blah")
  diff_ll <- abs(mod1$nll - mod2$nll)
  diff_k <- abs(mod1$k - mod2$k)

  ## Should do a check to make sure data is the same (l and m matrices, and age-length stuff)

  c(logLik = diff_ll, p_val = stats::pchisq(q = diff_ll,df = diff_k,lower.tail = FALSE))
}

#' AIC table for ageclass model fits
#'
#' @description
#' Compute analysis of variance table.
#'
#' @param object object of class `s4t_cjs`
#' @param ... object of class `s4t_cjs`
#' @returns An object of class `anova` inheriting form class `data.frame` ???
#'
#' @export
anova.s4t_ageclass_model <- function(object, ...) {
  diff_ll <- abs(mod1$nll - mod2$nll)
  diff_k <- abs(mod1$k - mod2$k)

  ## Should do a check to make sure data is the same (l and m matrices, and age-length stuff)

  c(logLik = diff_ll, p_val = stats::pchisq(q = diff_ll,df = diff_k,lower.tail = FALSE))
}

