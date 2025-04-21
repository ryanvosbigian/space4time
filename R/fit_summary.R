
#' summary of s4t_cjs object
#'
#' @param object A s4t_cjs object
#' @returns estimated parameters table and AIC
#' @export
summary.s4t_cjs <- function(object, ...) {
  print(object$estimated_parameters,digits = 4)
  print(data.frame(AIC = round(object$AIC,2)))
  return(list(estimated_parameters = object$estimated_parameters),
         AIC = data.frame(AIC = round(object$AIC,2)))
}

#' summary of s4t_ageclass_cjs object
#'
#' @param object A s4t_ageclass_cjs object
#' @returns estimated parameters table and AIC
#' @export
summary.s4t_ageclass_model <- function(object, ...) {
  print(object$estimated_parameters,digits = 4)
  # data.frame(AIC = object$AIC)
  print(data.frame(AIC = round(object$AIC,2)))

}

#' @export
anova.s4t_cjs <- function(mod1,mod2, ...) {
  diff_ll <- abs(mod1$nll - mod2$nll)
  diff_k <- abs(mod1$k - mod2$k)

  ## Should do a check to make sure data is the same (l and m matrices, and age-length stuff)

  c(logLik = diff_ll, p_val = stats::pchisq(q = diff_ll,df = diff_k,lower.tail = FALSE))
}

#' @export
anova.s4t_ageclass_model <- function(mod1,mod2, ...) {
  diff_ll <- abs(mod1$nll - mod2$nll)
  diff_k <- abs(mod1$k - mod2$k)

  ## Should do a check to make sure data is the same (l and m matrices, and age-length stuff)

  c(logLik = diff_ll, p_val = stats::pchisq(q = diff_ll,df = diff_k,lower.tail = FALSE))
}

