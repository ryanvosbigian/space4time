
#' summary of s4t_cjs object
#'
#' @param object A s4t_cjs object
#' @param ... Not used
#' @returns estimated parameters table and AIC
#' @export
summary.s4t_cjs <- function(object, ...) {
  # print(object$estimated_parameters,digits = 4)
  # print(data.frame(AIC = round(object$AIC,2)))

  # row.names(object$estimated_parameters)
  if (object$call$fixed_age == TRUE) {

    coefficients <- list(cjs = object$estimated_parameters,
                         age = object$fit$fixed_age$ageclass_fit$estimated_parameters)
  } else {
    coefficients <- list(all = object$estimated_parameters)
  }


  summary_list <- list(
    call = object$call,
    terms = object$estimated_parameters$parameter,
    coefficients = coefficients,
    AIC = object$AIC,
    vcov = object$vcov
  )

  new_summary_list <- structure(summary_list, class = paste("summary", class(object), sep = "."))
  new_summary_list


  # return(list(estimated_parameters = object$estimated_parameters,
  #        AIC = data.frame(AIC = round(object$AIC,2))))
}


#' print summary of s4t_cjs
#'
#' @export
#' @param x summary.s4t_cjs object
#' @param digits number of digits
#' @param ... passed to
#'
print.summary.s4t_cjs <- function(x,
                                  digits = max(3L, getOption("digits") - 3L),
                                  ...) {

  # pasting the formula call
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")

  # pasting the residual summary
  cat("\nAIC:\n")
  resAIC <- x$AIC
  # names(resAIC) <- ""
  print(resAIC)

  if (length(x$coefficients) == 1) {
    cat("\nCoefficients (all):\n")
    coefs_all <- as.matrix(x$coefficients$all[,2:ncol(x$coefficients$all)])
    # colnames(coefs_fixed) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    colnames(coefs_all)[1:3] <- c("Estimate", "Std. Error", "z value")
    printCoefmat(coefs_all, digits = digits, na.print = "NA",P.values = FALSE,has.Pvalue = FALSE, ...)

  } else {
    cat("\nCoefficients (age class):\n")
    coefs_age <- as.matrix(x$coefficients$age[,2:ncol(x$coefficients$age)])
    # colnames(coefs_fixed) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    colnames(coefs_age)[1:3] <- c("Estimate", "Std. Error", "z value")
    printCoefmat(coefs_age, digits = digits, na.print = "NA",P.values = FALSE,has.Pvalue = FALSE, ...)

    cat("\nCoefficients (mark-recapture):\n")
    coefs_cjs <- as.matrix(x$coefficients$cjs[,2:ncol(x$coefficients$cjs)])
    # colnames(coefs_fixed) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    colnames(coefs_cjs)[1:3] <- c("Estimate", "Std. Error", "z value")
    printCoefmat(coefs_cjs, digits = digits, na.print = "NA",P.values = FALSE,has.Pvalue = FALSE, ...)
  }


  invisible(x)
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
  # I would like to rename all the parameters in the rstan summary...

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

  if (any(modp)) { # compare multiple models
    mods <- c(list(object),dots[modp])
    # nobs.vec <- vapply(mods, nobs, 1L) #
    if (length(mods) > 2) stop("Not implemented with more than 2 model objects")
    ## CHECK THAT fixed_age has the same value.

    # compare based on
    liks = sapply(mods,FUN = function(x) x$logLik)
    ks = sapply(mods,FUN = function(x) x$k)



    if (ks[2] < ks[1]) {
      full_name <- deparse(substitute(object)) # replace as.character with deparse
      reduced_name <- as.character(as.list(substitute(list(...)))[-1])
    } else {
      reduced_name <- deparse(substitute(object)) # replace as.character with deparse
      full_name <- as.character(as.list(substitute(list(...)))[-1])
    }



    anova_summary <- data.frame(ord = 1:length(mods),
                           logLik = liks,
                           k = ks)
    anova_summary <- anova_summary %>%
      dplyr::arrange(k) %>%
      dplyr::mutate(diff_k = k - first(k),
                    diff_logLik = logLik - first(logLik),
                    Chi2_stat =  abs(-2 * (logLik - first(logLik))),
                    p_value = stats::pchisq(Chi2_stat, diff_k, lower.tail = FALSE)) %>%
      dplyr::select(Df = diff_k,
                    Chi2 = Chi2_stat,
                    p.value = p_value) %>%
      as.data.frame()

    rownames(anova_summary) <- paste(full_name, "vs", reduced_name)
    attr(anova_summary, "full") <- full_name
    attr(anova_summary, "reduced") <- reduced_name

    anova_summary <- structure(anova_summary, heading = c("Likelihood Ratio Test\n", paste("Response:", deparse(object$formula[[2L]]))))

  } else {
    # object
    stop("Not yet implemented")

  }


  ## Should do a check to make sure data is the same (l and m matrices, and age-length stuff)

  structure(anova_summary, class = c(paste("anova", class(object), sep = "."), "data.frame"))
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

