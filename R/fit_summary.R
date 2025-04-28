
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
    mods <- list(object,dots[modp])
    # nobs.vec <- vapply(mods, nobs, 1L) #
    if (length(mods) > 2) stop("Not implemented with more than 2 model objects")
    ## CHECK THAT fixed_age has the same value.

    # compare based on
    liks = sapply(mods,FUN = function(x) x$nll)
    ks = c(mods[[1]]$k,mods[[2]]$k)



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
    anova_summary <- anova_summary[2,]

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
#' @param object object of class `s4t_ageclass_model`
#' @param ... object of class `s4t_ageclass_model`
#' @returns An object of class `anova` inheriting from class `data.frame` ???
#'
#' @export
anova.s4t_ageclass_model <- function(object, ...) {
  mCall <- match.call(expand.dots = TRUE)
  dots <- list(...)
  .sapply <- function(L, FUN, ...) unlist(lapply(L, FUN, ...))
  modp <- (as.logical(vapply(dots, is, NA, "s4t_ageclass_model")))

  if (any(modp)) { # compare multiple models
    mods <- list(object,dots[modp])
    # nobs.vec <- vapply(mods, nobs, 1L) #
    if (length(mods) > 2) stop("Not implemented with more than 2 model objects")
    ## CHECK THAT fixed_age has the same value.

    # compare based on
    liks = sapply(mods,FUN = function(x) x$nll)
    ks = c(mods[[1]]$k,mods[[2]]$k)



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
    anova_summary <- anova_summary[2,]

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

