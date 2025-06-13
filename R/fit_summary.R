
#' summary of s4t_cjs_ml object
#'
#' @param object A s4t_cjs_ml object
#' @param ... Not used
#' @returns estimated parameters table and AIC
#' @export
summary.s4t_cjs_ml <- function(object, ...) {
  # print(object$estimated_parameters,digits = 4)
  # print(data.frame(AIC = round(object$AIC,2)))

  # row.names(object$estimated_parameters)
  if (object$fit$fixed_age$fixed_age == TRUE) {

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
#' @param pars regular expressions
#' @param probs a numerical vector of quantiles of interest.
#' @param ... Additional arguments to pass to `summary,stanfit-method`.
#'
#' @returns Returns a named list with ...
#' @export
summary.s4t_cjs_rstan <- function(object,
                                  pars = NULL,
                                  probs = c(0.025, 0.25, 0.50, 0.75, 0.975),
                                  ...) {
  # I would like to rename all the parameters in the rstan summary...

  if (is.null(pars)) {
    selected_pars <- rep(TRUE,length(object$original_units$compare_parnames[,"interp_parnames"]))
  } else {
    stopifnot(length(pars) == 1)
    stopifnot(is(pars,"character"))

    selected_pars <- grepl(pars,object$original_units$compare_parnames[,"interp_parnames"])
  }

  s <- rstan::summary(object$res,
                      pars = object$original_units$compare_parnames[,"parnames"],
                      probs = probs, ...)

  rownames(s$summary) <- object$original_units$compare_parnames[,"interp_parnames"]
  dimnames(s$c_summary)$parameter <- object$original_units$compare_parnames[,"interp_parnames"]
  # rownames(s$summary)
  return(s)
  # return(invisible(s))
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
#' @param object object of class `s4t_cjs_ml`
#' @param ... object of class `s4t_cjs_ml`
#' @returns An object of class `anova` inheriting form class `data.frame` ???
#'
#' @export
anova.s4t_cjs_ml <- function(object, ...) {
  mCall <- match.call(expand.dots = TRUE)
  dots <- list(...)
  .sapply <- function(L, FUN, ...) unlist(lapply(L, FUN, ...))
  modp <- (as.logical(vapply(dots, is, NA, "s4t_cjs_ml")))

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
      dplyr::mutate(diff_k = k - dplyr::first(k),
                    diff_logLik = logLik - dplyr::first(logLik),
                    Chi2_stat =  abs(-2 * (logLik - dplyr::first(logLik))),
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
#' @returns An object of class `anova` inheriting from class `data.frame` ??
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
      dplyr::mutate(diff_k = k - dplyr::first(k),
                    diff_logLik = logLik - dplyr::first(logLik),
                    Chi2_stat =  abs(-2 * (logLik - dplyr::first(logLik))),
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



#' Compute and return AIC of fitted model objects
#'
#' @description
#' Compute AIC of one or more fitted `s4t_cjs_ml` model objects
#'
#' @param object A `s4t_cjs_ml` model object
#' @param ... Optionally more fitted model objects
#' @param k, The penalty parameter, taken to be 2. Not used but needed for
#' generic consistency
#' @return If one object is provided, just returns a numeric value with the
#'     corresponding AIC. If more than one is provided, it returns a `data.frame`
#'     with rows corresponding to the objects and columns representing
#'     the number of parameters estimated (`df`), and the AIC
#'
#' @export
AIC.s4t_cjs_ml <- function(object, ..., k = 2) {
  object_list <- list(object, ...)
  if (length(object_list) == 1) {
    object$AIC
  } else {

    object_list_names <- as.character(c(substitute(object), (as.list(substitute(list(...)))[-1])))

    dfs <- sapply(object_list,FUN = function(x) x$k)
    AICs <- sapply(object_list,FUN = function(x) x$AIC)
    data.frame(df = dfs, AIC = AICs,row.names = object_list_names)
  }
}

#' Compute and return AIC of fitted model objects
#'
#' @description
#' Compute AIC of one or more fitted `s4t_ageclass_model` model objects
#'
#' @param object A `s4t_ageclass_model` model object
#' @param ... Optionally more fitted model objects
#' @param k, The penalty parameter, taken to be 2. Not used but needed for
#' generic consistency
#' @return If one object is provided, just returns a numeric value with the
#'     corresponding AIC. If more than one is provided, it returns a `data.frame`
#'     with rows corresponding to the objects and columns representing
#'     the number of parameters estimated (`df`), and the AIC
#'
#' @export
AIC.s4t_ageclass_model <- function(object, ..., k = 2) {
  object_list <- list(object, ...)
  if (length(object_list) == 1) {
    object$AIC
  } else {

    object_list_names <- as.character(c(substitute(object), (as.list(substitute(list(...)))[-1])))

    dfs <- sapply(object_list,FUN = function(x) x$k)
    AICs <- sapply(object_list,FUN = function(x) x$AIC)
    data.frame(df = dfs, AIC = AICs,row.names = object_list_names)
  }
}


#' Extract s4t_ageclass_model object from s4t_cjs_rstan or s4t_cjs_ml
#'
#' Add description
#'
#' @export
#' @param object a s4t_cjs_rstan or s4t_cjs_ml object
#' @returns a s4t_ageclass_model object
extract_ageclass_fit <- function(object) {
  if (object$call$fixed_age == FALSE) {
    stop("Cannot extract ageclass fit when ageclass\n is integrated into the mark-recapture model.")
  }
  object$fit$fixed_age$ageclass_fit
}

# fix no visible binding note
k <- Chi2_stat <-logLik <- logLik <- diff_k <- p_value <- NULL
