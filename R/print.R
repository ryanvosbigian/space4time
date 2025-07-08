

#' print summary of s4t_cjs_ml
#'
#' @export
#' @param x summary.s4t_cjs_ml object
#' @param digits number of digits
#' @param ... passed to
#'
print.summary.s4t_cjs_ml <- function(x,
                                  digits = max(3L, getOption("digits") - 3L),
                                  ...) {

  # pasting the formula call
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")

  # pasting the residual summary
  cat("\nAIC:\n")
  resAIC <- x$AIC
  names(resAIC) <- NULL
  print(resAIC)

  if (length(x$coefficients) == 1) {
    cat("\nCoefficients (all):\n")
    coefs_all <- as.matrix(x$coefficients$all[,2:ncol(x$coefficients$all)])
    # colnames(coefs_fixed) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    colnames(coefs_all)[1:3] <- c("Estimate", "Std. Error", "z value")
    stats::printCoefmat(coefs_all, digits = digits, na.print = "NA",P.values = FALSE,has.Pvalue = FALSE, ...)

  } else {
    cat("\nCoefficients (age class):\n")
    coefs_age <- as.matrix(x$coefficients$age[,2:ncol(x$coefficients$age)])
    # colnames(coefs_fixed) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    colnames(coefs_age)[1:3] <- c("Estimate", "Std. Error", "z value")
    stats::printCoefmat(coefs_age, digits = digits, na.print = "NA",P.values = FALSE,has.Pvalue = FALSE, ...)

    cat("\nCoefficients (mark-recapture):\n")
    coefs_cjs <- as.matrix(x$coefficients$cjs[,2:ncol(x$coefficients$cjs)])
    # colnames(coefs_fixed) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    colnames(coefs_cjs)[1:3] <- c("Estimate", "Std. Error", "z value")
    stats::printCoefmat(coefs_cjs, digits = digits, na.print = "NA",P.values = FALSE,has.Pvalue = FALSE, ...)
  }


  invisible(x)
}


#' print s4t_cjs_ml
#'
#' @export
#' @param x s4t_cjs_ml object
#' @param ... for generic consistency
print.s4t_cjs_ml <- function(x, ...) {
  stopifnot(is(x,"s4t_cjs_ml"))
  summary(x)
}


# print.summary.s4t_cjs_rstan <- function(x,
#                                         digits = max(3L, getOption("digits") - 3L),
#                                         ...) {
#   # pasting the formula call
#   cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
#
#   # pasting the residual summary
#   cat("\nAIC:\n")
#   resAIC <- x$AIC
#   names(resAIC) <- NULL
#   print(resAIC)
#
#   if (length(x$coefficients) == 1) {
#     cat("\nCoefficients (all):\n")
#     coefs_all <- as.matrix(x$coefficients$all[,2:ncol(x$coefficients$all)])
#     # colnames(coefs_fixed) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
#     colnames(coefs_all)[1:3] <- c("Estimate", "Std. Error", "z value")
#     stats::printCoefmat(coefs_all, digits = digits, na.print = "NA",P.values = FALSE,has.Pvalue = FALSE, ...)
#
#   } else {
#     cat("\nCoefficients (age class):\n")
#     coefs_age <- as.matrix(x$coefficients$age[,2:ncol(x$coefficients$age)])
#     # colnames(coefs_fixed) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
#     colnames(coefs_age)[1:3] <- c("Estimate", "Std. Error", "z value")
#     stats::printCoefmat(coefs_age, digits = digits, na.print = "NA",P.values = FALSE,has.Pvalue = FALSE, ...)
#
#     cat("\nCoefficients (mark-recapture):\n")
#     coefs_cjs <- as.matrix(x$coefficients$cjs[,2:ncol(x$coefficients$cjs)])
#     # colnames(coefs_fixed) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
#     colnames(coefs_cjs)[1:3] <- c("Estimate", "Std. Error", "z value")
#     stats::printCoefmat(coefs_cjs, digits = digits, na.print = "NA",P.values = FALSE,has.Pvalue = FALSE, ...)
#   }
#
#
#   invisible(x)
# }


#' print summary of s4t_cjs_rstan
#'
#' @export
#' @param x s4t_cjs_rstan object
#' @param pars regular expressions ... NVM go back to rstan
#' @param probs number of digits.
#' @param digits_summary the number of significant digits to use when printing the summary
#' @param include logical scalar indicating whether to include the parameters named by the
#'      pars argument.
#' @param ... passed to `summary` method for s4t_cjs_rstan object.
#'
print.s4t_cjs_rstan <- function(x,
                                pars = NULL,
                                probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                                digits_summary = 2,
                                include = TRUE,
                                ...) {
  obj <- x$res



  if (obj@mode == 1L) {
    cat("Stan model '", obj@model_name, "' is of mode 'test_grad';\n",
        "sampling is not conducted.\n", sep = "")
    return(invisible(NULL))
  } else if (obj@mode == 2L) {
    cat("Stan model '", obj@model_name, "' does not contain samples.\n",
        sep = "")
    return(invisible(NULL))
  }

  if (!include) pars <- setdiff(obj@sim$pars_oi, pars)

  s <- summary(x, pars, probs, ...)


  if (is.null(s)) return(invisible(NULL))

  n_kept <- obj@sim$n_save - obj@sim$warmup2
  cat("Inference for Stan model: ", obj@model_name, ".\n", sep = "")
  cat(obj@sim$chains, " chains, each with iter=", obj@sim$iter,
      "; warmup=", obj@sim$warmup, "; thin=", obj@sim$thin, "; \n",
      "post-warmup draws per chain=", n_kept[1], ", ", "total post-warmup draws=",
      sum(n_kept), ".\n\n", sep = "")

  s$summary[, "n_eff"] <- round(s$summary[, "n_eff"], 0)


  print(round(s$summary, digits_summary), ...)

  sampler <- attr(obj@sim$samples[[1]], "args")$sampler_t

  if (!is.null(obj@stan_args[[1]]$method) && isTRUE(obj@stan_args[[1]]$method ==
                                                  "variational")) {
    if ("diagnostics" %in% names(obj@sim) & "ir_idx" %in% names(obj@sim$diagnostics) &
        !is.null(obj@sim$diagnostics$ir_idx)) {
      cat("\nApproximate samples were drawn using VB(",
          obj@stan_args[[1]]$algorithm, ") + PSIS at ", obj@date,
          ".\n", sep = "")
    } else {
      cat("\nApproximate samples were drawn using VB(",
          obj@stan_args[[1]]$algorithm, ") at ", obj@date,
          ".\n", sep = "")
      message("We recommend genuine 'sampling' from the posterior distribution for final inferences!")
    }
    return(invisible(NULL))
  } else {
    cat("\nSamples were drawn using ", sampler, " at ", obj@date,
        ".\n", "For each parameter, n_eff is a crude measure of effective sample size,\n",
        "and Rhat is the potential scale reduction factor on split chains (at \n",
        "convergence, Rhat=1).\n", sep = "")
    return(invisible(NULL))
  }
}

#' Print summary of `clean_s4t_ch` object
#'
#' Returns the observations of individuals that had observations dropped
#'     for data cleaning.
#'
#' @export
#' @param x `clean_s4t_ch` object from `clean_s4t_ch()`
#' @param ... passed to `print.tibble()`
#'
#' @export
print.clean_s4t_ch <- function(x, ...) {
  stopifnot(is(x,"clean_s4t_ch"))

  x$intermediate_ch_df %>%
    as.data.frame() %>%
    dplyr::group_by(id) %>%
    dplyr::filter(TRUE %in% drop_obs) %>%
    dplyr::arrange(id) %>%
    print(,...)

  # return(invisible(NULL))
}



drop_obs <- NULL



#' Print summary of `cov_s4t_ch` object
#'
#' Returns the indices
#'
#'
#' @param x `s4t_ch` object
#' @param ... passed to `print.tibble()`
#'
#' @export
print.cov_s4t_ch <- function(x, ...) {
  stopifnot(is(x,"cov_s4t_ch"))

  cat("Indices (and covariates) for state transitions:\n")
  x$indices_theta %>%
    as.data.frame() %>%
    print(,...)

  cat("\nIndices (and covariates) for detection probability:\n")
  x$indices_p_obs %>%
    as.data.frame() %>%
    print(,...)

  return(invisible(NULL))
}

#' Print summary of `s4t_config` object
#'
#' Print summary of `s4t_config` object.
#'
#' @export
#' @param x `s4t_config` object
#' @param ... passed to `print()`
#'
#' @export
print.s4t_config <- function(x, ...) {
  stopifnot(is(x,"s4t_config"))

  cat("Site and age transition configuration object\n")
  cat("\nThere are N = ",length(x$sites_names)," with N = ",
      sum(x$holdover_config)," sites with holdovers\n",sep="")
  cat("\nSites: ",paste0(x$sites_names,collapse = ", "),"\n",sep="")

  cat("\nSites with holdovers: ",
      paste0(x$sites_names[which(rowSums(x$holdover_config) > 0)],collapse = ", "),
      "\n",
      sep="")

  if (!is.null(x$sites_to_pool)) {
    cat("\nSites pooled:\n")
    for (i in 1:length(x$sites_to_pool)) {
      cat(names(x$sites_to_pool)[i]," include: ",paste0(x$sites_to_pool[[i]],collapse =", "),"\n",
          sep = "")
    }

  }

  cat("\nSite -> site:\n")

  cat(paste0(rownames(x$sites_config)," -> ",
             apply(x$sites_config,1,function(y) ifelse(length(colnames(x$sites_config)[which(y == 1)]) == 1,
                                                       colnames(x$sites_config)[which(y == 1)],
                                                       ""
                                                       )
                   ),
             "\n",collapse = "")
      )



  cat("\nAge range per site:\n")
  cat(paste0(x$sites_names,": ",x$obs_min_a,"-",x$obs_max_a,"\n",collapse = ""))

  return(invisible(NULL))
}



#' Print summary of `s4t_ch` object
#'
#' Print summary of `s4t_ch` object.
#'
#'
#' @param x `s4t_ch` object
#' @param ... passed to `print()`
#'
#'
print.s4t_ch <- function(x, ...) {
  stopifnot(is(x,"s4t_ch"))

  cat("Capture history object\n")
  cat("\nThere are N = ",length(x$sites_names)," with N = ",
      sum(x$holdover_config)," sites with holdovers\n",sep="")
  cat("\nSites: ",paste0(x$sites_names,collapse = ", "),"\n",sep="")

  cat("\nSites with holdovers: ",
      paste0(x$sites_names[which(rowSums(x$holdover_config) > 0)],collapse = ", "),
      "\n",
      sep="")

  if (!is.null(x$sites_to_pool)) {
    cat("\nSites pooled:\n")
    for (i in 1:length(x$sites_to_pool)) {
      cat(names(x$sites_to_pool)[i]," include: ",paste0(x$sites_to_pool[[i]],collapse =", "),"\n",
          sep = "")
    }

  }

  cat("\nSite -> site:\n")

  cat(paste0(rownames(x$sites_config)," -> ",
             apply(x$sites_config,1,function(y) ifelse(length(colnames(x$sites_config)[which(y == 1)]) == 1,
                                                       colnames(x$sites_config)[which(y == 1)],
                                                       ""
             )
             ),
             "\n",collapse = "")
  )



  cat("\nAge range per site:\n")
  cat(paste0(x$sites_names,": ",x$obs_min_a,"-",x$obs_max_a,"\n",collapse = ""))

  return(invisible(NULL))
}

