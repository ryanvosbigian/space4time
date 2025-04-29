

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


print.summary.s4t_cjs_rstan <- function(x,
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


#' print summary of s4t_cjs_rstan
#'
#' @export
#' @param x s4t_cjs_rstan object
#' @param pars a character vector of parameter names. Default is the parameters to save,
#'     but if `include == FALSE`, then the specified parameters are excluded.
#' @param probs number of digits.
#' @param digits_summary the number of significant digits to use when printing the summary
#' @param include logical scalar indicating whether to include the parameters named by the
#'      pars argument.
#' @param ... passed to `summary` method for s4t_cjs_rstan object.
#'
print.s4t_cjs_rstan <- function(x,
                                pars = x$res@sim$pars_oi,
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

  s <- rstan::summary(obj, pars, probs, ...)


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
