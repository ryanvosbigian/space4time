



#' Plot transition probability estimates from fitted `s4t_cjs_ml` or `s4t_cjs_rstan`objects
#'
#' @description
#' Plot transition probability estimates (which include movement and survival
#'     probabilities) from fitted `s4t_cjs_ml` or `s4t_cjs_rstan`objects.
#'     If the package `geomtextpath` is installed, the figure is more visually appealing.
#'
#'
#' @export
#'
#' @param x a `s4t_cjs` or `s4t_cjs_rstan` object
#' @param textsize an integer for the font size to pass to `geomtextpath::geom_textsegment()`
#' @param ... passed to `dplyr::filter()` for selecting particular survivals. Filters
#'     apparent survivals from `s4t_cjs_ml$apparent_surv` or `s4t_cjs_rstan$apparent_surv`
#' @returns a `ggplot2` figure showing estimated transition probabilities.
#'
#' @examples
#' \dontrun{
#'  sim.dat <- sim_simple_s4t_ch(N = 2000)
#'  m1 <- fit_s4t_cjs_rstan(p_formula = ~ t,
#'                          theta_formula = ~ a1 * a2 * s * j,
#'                          ageclass_formula = ageclass ~ FL,
#'                          fixed_age = TRUE,
#'                          s4t_ch = sim.dat$s4t_ch)
#' plotTransitions(m1)
#' }
#'
plotTransitions <- function(x,textsize = 3, ...) {
  if (!is(x,"s4t_cjs_ml") & !is(x,"s4t_cjs_rstan")) stop("x must be class s4t_cjs or s4t_cjs_rstan")

  cohort_transitions <- as.data.frame(x$cohort_transitions)

  # if ("mean" %in% colnames(cohort_transitions)) {
  #   cohort_transitions$estimate <- cohort_transitions$mean
  # }

  if (is(x,"s4t_cjs_rstan")) {
    cohort_transitions$estimate <- cohort_transitions$mean
  } else {
    cohort_transitions$estimate <- cohort_transitions$estimate
  }

  # suggest geomtextpath and then use geom_textsegment if available
  p <- cohort_transitions %>%
    as.data.frame() %>%
    dplyr::mutate(site_diff = as.integer(as.character(k)) - as.integer(as.character(j)),
           age_diff = as.integer(as.character(a2)) - as.integer(as.character(a1)),
           j = factor(j),
           k = factor(k),
           a1 = factor(a1,levels = rev(unique(a1))),
           a2 = factor(a2,levels = rev(unique(a2))),
           # site_rel = factor(site_rel),
           # site_rec = factor(site_rec),
           # age_rel = factor(age_rel,levels = rev(unique(age_rel))),
           # age_rec = factor(age_rec,levels = rev(unique(age_rec))),
           Theta = round(estimate,2),
           time_rel_label = paste0("Release time: ",s)) %>%
    dplyr::filter(estimate > 0,
                  ...) %>%
    ggplot2::ggplot(ggplot2::aes(x = j, y = a1)) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Sites",
         y = "Age",
         color = "Theta") +
    ggplot2::geom_point(ggplot2::aes(x = j, y = a1),color = "black") +
    ggplot2::geom_point(ggplot2::aes(x = k, y = a1),color = "black")

  if (requireNamespace("geomtextpath", quietly = TRUE)) {
    p +
      geomtextpath::geom_textsegment(ggplot2::aes(x = j,xend = k,
                                         y = a1, yend = a2,
                                         color = estimate,
                                         label = Theta),hjust = 0.25,
                                     size = textsize) +
      ggplot2::facet_wrap(~time_rel_label)
  } else {
    p +
      ggplot2::geom_segment(ggplot2::aes(x = j,xend = k,
                                y = a1, yend = a2,
                                color = estimate)) +
      ggplot2::geom_text(ggplot2::aes(label = Theta,
                             hjust = -site_diff,
                             vjust = -0.1 + 1.25*age_diff

      )) +
      ggplot2::facet_wrap(~time_rel_label)
  }


}




#' Plot apparent survival estimates from fitted `s4t_cjs_ml` or `s4t_cjs_rstan`objects
#'
#' @description
#' Plot apparent survival estimates from fitted `s4t_cjs_ml` or `s4t_cjs_rstan`objects.
#'     If the package `geomtextpath` is installed, the figure is more visually appealing.
#'
#'
#' @export
#'
#' @param x a `s4t_cjs` or `s4t_cjs_rstan` object
#' @param textsize an integer for the font size to pass to `geomtextpath::geom_textsegment()`
#' @param ... passed to `dplyr::filter()` for selecting particular survivals. Filters
#'     apparent survivals from `s4t_cjs_ml$apparent_surv` or `s4t_cjs_rstan$apparent_surv`
#' @returns a `ggplot2` figure showing estimated apparent survival
#'
#' @examples
#' \dontrun{
#'  sim.dat <- sim_simple_s4t_ch(N = 2000)
#'  m1 <- fit_s4t_cjs_rstan(p_formula = ~ t,
#'                          theta_formula = ~ a1 * a2 * s * j,
#'                          ageclass_formula = ageclass ~ FL,
#'                          fixed_age = TRUE,
#'                          s4t_ch = sim.dat$s4t_ch)
#' plotSurvival(m1)
#' }
#'
plotSurvival <- function(x,textsize = 3, ...) {

  if (!is(x,"s4t_cjs_ml") & !is(x,"s4t_cjs_rstan")) stop("x must be class s4t_cjs or s4t_cjs_rstan")

  apparent_surv <- as.data.frame(x$apparent_surv)

  if (is(x,"s4t_cjs_rstan")) {
    apparent_surv$estimate <- apparent_surv$mean
  } else {
    apparent_surv$estimate <- apparent_surv$estimate
  }

  # suggest geomtextpath and then use geom_textsegment if available
  p <- apparent_surv %>%
    as.data.frame() %>%
    dplyr::mutate(site_diff = as.integer(as.character(k)) - as.integer(as.character(j)),
                  # age_diff = as.integer(as.character(a2)) - as.integer(as.character(a1)),
                  j = factor(j),
                  k = factor(k),
                  a1 = factor(a1,levels = rev(unique(a1))),
                  # a2 = factor(a2),
                  # site_rel = factor(site_rel),
                  # site_rec = factor(site_rec),
                  # age_rel = factor(age_rel,levels = rev(unique(age_rel))),
                  # age_rec = factor(age_rec,levels = rev(unique(age_rec))),
                  Theta = round(estimate,2),
                  time_rel_label = paste0("Release time: ",s)) %>%
    # dplyr::mutate(site_diff = as.integer(as.character(k)) - as.integer(as.character(j)),
    #        # age_diff = a1,
    #        j = factor(j),
    #        k = factor(k),
    #        a1 = factor(a1,levels = rev(unique(a1))),
    #        Theta = round(estimate,2)) %>%
    dplyr::filter(estimate > 0,
                  ...) %>%
    ggplot2::ggplot(ggplot2::aes(x = j, y = (a1))) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Sites",
                 y = "Age",
                 color = "Survival") +
    ggplot2::geom_point(ggplot2::aes(x = j, y = a1),color = "black") +
    ggplot2::geom_point(ggplot2::aes(x = k, y = a1),color = "black")

  if (requireNamespace("geomtextpath", quietly = TRUE)) {
    p +
      geomtextpath::geom_textsegment(ggplot2::aes(x = j,xend = k,
                                         y = a1,yend = a1,
                                         color = estimate,
                                         label = Theta),
                                     size = textsize) +
      ggplot2::facet_wrap(~time_rel_label)
  } else {
    message("Install package geomtextpath for better figures")
    p +
      ggplot2::geom_segment(ggplot2::aes(x = j,xend = k,
                                y = a1,
                                color = estimate)) +
      ggplot2::geom_text(ggplot2::aes(label = Theta,
                             hjust = -site_diff
                             # vjust = -0.1 + 1.25*age_diff

      )) +
      ggplot2::facet_wrap(~time_rel_label)
  }


}




#' Plot summary of s4t_cjs_ch object
#'
#' add description
#'
#'
#' @param x a `s4t_cjs_ch` object
#' @param ... not used
#' @returns a ggplot2 type figure of the observed transition
plotCH <- function(x, ...) {
  if (!is(x,"s4t_cjs_ch")) stop("x must be class s4t_cjs_ch")

  return(plot(1:5,2:6))

#
#   s4t_ch$obs_data$obs_ch %>%
#     as.data.frame() %>%
#     dplyr::mutate(ageclass = s4t_ch$obs_data$obs_aux$ageclass,
#                   obs_time = s4t_ch$obs_data$obs_aux$obs_time) %>%
#
#
#
#   cohort_surv <- as.data.frame(x$cohort_surv)
#
#   # if ("mean" %in% colnames(cohort_surv)) {
#   #   cohort_surv$estimate <- cohort_surv$mean
#   # }
#
#   if (is(x,"s4t_cjs_rstan")) {
#     cohort_surv$estimate <- cohort_surv$mean
#   } else {
#     cohort_surv$estimate <- cohort_surv$estimate
#   }
#
#   # suggest geomtextpath and then use geom_textsegment if available
#   p <- cohort_surv %>%
#     as.data.frame() %>%
#     dplyr::mutate(site_diff = as.integer(as.character(k)) - as.integer(as.character(j)),
#                   age_diff = as.integer(as.character(a2)) - as.integer(as.character(a1)),
#                   # j = factor(j),
#                   # k = factor(k),
#                   site_rel = factor(site_rel),
#                   site_rec = factor(site_rec),
#                   age_rel = factor(age_rel,levels = rev(unique(age_rel))),
#                   age_rec = factor(age_rec,levels = rev(unique(age_rec))),
#                   Theta = round(estimate,2),
#                   time_rel_label = paste0("Release time: ",time_rel)) %>%
#     dplyr::filter(estimate > 0,
#                   ...) %>%
#     ggplot2::ggplot(aes(x = site_rel, y = age_rel)) +
#     ggplot2::theme_bw() +
#     ggplot2::labs(x = "Sites",
#                   y = "Age") +
#     geom_point(aes(x = site_rel, y = age_rel),color = "black") +
#     geom_point(aes(x = site_rec, y = age_rel),color = "black")
#
#   if (requireNamespace("geomtextpath", quietly = TRUE)) {
#     p +
#       geomtextpath::geom_textsegment(aes(x = site_rel,xend = site_rec,
#                                          y = age_rel, yend = age_rec,
#                                          color = estimate,
#                                          label = Theta),hjust = 0.25) +
#       ggplot2::facet_wrap(~time_rel_label)
#   } else {
#     message("Install package geomtextpath for better figures")
#
#     p +
#       ggplot2::geom_segment(aes(x = site_rel,xend = site_rec,
#                                 y = age_rel, yend = age_rec,
#                                 color = estimate)) +
#       ggplot2::geom_text(aes(label = Theta,
#                              hjust = -site_diff,
#                              vjust = -0.1 + 1.25*age_diff
#
#       )) +
#       ggplot2::facet_wrap(~time_rel_label)
#   }
#

}

#' Create traceplots from `s4t_cjs_rstan` object
#'
#' add description
#' @export
#'
#' @importFrom rstan traceplot
#'
#' @param object a `s4t_cjs_rstan` object
#' @param pars a `character` object that parses parameters using regular
#'      expressions
#' @param ... not used
#' @returns a ggplot2 object of traceplots
traceplot.s4t_cjs_rstan <- function(object,
                      pars = NULL,
                      ...) {
  stopifnot(is(object,"s4t_cjs_rstan"))

  object$res@sim$fnames_oi

  compare_parnames <- object$original_units$compare_parnames

  if (is.null(pars)) {
    keep <- rep(FALSE, nrow(compare_parnames))
    if (length(keep) < 10) {
      keep[1:length(keep)] <- TRUE
    } else {
      keep[1:10] <- TRUE
    }


  } else {
    stopifnot(is(pars,"character"))
    stopifnot(length(pars) == 1)



    keep <- grepl(pars,compare_parnames[,"interp_parnames"])




  }
  kept_parnames <- as.character(compare_parnames[keep,2])
  names(kept_parnames) <- compare_parnames[keep,1]


  p <- rstan::traceplot(object = object$res,
                   pars = compare_parnames[keep,"parnames"],
                   ...)

  p$facet$params$labeller <- ggplot2::as_labeller(kept_parnames)

  p

}

# fix no visible binding note
b <- k <- j <- r <- site_rel <- site_rec <- age_rel <- estimate <- time_rel <-
  Theta <- site_diff <- age_rec <- age_diff <- NULL
