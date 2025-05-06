



#' Plot s4t_cjs or s4t_cjs_rstan object results
#'
#' add description
#' @export
#'
#' @param x a `s4t_cjs` or `s4t_cjs_rstan` object
#' @param ... not used
#' @returns a ggplot2 type figure of the transition probabilities
plotTheta <- function(x, ...) {
  if (!is(x,"s4t_cjs") & !is(x,"s4t_cjs_rstan")) stop("x must be class s4t_cjs or s4t_cjs_rstan")

  cohort_surv <- as.data.frame(x$cohort_surv)

  # if ("mean" %in% colnames(cohort_surv)) {
  #   cohort_surv$estimate <- cohort_surv$mean
  # }

  if (is(x,"s4t_cjs_rstan")) {
    cohort_surv$estimate <- cohort_surv$mean
  } else {
    cohort_surv$estimate <- cohort_surv$estimate
  }

  # suggest geomtextpath and then use geom_textsegment if available
  p <- cohort_surv %>%
    as.data.frame() %>%
    dplyr::mutate(site_diff = as.integer(as.character(k)) - as.integer(as.character(j)),
           age_diff = as.integer(as.character(a2)) - as.integer(as.character(a1)),
           # j = factor(j),
           # k = factor(k),
           site_rel = factor(site_rel),
           site_rec = factor(site_rec),
           age_rel = factor(age_rel,levels = rev(unique(age_rel))),
           age_rec = factor(age_rec,levels = rev(unique(age_rec))),
           Theta = round(estimate,2),
           time_rel_label = paste0("Release time: ",time_rel)) %>%
    dplyr::filter(estimate > 0,
                  ...) %>%
    ggplot2::ggplot(ggplot2::aes(x = site_rel, y = age_rel)) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Sites",
         y = "Age",
         color = "Theta") +
    ggplot2::geom_point(ggplot2::aes(x = site_rel, y = age_rel),color = "black") +
    ggplot2::geom_point(ggplot2::aes(x = site_rec, y = age_rel),color = "black")

  if (requireNamespace("geomtextpath", quietly = TRUE)) {
    p +
      geomtextpath::geom_textsegment(ggplot2::aes(x = site_rel,xend = site_rec,
                                         y = age_rel, yend = age_rec,
                                         color = estimate,
                                         label = Theta),hjust = 0.25) +
      ggplot2::facet_wrap(~time_rel_label)
  } else {
    p +
      ggplot2::geom_segment(ggplot2::aes(x = site_rel,xend = site_rec,
                                y = age_rel, yend = age_rec,
                                color = estimate)) +
      ggplot2::geom_text(ggplot2::aes(label = Theta,
                             hjust = -site_diff,
                             vjust = -0.1 + 1.25*age_diff

      )) +
      ggplot2::facet_wrap(~time_rel_label)
  }


}




#' Plot s4t_cjs or s4t_cjs_rstan object results
#'
#' add description
#' @export
#'
#' @param x a `s4t_cjs` or `s4t_cjs_rstan` object
#' @param ... not used
#' @returns a ggplot2 type figure of the apparent survivals
plotSurvival <- function(x, ...) {

  if (!is(x,"s4t_cjs") & !is(x,"s4t_cjs_rstan")) stop("x must be class s4t_cjs or s4t_cjs_rstan")

  overall_surv <- as.data.frame(x$overall_surv)

  if (is(x,"s4t_cjs_rstan")) {
    overall_surv$estimate <- overall_surv$mean
  } else {
    overall_surv$estimate <- overall_surv$estimate
  }

  # suggest geomtextpath and then use geom_textsegment if available
  p <- overall_surv %>%
    as.data.frame() %>%
    dplyr::mutate(site_diff = as.integer(as.character(k)) - as.integer(as.character(j)),
                  # age_diff = as.integer(as.character(a2)) - as.integer(as.character(a1)),
                  j = factor(j),
                  k = factor(k),
                  site_rel = factor(site_rel),
                  site_rec = factor(site_rec),
                  age_rel = factor(age_rel,levels = rev(unique(age_rel))),
                  # age_rec = factor(age_rec,levels = rev(unique(age_rec))),
                  Theta = round(estimate,2),
                  time_rel_label = paste0("Release time: ",time_rel)) %>%
    # dplyr::mutate(site_diff = as.integer(as.character(k)) - as.integer(as.character(j)),
    #        # age_diff = a1,
    #        j = factor(j),
    #        k = factor(k),
    #        a1 = factor(a1,levels = rev(unique(a1))),
    #        Theta = round(estimate,2)) %>%
    dplyr::filter(estimate > 0,
                  ...) %>%
    ggplot2::ggplot(ggplot2::aes(x = site_rel, y = (age_rel))) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Sites",
                 y = "Age",
                 color = "Survival") +
    ggplot2::geom_point(ggplot2::aes(x = site_rel, y = age_rel),color = "black") +
    ggplot2::geom_point(ggplot2::aes(x = site_rec, y = age_rel),color = "black")

  if (requireNamespace("geomtextpath", quietly = TRUE)) {
    p +
      geomtextpath::geom_textsegment(ggplot2::aes(x = site_rel,xend = site_rec,
                                         y = age_rel,yend = age_rel,
                                         color = estimate,
                                         label = Theta)) +
      ggplot2::facet_wrap(~time_rel_label)
  } else {
    message("Install package geomtextpath for better figures")
    p +
      ggplot2::geom_segment(ggplot2::aes(x = site_rel,xend = site_rec,
                                y = age_rel,
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
#' @export
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
#' @param object a `s4t_cjs_rstan` object
#' @param pars a `character` object that parses parameters using regular
#'      expressions
#' @param ... not used
#' @returns a ggplot2 object of traceplots
traceplot <- function(object,
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
k <- j <- site_rel <- site_rec <- age_rel <- estimate <- time_rel <-
  Theta <- site_diff <- age_rec <- age_diff <- NULL
