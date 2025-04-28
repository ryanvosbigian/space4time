



#' Plot s4t_cjs_rstan object results
plotTheta <- function(x, ...) {

  # suggest geomtextpath and then use geom_textsegment if available

  x$cohort_surv %>%
    as.data.frame() %>%
    filter(!(mean == 0)) %>%
    mutate(site_diff = k - j,
           age_diff = a2 - a1,
           j = factor(j),
           k = factor(k),
           a1 = factor(a1,levels = rev(unique(a1))),
           a2 = factor(a2),
           Theta = printable_round(mean,2)) %>%
    ggplot2::ggplot(aes(x = j, y = (a1))) +
    ggplot2::geom_segment(aes(x = j,xend = k,
                              y = a1, yend = a2,
                              color = mean)) +
    ggplot2::geom_text(aes(label = Theta,
                           hjust = -site_diff,
                           vjust = -0.1 + 1.25*age_diff

                       )) +
    facet_wrap(~s) +
    theme_bw() +
    labs(x = "Sites",
         y = "Age")



}




#' Plot s4t_cjs_rstan object results
plotSurvival <- function(x, ...) {

  # suggest geomtextpath and then use geom_textsegment if available

  x$overall_surv %>%
    as.data.frame() %>%
    filter(!(mean == 0)) %>%
    mutate(site_diff = k - j,
           # age_diff = a1,
           j = factor(j),
           k = factor(k),
           a1 = factor(a1,levels = rev(unique(a1))),
           Theta = printable_round(mean,2)) %>%
    ggplot2::ggplot(aes(x = j, y = (a1))) +
    ggplot2::geom_segment(aes(x = j,xend = k,
                              y = a1,
                              color = mean)) +
    ggplot2::geom_text(aes(label = Theta,
                           hjust = -site_diff
                           # vjust = -0.1 + 1.25*age_diff

    )) +
    facet_wrap(~t) +
    theme_bw() +
    labs(x = "Sites",
         y = "Age")



}
