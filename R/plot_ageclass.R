
#' Plot confusion matrix for fit_ageclass fits
#'
#' add description. Currently only does confusion matrix.
#'
#' @export
#'
#' @param x a `s4t_ageclass_model` object
#' @param param not used.
#' @param ... Not used.
#' @returns ggplot2 style figure
plot.s4t_ageclass_model <- function(x, param,...) {


  # check that param is a character and length 1

  info_ageclass <- ageclass_call(age_formula = x$call$age_formula,
                                 obs_aux = x$s4t_ch$ch$obs_aux,#max_a = max(x$s4t_ch$ch_info$set_max_a),
                                 ll = FALSE)

  prob <- ageclass_nll(par = x$res$par,
                       max_a = max(x$s4t_ch$ch_info$set_max_a),
                       mod_mat_a_beta = info_ageclass$mod_mat_a_beta,
                       ll = FALSE)

  min_a <- x$s4t_ch$ch_info$observed_relative_min_max$min_obs_age
  max_a <-  x$s4t_ch$ch_info$observed_relative_min_max$max_obs_age

  age_cols <- paste0("Age",min_a:(ncol(prob) + min_a - 1))

  colnames(prob) = age_cols

  pred_age <- apply(prob,MARGIN = 1,which.max) + min_a - 1

  compare_obs_pred <- x$s4t_ch$ch$obs_aux %>%
    cbind(prob) %>%
    cbind(pred_ageclass = pred_age) %>%
    as.data.frame() %>%
    dplyr::filter(!is.na(ageclass))


  # plot confusion matrix
  p1 <- suppressWarnings(suppressMessages(compare_obs_pred %>%
    dplyr::group_by(ageclass) %>%
    dplyr::mutate(count_ageclass = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(pred_ageclass) %>%
    dplyr::mutate(count_predageclass = dplyr::n()) %>%
    dplyr::group_by(ageclass,pred_ageclass) %>%
    dplyr::summarize(count = dplyr::n(),
                     Proportion.observed = count / count_ageclass,
                     Proportion.predicted = count / count_predageclass) %>%
    tidyr::pivot_longer(cols = c("Proportion.observed","Proportion.predicted"),names_to = "Prop",values_to = "val") %>%
    ggplot2::ggplot(ggplot2::aes(ageclass,pred_ageclass)) +
    ggplot2::geom_raster(ggplot2::aes(fill = val)) +
    ggplot2::scale_fill_gradient(name = "Proportion",low = "white",high = "black",
                                 limits = c(0,1)) +
    ggplot2::facet_wrap( ~ Prop) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Observed ageclass",
                  y = "Predicted ageclass")))

  return(p1)

  # skipping the rest for now

  # table(ageclass = compare_obs_pred$ageclass,pred_ageclass = compare_obs_pred$pred_ageclass)
  #
  # prop.table(table(compare_obs_pred$ageclass,compare_obs_pred$pred_ageclass),
  #            margin = 2)
  #
  #
  # compare_obs_pred$PARAM <-compare_obs_pred[,param]
  #
  #
  # compare_obs_pred %>%
  #
  #
  # compare_obs_pred %>%
  #   dplyr::mutate(tmp_id = 1:dplyr::n()) %>%
  #   tidyr::pivot_longer(cols = dplyr::matches(age_cols),names_to = "AgeProb",
  #                       values_to = "Prob") %>%
  #   group_by(obs_time,ageclass,tmp_id) %>%
  #   summarize(Age1 = mean(ageclass == 1),
  #             Age2 = mean(ageclass == 2),
  #             Age3 = mean(ageclass == 3)) %>%
  #   # ungroup() %>%
  #   select(Bin,obs_time,Age1,Age2,Age3) %>%
  #   right_join(age_join_by_df[,c("PITCode","Bin", "ageclass",  "obs_time" )])
  #
  #
  # obs_ageclass_df <- Prop_by_bin %>%
  #   ungroup() %>%
  #   mutate(MyALK_age = apply(Prop_by_bin[,3:5],1,function(x) which(x == max(x))[1])) %>%
  #   mutate(Est_ageclass = ifelse(!is.na(ageclass),ageclass, MyALK_age)) %>%
  #   rename(id = PITCode)




}

# to fix no visible binding note
pred_ageclass <- count <- count_ageclass <- count_predageclass <-
  val <- NULL
