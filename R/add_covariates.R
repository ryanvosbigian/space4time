

#' Add covariates to a space4time capture history object
#'
#' @description
#' A short description...
#'
#' @param cov_df a `data.frame` or `list` of `data.frame`'s containing the covariates
#'     for theta and p `a1,a2,j,k,s,t,r,g` indices. See details.
#' @param s4t_ch a `s4t_ch` object.
#' @returns a `s4t_ch` object with covariates added
#' @details
#' Additional details...
#'
#'
#' @export
add_covariates <- function(cov_df,s4t_ch) {
  if (is(cov_df,"data.frame") | is(cov_df,"tbl_df")) {
    cov_df <- list(cov_df)
  }




  tmp_format <- format_s4t_cjs(p_formula = ~1,
                               theta_formula = ~1,
                               ageclass_formula =  ~1,
                               s4t_ch = s4t_ch)


  if (is.null(s4t_ch$cov_df$cov_p)) {
    s4t_ch$cov_df$cov_p <- list()
  }

  if (is.null(s4t_ch$cov_df$cov_theta)) {
    s4t_ch$cov_df$cov_theta <- list()
  }

  indices_theta_original <- tmp_format$indices_theta_original

  indices_p_obs_original <- tmp_format$indices_p_obs_original

  # err <- FALSE

  for (i in 1:length(cov_df)) {

    ## start with p_obs
    # determine if the extent of the levels are enough:
    tmp_df <- cov_df[[i]]
    # start with p_obs
    match_col <- intersect(colnames(tmp_df),colnames(indices_p_obs_original))
    cov_col <- setdiff(colnames(tmp_df),match_col)

    # test just one set of levels:
    for (j in match_col) {
      tmp_df[,j] <- as.factor(tmp_df[,j])
      lvls_cov <- levels(tmp_df[,j])
      lvls_ind <- levels(indices_p_obs_original[,which(colnames(indices_p_obs_original) == j)])
      lvls_diff <- setdiff(lvls_ind,lvls_cov)
      if(length(lvls_diff) > 0) {
        # err <- TRUE
        message(paste0("Following levels are missing from index ",j," when joining ",
                       "\ncovariate(s) ",paste0(cov_col,collapse = ", "),
                       ": \nLevel(s):",paste0(lvls_diff,collapse = ", ")))
      }

    }
    # test coverage of the levels put together:
    tmp_indices_p_obs_original <- dplyr::left_join(indices_p_obs_original,stats::na.fail(tmp_df))
    na_rows <- !stats::complete.cases(tmp_indices_p_obs_original[,cov_col])

    added_col <- NULL
    if (sum(na_rows) > 0) {
      # then some of the levels are incomplete
      tmp_indices_p_obs_original$INDICATORVARIABLE <- ifelse(na_rows,0,1)

      added_col <- paste0("OFF_",cov_col,collapse = "_")


      colnames(tmp_indices_p_obs_original)[ncol(tmp_indices_p_obs_original)] <-
        added_col


      tmp_indices_p_obs_original[,cov_col] <- ifelse(is.na(tmp_indices_p_obs_original[,cov_col] ),0,tmp_indices_p_obs_original[,cov_col] )
      message("Some levels from covariate(s) ",paste0(cov_col,collapse = ", "),
              "\n    were missing. Added a column to cov_p to split observations between \n    observed and unobserved values: ",added_col)
    }

    red_p_obs_df <- tmp_indices_p_obs_original[,c(colnames(tmp_df),added_col)]
    new_cov_p <- dplyr::distinct(red_p_obs_df)

    if (sum(na_rows) > 0) {
      missing_levels <- dplyr::distinct(new_cov_p[new_cov_p[,added_col] == 1,match_col])

      message("Missing levels: ")
      print(missing_levels)
    }

    # add to s4t_ch$cov_df$cov_p
    s4t_ch$cov_df$cov_p[[length(s4t_ch$cov_df$cov_p) + 1]] <- new_cov_p


    ## now do theta

    # determine if the extent of the levels are enough:
    tmp_df <- cov_df[[i]]
    # start with p_obs
    match_col <- intersect(colnames(tmp_df),colnames(indices_theta_original))
    cov_col <- setdiff(colnames(tmp_df),match_col)

    # test just one set of levels:
    for (j in match_col) {
      tmp_df[,j] <- as.factor(tmp_df[,j])
      lvls_cov <- levels(tmp_df[,j])
      lvls_ind <- levels(indices_theta_original[,which(colnames(indices_theta_original) == j)])
      lvls_diff <- setdiff(lvls_ind,lvls_cov)
      if(length(lvls_diff) > 0) {
        # err <- TRUE
        message(paste0("Following levels are missing from index ",j," when joining ",
                       "\ncovariate(s) ",paste0(cov_col,collapse = ", "),
                       ": \nLevel(s):",paste0(lvls_diff,collapse = ", ")))
      }

    }
    # test coverage of the levels put together:
    tmp_indices_theta_original <- dplyr::left_join(indices_theta_original,stats::na.fail(tmp_df))
    na_rows <- !stats::complete.cases(tmp_indices_theta_original[,cov_col])
    added_col <- NULL
    if (sum(na_rows) > 0) {
      # then some of the levels are incomplete
      tmp_indices_theta_original$INDICATORVARIABLE <- ifelse(na_rows,1,0)
      colnames(tmp_indices_theta_original)[ncol(tmp_indices_theta_original)] <-
        paste0("OFF_",cov_col,collapse = "_")

      added_col <- paste0("OFF_",cov_col,collapse = "_")

      tmp_indices_theta_original[,cov_col] <- ifelse(is.na(tmp_indices_theta_original[,cov_col] ),0,tmp_indices_theta_original[,cov_col] )
      message("Some levels from covariate(s) ",paste0(cov_col,collapse = ", "),
              "\n    were missing. Added a column to cov_theta to split observations between \n    observed and unobserved values: ",added_col)
    }

    red_theta_df <- tmp_indices_theta_original[,c(colnames(tmp_df),added_col)]
    new_cov_theta <- dplyr::distinct(red_theta_df)

    if (sum(na_rows) > 0) {
      missing_levels <- dplyr::distinct(new_cov_theta[new_cov_theta[,added_col] == 1,match_col])

      message("Missing levels: ")
      print(missing_levels)
    }

  }

  # add to s4t_ch$cov_df$cov_theta
  s4t_ch$cov_df$cov_theta[[length(s4t_ch$cov_df$cov_theta) + 1]] <- new_cov_theta



  return(s4t_ch)
}


#' Add covariates to a space4time capture history object
#'
#' @description
#' A short description...
#'
#' @param cov_df a `data.frame` or `list` of `data.frame`'s containing the covariates
#'     for theta and p `a1,a2,j,k,s,t,r,g` indices. See details.
#' @param s4t_ch a `s4t_ch` object.
#' @returns a `s4t_ch` object with covariates added
#' @details
#' Additional details...
#'
#'
#' @export
extract_covariates <- function(s4t_ch) {

  tmp_format <- suppressMessages(format_s4t_cjs(p_formula = ~1,
                               theta_formula = ~1,
                               ageclass_formula =  ~1,
                               s4t_ch = s4t_ch))


  if (is.null(s4t_ch$cov_df$cov_p)) {
    s4t_ch$cov_df$cov_p <- list()
  }

  if (is.null(s4t_ch$cov_df$cov_theta)) {
    s4t_ch$cov_df$cov_theta <- list()
  }

  indices_theta_original <- tmp_format$indices_theta_original

  indices_p_obs_original <- tmp_format$indices_p_obs_original

  obj <- list(indices_theta = indices_theta_original,
       indices_p_obs = indices_p_obs_original)

  class(obj) = "cov_s4t_ch"

  return(obj)
}
