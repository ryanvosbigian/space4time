



process_site_config <- function(DART_config,configdate = Sys.Date()) {
  lines <- tryCatch(readLines(DART_config),
                    error = function(err) stop(paste0("DART site configuration file failed: ",DART_config)))

  message(paste0("Obtained the PTAGIS site configuration formatted by DART from:\n",DART_config))

  parsed_df <- data.frame(sitecode = as.character(NA),
                          arrayname = as.character(NA),
                          arraycodes = as.character(NA),
                          date_start = as.Date(NA),
                          date_end = as.Date(NA),
                          releasecode1 = as.character(NA),
                          releasecode2 = as.character(NA),
                          exception = as.logical(NA),
                          stringsAsFactors = FALSE)[0,]

  # Find the start and end of the desired section
  suppressWarnings(start_line <- which(grepl("# Detector Configuration", lines)) + 2)
  suppressWarnings(end_line <- which(grepl("# Interrogation Site Configuration", lines)) - 2)

  detector_config_lines <- lines[start_line:end_line]

  code_start <- grep("code:",detector_config_lines)
  code_end <- c(code_start[-1]-1,length(detector_config_lines))


  for (c in 1:length(code_start)) {



    tmp_code_lines <- detector_config_lines[code_start[c]:code_end[c]]
    tmp_openparanth <- grep("[{]$",tmp_code_lines)
    tmp_cloparanth <- grep("[}]$",tmp_code_lines)
    tmp_range <- grep("range:",tmp_code_lines)



    # or length tmp_range
    for (r in 1:(length(tmp_range))) {
      # Present <- Sys.Date()


      # tmp_code_lines[tmp_range[p-1]]
      tmpdate_range <- gsub("^.*range: ","",tmp_code_lines[tmp_range[r]])

      tmp_startdate <- as.Date(stringr::str_split_i(tmpdate_range," ",1),format = "%d-%b-%y")

      #
      tmp_enddate <- ifelse(stringr::str_split_i(tmpdate_range," ",2) == "Present",as.Date(configdate),as.Date(stringr::str_split_i(tmpdate_range," ",2),format = "%d-%b-%y"))

      next_closed_paran <- min(tmp_cloparanth[which(tmp_cloparanth-tmp_range[r] > 0)])

      tmp_array_lines <- tmp_code_lines[seq(tmp_range[r] + 2,next_closed_paran-1)]

      arraynames <- stringr::str_split_i(tmp_array_lines,pattern = " : ",i = 3)
      arraycodes <- stringr::str_split_i(tmp_array_lines,pattern = " : ",i = 4)

      releasecode1 <- stringr::str_split_i(tmp_array_lines,pattern = ": ",i = 1)
      releasecode1 <- gsub(" ","",releasecode1)

      releasecode2 <- stringr::str_split_i(tmp_array_lines,pattern = ":",i = 2)

      releasecode2 <- gsub(" ","",releasecode2)

      new_row <- data.frame(sitecode = gsub("code: ","",tmp_code_lines[1]),
                            arrayname = arraynames,
                            arraycodes = arraycodes,
                            date_start = tmp_startdate,
                            date_end =as.Date(tmp_enddate),
                            releasecode1 = releasecode1,
                            releasecode2 = releasecode2,
                            exception = FALSE,
                            stringsAsFactors = FALSE)

      parsed_df <- rbind(parsed_df, new_row)
    }

    # deal with exceptions:
    if (sum(grepl("exception",tmp_code_lines)) > 0) {
      # stop()
      tmp_exception_lines <- tmp_code_lines[grepl("exception",tmp_code_lines)]

      format_tmp_exception_lines <- gsub("^.*exception: ","",tmp_exception_lines)

      tmp_date_start <- stringr::str_split_i(format_tmp_exception_lines," ",1)
      tmp_date_end <- stringr::str_split_i(format_tmp_exception_lines," ",2)


      tmp_date_start <- as.Date(tmp_date_start,format = "%d-%b-%y")
      tmp_date_end <- as.Date(ifelse(tmp_date_end == "Present",as.Date(configdate),as.Date(tmp_date_end,format = "%d-%b-%y")))

      format2_tmp_exc_lines <- gsub(".*[{] | [}]","",format_tmp_exception_lines)

      arraynames <- stringr::str_split_i(format2_tmp_exc_lines,pattern = " : ",i = 3)
      arraycodes <- stringr::str_split_i(format2_tmp_exc_lines,pattern = " : ",i = 4)

      releasecode1 <- stringr::str_split_i(format2_tmp_exc_lines,pattern = ": ",i = 1)

      releasecode2 <- stringr::str_split_i(format2_tmp_exc_lines,pattern = ":",i = 2)

      releasecode2 <- gsub(" ","",releasecode2)

      new_row <- data.frame(sitecode = gsub("code: ","",tmp_code_lines[1]),
                            arrayname = arraynames,
                            arraycodes = arraycodes,
                            date_start = tmp_date_start,
                            date_end =tmp_date_end,
                            releasecode1 = releasecode1,
                            releasecode2 = releasecode2,
                            exception = TRUE,
                            stringsAsFactors = FALSE)

      parsed_df <- rbind(parsed_df, new_row)


    }


  }

  parsed_df$date_start <- as.Date(parsed_df$date_start)
  parsed_df$date_end <- as.Date(parsed_df$date_end)


  return(parsed_df)
}





# read_DART_file <- function(filepath, ...) {
#   dots <- list(...)
#   if (length(dots)  == 0) {
#     files = list(filepath)
#   } else {
#     # stopifnot(is.character(...))
#     files = list(filepath,...)
#   }
#
#
#   # find first row:
#   skip_df = read.csv(files[[1]])
#
#   skip = grep("[#]RelGrpStartDate",skip_df[,1])
#
#   capture_data <- readr::read_csv(files[[1]],skip = skip) %>%
#     as.data.frame() %>%
#     dplyr::filter(!is.na(TagID))
#
#
#
#   if (length(files) > 1) {
#     for (i in 2:length(files)) {
#       # find first row:
#       skip_df = read.csv(files[[1]])
#
#       skip = grep("[#]RelGrpStartDate",skip_df[,1])
#
#       tmp <- readr::read_csv(filepath,comment = "##") %>%
#         as.data.frame() %>%
#         dplyr::filter(!is.na(TagID))
#
#       capture_data <- rbind(capture_data,tmp)
#
#     }
#   }
#
#   capture_data$time <- tryCatch(lubridate::as_datetime(capture_data$RelVTime),
#                                 warning = function(e){
#                                   tryCatch(lubridate::mdy(capture_data$RelVTime),
#                                            warning = function(e2){
#                                              lubridate::as_date(
#                                                stringr::str_split_i(capture_data$RelVTime,
#                                                                     " ",
#                                                                     1),
#                                                format = "%m/%d/%Y")
#                                            })
#
#
#                                 })
#
#
#   capture_data$RecTime <- tryCatch(lubridate::as_datetime(capture_data$ObsDateLast),
#                                 warning = function(e){
#                                   tryCatch(lubridate::mdy(capture_data$ObsDateLast),
#                                            warning = function(e2){
#                                              lubridate::as_date(
#                                                stringr::str_split_i(capture_data$ObsDateLast,
#                                                                     " ",
#                                                                     1),
#                                                format = "%m/%d/%Y")
#                                            })
#
#
#                                 })
#
#
#   initialreleases <- capture_data %>%
#     dplyr::mutate(#time = lubridate::as_datetime(RelVTime),
#            time = lubridate::year(time),
#            # RelSite = ifelse(RelSite == "BBCTRP","BIGBEC",RelSite),
#            removed = FALSE) %>%
#     dplyr::select(id = TagID,site = RelSite,time = time, removed = removed)  %>%
#     dplyr::distinct(id,site,time,removed)
#
#   recaps <- capture_data %>%
#     dplyr::filter(!is.na(ObsSite)) %>%
#     # dplyr::mutate(ObsSite = case_when(ObsSite %in% c("B2J","BCC") ~ "BON",
#     #                            TRUE ~ ObsSite)) %>%
#     dplyr::filter(TagID %in% initialreleases$id) %>%
#     dplyr::mutate(removed = FALSE) %>%
#     dplyr::mutate(time = RecTime,#lubridate::as_datetime(ObsDateLast),
#            time = lubridate::year(time)) %>%
#     dplyr::select(id = TagID,site = ObsSite,time = time, removed = removed)
#
#   combined_capture_data <- rbind(initialreleases,recaps) %>%
#     dplyr::arrange(id, site, time)
#
#   return(combined_capture_data)
# }


# fix no visible binding note
ObsDateLast <- ObsSite <- RelSite <- RelVTime <- TagID <- NULL


# source("C:/Users/rvosbigian/OneDrive - University of Idaho/Post_masters/cjs_s4t/process_site_info_01_06_2025.R")
identify_barged_fish <- function(capture_data,parsed_df) {

  N_nositeinfo <- 0

  # investigate removed individuals
  for (i in 1:nrow(capture_data)) {
    if (is.na(capture_data$ObsSite[i])) next()

    tmp_obssite <- capture_data$ObsSite[i]
    tmp_obsmonitor <- capture_data$ObsMonitor[i]
    tmp_time <- capture_data$time[i]

    if (!(tmp_obssite %in% c("GOJ","GRJ","LMJ"))) next()

    tmp_siteinfo <- parsed_df %>%
      dplyr::filter(sitecode == tmp_obssite,
             grepl(tmp_obsmonitor,arrayname),
             date_start < tmp_time,
             date_end > tmp_time)

    if (nrow(tmp_siteinfo) == 0) {
      capture_data$removed[i] <- FALSE
      N_nositeinfo <- N_nositeinfo + 1
      # stop("1")
    }

    if (length(unique(tmp_siteinfo$releasecode1)) == 1) {
      capture_data$removed[i] <- ifelse(tmp_siteinfo$releasecode1[1] == 'T',TRUE,FALSE)
    } else {

      if (nrow(tmp_siteinfo) > 1) {
        tmp_siteinfo <- tmp_siteinfo %>%
          dplyr::filter(exception == TRUE)

        if (nrow(tmp_siteinfo) == 1) {
          capture_data$removed[i] <- ifelse(tmp_siteinfo$releasecode1 == 'T',TRUE,FALSE)
        } else {
          tmp_releasecode <-unique(tmp_siteinfo$releasecode1)
          if (length(tmp_releasecode) == 1) {
            capture_data$removed[i] <- ifelse(tmp_releasecode == 'T',TRUE,FALSE)
          } else {
            stop("2")
          }

        }


      }

      # capture_data$removed[i] <- ifelse(tmp_siteinfo$releasecode1 == 'T',TRUE,FALSE)


    }

  }

  if (N_nositeinfo > 0) {
    warning(paste0("Number of observations without info. Assuming not transported. N = ",N_nositeinfo))
  }


  return(capture_data$removed)


} # end function




#' Read DART Basin TribPit observation files and combine with auxiliary age data
#'
#' @description
#' Read DART Basin TribPit observation files and combine with auxiliary age data
#'
#'
#'
#'
#' @param filepath the path to the DART Basin TribPit observation file
#' @param aux_age_df A data frame object in the same format as aux_age_df
#'     (see documentation for `s4t_ch`)
#' @param DART_config the path to a site configuration file `"site_config.txt"` generated by
#'     Columbia Basin Research. If blank, defaults to `www.cbr.washington.edu/downloads/paramest/sites_config.txt`
#' @param ... any additional file paths of DART files to include
#'
#' @returns A list with two elements. The first is the formatted capture history data (`ch_df`)
#'     and the second is the age auxiliary data (`aux_age_df`).
#'
#'
#' @details
#' Currently, only Dams in the Snake River are checked for transported fish.
#'
#'
#' @export
#'
#'
read_DART_file <- function(filepath,aux_age_df,DART_config = "https://www.cbr.washington.edu/downloads/paramest/sites_config.txt", ...) {
  dots <- list(...)
  if (length(dots)  == 0) {
    files = list(filepath)
  } else {
    # stopifnot(is.character(...))
    files = list(filepath,...)
  }


  age_df <- aux_age_df

  # find first row:
  skip_df = utils::read.csv(files[[1]])

  skip = grep("[#]RelGrpStartDate",skip_df[,1])

  suppressMessages(capture_data <- readr::read_csv(files[[1]],skip = skip) %>%
    as.data.frame() %>%
    dplyr::filter(!is.na(TagID))
    )


  if (length(files) > 1) {
    for (i in 2:length(files)) {
      # find first row:
      skip_df = utils::read.csv(files[[i]])

      skip = grep("[#]RelGrpStartDate",skip_df[,1])

      # skip_df = readLines(files[[i]])
      #
      # skip = grep("[#]RelGrpStartDate",skip_df)

      suppressWarnings(suppressMessages(tmp <- readr::read_csv(files[[i]],comment = "##",skip = skip) %>%
        as.data.frame() %>%
        dplyr::filter(!is.na(TagID))))

      capture_data <- rbind(capture_data,tmp)

    }
  }




  capture_data$time <- tryCatch(lubridate::as_datetime(capture_data$RelVTime),
                                warning = function(e){
                                  tryCatch(lubridate::mdy(capture_data$RelVTime),
                                           warning = function(e2){
                                             lubridate::as_date(
                                               stringr::str_split_i(capture_data$RelVTime,
                                                                    " ",
                                                                    1),
                                               format = "%m/%d/%Y")
                                           })


                                })





  capture_data$RecTime <- tryCatch(lubridate::as_datetime(capture_data$ObsDateLast),
                                   warning = function(e){
                                     tryCatch(lubridate::mdy(capture_data$ObsDateLast),
                                              warning = function(e2){
                                                lubridate::as_date(
                                                  stringr::str_split_i(capture_data$ObsDateLast,
                                                                       " ",
                                                                       1),
                                                  format = "%m/%d/%Y")
                                              })


                                   })

  # if (is.null(DART_config)) {
  #   DART_config <- "https://www.cbr.washington.edu/downloads/paramest/sites_config.txt"
  # }

  parsed_df <- process_site_config(DART_config)

  capture_data$removed <- FALSE
  capture_data$removed <- (identify_barged_fish(capture_data,parsed_df))

  initialreleases <- capture_data %>%
    dplyr::mutate(#time = lubridate::as_datetime(RelVTime),
      time = lubridate::year(time),
      # RelSite = ifelse(RelSite == "BBCTRP","BIGBEC",RelSite),
      removed = FALSE) %>%
    dplyr::select(id = TagID,site = RelSite,time = time, removed = removed)  %>%
    dplyr::distinct(id,site,time,removed)

  recaps <- capture_data %>%
    dplyr::filter(!is.na(ObsSite)) %>%
    # dplyr::mutate(ObsSite = case_when(ObsSite %in% c("B2J","BCC") ~ "BON",
    #                            TRUE ~ ObsSite)) %>%
    dplyr::filter(TagID %in% initialreleases$id) %>%
    # dplyr::mutate(removed = FALSE) %>%
    dplyr::mutate(time = RecTime,#lubridate::as_datetime(ObsDateLast),
                  time = lubridate::year(time)) %>%
    dplyr::select(id = TagID,site = ObsSite,time = time, removed = removed)

  ### keep the first observation (typically kelts make up the later observations)
  recaps <- recaps %>%
    dplyr::group_by(id, site) %>%
    dplyr::arrange(id, site, time) %>%
    dplyr::mutate(N = 1:dplyr::n()) %>%
    dplyr::filter(N == 1) %>%
    dplyr::select(id, site, time, removed)


  combined_capture_data <- rbind(initialreleases,recaps) %>%
    dplyr::arrange(id, site, time)



  # adding rows to age_df if missing
  missing_obs_aux <- setdiff(capture_data$TagID,age_df$id)
  if (length(missing_obs_aux) > 0) {
    needed_cols <- c("obs_time","obs_site","id","ageclass","FL")
    check_cols <- setdiff(needed_cols,colnames(age_df))

    if (length(check_cols) > 0) {
      message(paste0("These columns are missing from age_df:",
                     paste0(check_cols,collapse = ", ")))
      stop("aux_age_df does not contain necessary columns")
    }

    needed_cap_cols <- c("RelSite","TagID","Lgth")
    check_cap_cols <- setdiff(needed_cap_cols,colnames(capture_data))
    if (length(check_cap_cols) > 0) {
      message(paste0("These columns are missing from DART file:",
                     paste0(check_cols,collapse = ", ")))
      stop(paste0("DART fils do not contain necessary columns"))
    }

  row_na_df <- age_df[1,]
  row_na_df[,1:ncol(row_na_df)] <- NA
  if (length(missing_obs_aux) > 1) {
    row_na_df <- row_na_df %>% tibble::add_row(id = rep(NA,length(missing_obs_aux)-1))
  }


  na_capture_data <- capture_data %>%
    dplyr::filter(TagID %in% missing_obs_aux) %>%
    dplyr::distinct(TagID,RelSite,time,.keep_all = TRUE)


  row_na_df$obs_site <- na_capture_data$RelSite

  row_na_df$id <- na_capture_data$TagID

  row_na_df$obs_time <- lubridate::year(na_capture_data$time)

  row_na_df$FL <- na_capture_data$Lgth

  message(paste0("Appending N = ",length(missing_obs_aux)," observations to age_df. \n    Using RelSite as obs_site, year of release as obs_time, and Lgth as FL"))

  message(paste0("Note: the following columns were populated with NAs: \n    ",
                 paste0(setdiff(colnames(row_na_df),c("obs_site",
                                                      "id",
                                                      "obs_time",
                                                      "FL")),
                        collapse = ", "
                        )
                 )
          )

  age_df <- rbind(age_df,row_na_df)
  }

  # capture_data


  return(list(ch_df = combined_capture_data,
              aux_age_df = age_df))
}


#' Remove observations of kelts
#'
#' @description
#' Removes any observations following an observation at specified sites. Note that
#'     the observations at the specified sites are retained.
#'
#'
#' @export
#'
#' @param ch_df a capture history data frame (see documentation for s4t_ch)
#' @param kelt_obssite a vector of sites that identify kelts (i.e. adult fish ladders)
#'
#' @returns a capture history data frame without observations following
#'     observations at the specified sites.
#'
#' @export
#'
#' @examples
#'
remove_kelt_obs <- function(ch_df, kelt_obssite) {
  kelt_sites <- intersect(as.character(kelt_obssite),unique(ch_df$site))

  message(paste0("Dropping observations after observations at the following sites: ",paste0(kelt_sites,collapse = ", ")))


  kelt_obs <- ch_df %>%
    dplyr::filter(site %in% kelt_sites) %>%
    dplyr::group_by(id) %>%
    dplyr::arrange(time) %>%
    dplyr::summarize(kelt_time = dplyr::first(time))

  remove_kelts_ch_df <- ch_df %>%
    dplyr::left_join(kelt_obs, by = "id") %>%
    dplyr::group_by(id) %>%
    dplyr::filter(is.na(kelt_time) | time <= kelt_time) %>%
    dplyr::ungroup() %>%
    dplyr::select(id, site, time, removed)

  return(as.data.frame(remove_kelts_ch_df))

}



estimate_cohort_surv_vc <- function(obj) {
  res <- obj$res
  s4t_ch <- obj$fit$s4t_ch
  obj$fit$indices_theta$g

  n_init_relsite <- s4t_ch$ch_info$n_init_relsite
  n_groups <- length(unique(obj$fit$indices_theta$g)) # may cause issues



  init_relsite_list <- s4t_ch$ch_info$init_relsite_list
  set_min_a <- s4t_ch$ch_info$observed_relative_min_max$set_min_a
  set_max_a <- s4t_ch$ch_info$observed_relative_min_max$set_max_a


  max_s_rel <- s4t_ch$ch_info$max_s_rel
  max_t_recap <- s4t_ch$ch_info$max_t_recap
  n_sites <- s4t_ch$ch_info$n_sites
  recap_sites <- s4t_ch$ch_info$recap_sites
  last_sites <- s4t_ch$ch_info$last_sites
  recap_sites_not_last <- s4t_ch$ch_info$recap_sites_not_last
  not_last_sites <- c(1:n_sites)[!(1:n_sites %in% s4t_ch$ch_info$last_sites)]
  sites_config <- s4t_ch$s4t_config$sites_config
  holdover_config <- s4t_ch$s4t_config$holdover_config



  mod_mat_theta <- obj$fit$mod_mat_theta
  indices_theta <- obj$fit$indices_theta

  min_ageclass_mat <- s4t_ch$ch_info$min_ageclass_mat
  max_ageclass_mat <- s4t_ch$ch_info$max_ageclass_mat

  ja <- numDeriv::jacobian(function(x, ...) stats::plogis(return_cohort_surv(par = x, ...)),
                           x = res$par,
                           n_init_relsite = n_init_relsite,
                           init_relsite_list = init_relsite_list,
                           n_groups = n_groups,
                           set_min_a = set_min_a,
                           set_max_a = set_max_a,
                           max_s_rel = max_s_rel,
                           max_t_recap = max_t_recap,
                           n_sites = n_sites,
                           recap_sites = recap_sites,
                           last_sites = last_sites,
                           recap_sites_not_last = recap_sites_not_last,
                           not_last_sites = not_last_sites,
                           sites_config = sites_config,
                           holdover_config = holdover_config,
                           mod_mat_theta = mod_mat_theta,
                           indices_theta = indices_theta,
                           min_ageclass_mat = min_ageclass_mat,
                           max_ageclass_mat = max_ageclass_mat)

  vc <- solve(res$hessian)

  der_vc <- ja %*% vc %*% t(ja)

  return(der_vc)

}




#' Calculate abundance estimates
#'
#' @description
#' Calculate abundance estimates of cohorts
#'
#'
#' @export
#'
#' @param obj a capture history data frame (see documentation for s4t_ch)
#' @param abund a vector of sites that identify kelts (i.e. adult fish ladders)
#' @param type the type of summarization for the abundance estimates. BroodYear
#'     summarizes by broodyear (s - a1), ReleaseYear summarizes by time (s),
#'     and None does not summarize.
#'
#' @returns a capture history data frame without observations following
#'     observations at the specified sites.
#'
#' @export
#'
#' @examples
#'
abundance_estimates <- function(obj, abund, type = c("BroodYear","ReleaseYear","None")) {

  abund <- as.data.frame(abund)

  tmp_mis <- setdiff(c("a1","s","j","abundance","abundance_se"),colnames(abund))
  if ("abundance_se" %in% tmp_mis) {
    abund$abundance_se <- 0
  }

  tmp_mis <- setdiff(c("a1","s","j","abundance","abundance_se"),colnames(abund))

  if (length(tmp_mis) > 0) {
    stop(paste0("Missing columns from `abund`: ",paste0(tmp_mis,collapse = ", ")))
  }

  tmp_add <- setdiff(colnames(abund),c("a1","s","j","r","g","abundance","abundance_se"))

  if (length(tmp_add) > 0) {
    message(paste0("Only the following columns can be included in `abund`: ",
                   paste0(c("a1","s","j","r","g","abundance","abundance_se"),collapse = ", ")))
    stop(paste0("Remove columns from `abund`: ",paste0(tmp_add,collapse = ", ")))
  }

  abund$j <- as.character(abund$j)



  cohort_transitions <- obj$cohort_transitions
  suppressMessages(transition_df <- dplyr::left_join(cohort_transitions,abund))

  if (is(obj,"s4t_cjs_ml")) {



    # cohort_vc <- obj$fit$cohort_surv_vc
    cohort_vc <- estimate_cohort_surv_vc(obj = obj)



  } else if (is(obj,"s4t_cjs_rstan")) {
    colnames(transition_df)[colnames(transition_df) == "mean"] <- "estimate"
    colnames(transition_df)[colnames(transition_df) == "mean"] <- "estimate_se"

    tmp_cohort_trans_parnames <- paste0("cohort_surv[",1:obj$res@sim$dims_oi$cohort_surv,"]")
    # rstan::as.matrix.stanfit(obj$res,pars = tmp_cohort_trans_parnames)
    tmp_cohort <- rstan::extract(obj$res,pars = tmp_cohort_trans_parnames)
    cohort_vc <- stats::cov(do.call(cbind,tmp_cohort))

    # stats::cov()
  }



  # transition_df$abundance_se <- runif(nrow(transition_df),20,50)

  abund_vc <- diag(transition_df$abundance_se^2)



  outer_N <- transition_df$abundance %*% t(transition_df$abundance)
  outer_Theta <- transition_df$estimate %*% t(transition_df$estimate)

  cohort_abund_vc <- (abund_vc * cohort_vc) +
    (abund_vc * outer_Theta) +
    (outer_N * cohort_vc)

  # for approximation, just set NaNs to zero:
  if (sum(is.nan(cohort_abund_vc) > 0)) {
    message("Setting Nan values in vcov matrix to 0, results are approximate.")
    cohort_abund_vc <- ifelse(is.nan(cohort_abund_vc),0,cohort_abund_vc)
  }
  if (sum(is.na(cohort_abund_vc) > 0)) {
    message("Setting Nan values in vcov matrix to 0, results are approximate.")
    cohort_abund_vc <- ifelse(is.na(cohort_abund_vc),0,cohort_abund_vc)
  }


  calc_se <- function(Index,cohort_abund_vc) {

    a = matrix(rep(0,nrow(cohort_abund_vc)),nrow=1)
    a[Index] <- 1
    as.numeric(sqrt(a %*% cohort_abund_vc %*% t(a)))
  }

  if (type == "BroodYear") {





    transition_df %>%
      dplyr::mutate(Index = 1:dplyr::n(),
                    broodyear = s - a1) %>%
      dplyr::mutate(cohort_abund = estimate * abundance,
             cohort_abund_se = sqrt((abundance_se^2 + abundance^2) * (diag(cohort_vc) + estimate^2) - ((abundance^2) * (estimate^2)))) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(broodyear,r,g,j) %>% #View()
      dplyr::summarize(abundance_broodyear = sum(cohort_abund),
                       # paste0(Index,collapse = " "),
                abundance_broodyear_se = calc_se(Index,cohort_abund_vc)
                )


  } else if (type == "None") {
    transition_df %>%
      dplyr::mutate(Index = 1:dplyr::n(),
                    broodyear = s - a1) %>%
      dplyr::mutate(estimate_se = sqrt(diag(cohort_vc)),
                    cohort_abund = estimate * abundance,
                    cohort_abund_se = sqrt((abundance_se^2 + abundance^2) * (diag(cohort_vc) + estimate^2) - ((abundance^2) * (estimate^2))))
  } else if (type == "ReleaseYear") {
    transition_df %>%
      dplyr::mutate(Index = 1:dplyr::n(),
                    broodyear = s - a1) %>%
      dplyr::mutate(cohort_abund = estimate * abundance,
                    cohort_abund_se = sqrt((abundance_se^2 + abundance^2) * (diag(cohort_vc) + estimate^2) - ((abundance^2) * (estimate^2)))) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(s,j,r,g) %>%
      dplyr::summarize(abundance_releaseyear = sum(cohort_abund),
                       abundance_releaseyear_se = calc_se(Index,cohort_abund_vc))

  }

}


arrayname <- sitecode <- date_start <- date_end <- exception <- RecTime <- kelt_time <-
  abundance <- abundance_se <- broodyear <- cohort_abund <- Index <-  NULL
