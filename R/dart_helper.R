


read_DART_file <- function(filepath, ...) {
  dots <- list(...)
  if (length(dots)  == 0) {
    files = list(filepath)
  } else {
    # stopifnot(is.character(...))
    files = list(filepath,...)
  }


  capture_data <- readr::read_csv(files[[1]],comment = "##") %>%
    as.data.frame()

  if (length(files) > 1) {
    for (i in 1:length(files)) {
      tmp <- readr::read_csv(filepath,comment = "##") %>%
        as.data.frame()
      capture_data <- rbind(capture_data,tmp)

    }
  }

  initialreleases <- capture_data %>%
    dplyr::mutate(time = lubridate::as_datetime(RelVTime),
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
    dplyr::mutate(removed = FALSE) %>%
    dplyr::mutate(time = lubridate::as_datetime(ObsDateLast),
           time = lubridate::year(time)) %>%
    dplyr::select(id = TagID,site = ObsSite,time = time, removed = removed)

  combined_capture_data <- rbind(initialreleases,recaps) %>%
    dplyr::arrange(id, site, time)

  return(combined_capture_data)
}


# fix no visible binding note
ObsDateLast <- ObsSite <- RelSite <- RelVTime <- TagID <- NULL
