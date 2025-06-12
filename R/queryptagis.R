
#'
query_ptagis_mmr_metadata <- function() {
  url_ptagis_mrr <- "https://api.ptagis.org/sites/mrr"

  req_mrr <- httr2::request(url_ptagis_mrr) %>%
    # httr2::req_url_query(siteCode = mmr_site) %>%
    httr2::req_progress()
    # httr2::req_perform()
    # httr2::resp_content_type()

  resp_mrr <- req_mrr %>%
    httr2::req_perform()

  resp_mrr_list <- resp_mrr %>%
    httr2::resp_body_json()

  df_list <- lapply(resp_mrr_list,as.data.frame)

  cnames <- unique(unlist(lapply(df_list,colnames)))

  template <- rep(NA,length(cnames)); names(template) = cnames

  df_list2 <- lapply(df_list,function(x) {template[names(x)] <- x; return(template)})

  parsed_df <- do.call(rbind.data.frame,df_list2)

  return(parsed_df)

}


query_ptagis_mmr_returnsite <- function(mmr_site) {
  # need to rewrite if want to add multiple mmr sites
    url_ptagis_mrr <- "https://api.ptagis.org/sites/mrr"

  stopifnot(length(mmr_site) == 1)

  url_ptagis_mrr <- paste0(url_ptagis_mrr,'/',mmr_site)

  req_mrr <- httr2::request(url_ptagis_mrr) %>%
    # httr2::req_url_query(siteCode = mmr_site) %>%
    httr2::req_progress()
  # httr2::req_perform()
  # httr2::resp_content_type()

  resp_mrr <- req_mrr %>%
    httr2::req_perform()

  resp_mrr_list <- resp_mrr %>%
    httr2::resp_body_json()

  as.data.frame(resp_mrr_list)


}

query_ptagis_obs <- function(pit_codes) {
  #DOESN'T WORK
  pit_codes = "3DD.0077B03C81"
  # need to rewrite if want to add multiple mmr sites
  url_ptagis_events <- "https://api.ptagis.org/data/events"
  # years = 2020
  # stopifnot(length(mmr_site) == 1)

  # url_ptagis_mrr <- paste0(url_ptagis_mrr,'/',mmr_site)

  req_events <- httr2::request(url_ptagis_events) %>%
    httr2::req_url_query(tagCode = pit_codes,
                         apiKey = "",
                         .multi = "explode") %>%
    httr2::req_progress()
  # httr2::req_perform()
  # httr2::resp_content_type()

  resp_events <- req_events %>%
    httr2::req_perform()

  resp_encoding(resp_mrr)
  resp_body_string(resp_mrr) %>% tail()
  resp_mrr_list <- resp_mrr %>%
    httr2::resp_body_string() %>%
    readr::read_csv() %>%
    dplyr::filter(!is.na(NFish)) %>%
    dplyr::mutate(Year = as.integer(Year),
                  Length = as.numeric(Length))




  return(resp_mrr_list)
}



query_dart_release_recap <- function(mmr_sites,recap_sites, years) {
  url_dart_mrr_recap <- "https://www.cbr.washington.edu/dart/cs/php/rpt/mg.php"


  recap_sites <- c("LGR juveniles" = "LWG",
                 "LGR adult ladder" = "GRA","GOJ",
                 "Lower Monumental juveniles" = "LMN",
                 "Little Goose juveniles" = "LGS",
                 "Ice Harbor comb juveniles and adults" = "IHA",
                 "McNary juveniles" = "MCN",
                 "McNary adults" = "MCA",
                 "John Day juveniles" = "JDA",
                 "BON juveniles" = "BON",
                 "BON adult fishway" = "B2A",
                 # "Columbia estuary rkm 62" = "PD5",
                 # "Columbia estuary rkm 68" = "PD6",
                 # "Estuary rkm 70" = "PD7",
                 # "Columbia estuary rkm 82" = "PD8",
                 # "Columbia estuary rkm 75 OR" = "PDO",
                 # "Columbia estuary rkm 75 WA" = "PDW",
                 "Estuary towed array" = "TWX"
  )

  req_mrr <- httr2::request(url_dart_mmr) %>%
    httr2::req_url_query(sc = 1,
                         mgconfig="pit_rel",
                         outputFormat = "csvSingle",
                         year = as.character(years),
                         rel_site = mmr_sites,
                         loc = recap_sites,
                         species = "3",
                         run = "Null",
                         rear_type = "Null",
                         stage = "J",
                         span = "yes", # yes
                         startdate = "1/1",
                         enddate="12/31",
                         syear = as.character(min(years)),
                         eyear = as.character(max(years)),
                         summary = "no",
                         .multi = "explode") %>%
    httr2::req_progress()
  # httr2::req_perform()
  # httr2::resp_content_type()

  resp_mrr <- req_mrr %>%
    httr2::req_perform()

  resp_encoding(resp_mrr)
  resp_body_string(resp_mrr) %>% tail()
  resp_mrr_list <- resp_mrr %>%
    httr2::resp_body_string() %>%
    readr::read_csv() %>%
    dplyr::filter(!is.na(NFish)) %>%
    dplyr::mutate(Year = as.integer(Year),
                  Length = as.numeric(Length))




  return(resp_mrr_list)
}


#
query_dart_release <- function(mmr_sites, years) {
  # need to rewrite if want to add multiple mmr sites
  url_dart_mmr <- "https://www.cbr.washington.edu/dart/cs/php/rpt/pit_rel_de.php"
  # years = 2020
  # stopifnot(length(mmr_site) == 1)

  # url_ptagis_mrr <- paste0(url_ptagis_mrr,'/',mmr_site)

  req_mrr <- httr2::request(url_dart_mmr) %>%
    httr2::req_url_query(sc = 1,
                         outputFormat = "csv",
                         year = as.character(years),
                         rel_site = mmr_sites,
                         species = "3",
                         run = "Null",
                         rear_type = "Null",
                         stage = "J",
                         span = "yes", # yes
                         startdate = "1/1",
                         enddate="12/31",
                         syear = as.character(min(years)),
                         eyear = as.character(max(years)),
                         summary = "no",
                         .multi = "explode") %>%
    httr2::req_progress()
  # httr2::req_perform()
  # httr2::resp_content_type()

  resp_mrr <- req_mrr %>%
    httr2::req_perform()

  resp_encoding(resp_mrr)
  resp_body_string(resp_mrr) %>% tail()
  resp_mrr_list <- resp_mrr %>%
    httr2::resp_body_string() %>%
    readr::read_csv() %>%
    dplyr::filter(!is.na(NFish)) %>%
    dplyr::mutate(Year = as.integer(Year),
                  Length = as.numeric(Length))




  return(resp_mrr_list)
}



#
query_dart_obs <- function(obs_sites, years) {
  # need to rewrite if want to add multiple mmr sites
  url_dart_obs <- "https://www.cbr.washington.edu/dart/cs/php/rpt/pitall_obs_de.php"
  # years = 2020
  # stopifnot(length(mmr_site) == 1)

  # url_ptagis_mrr <- paste0(url_ptagis_mrr,'/',mmr_site)

  req_obs <- httr2::request(url_dart_obs) %>%
    httr2::req_url_query(sc = 1,
                         queryName = "pit_obs_de",
                         outputFormat = "csv",
                         year = as.character(years),
                         proj = obs_sites[1],
                         species = "3",
                         run = "Null",
                         rear_type = "Null",
                         stage = "J",
                         span = "yes", # yes
                         startdate = "1/1",
                         enddate="6/30",
                         syear = as.character(min(years)),
                         eyear = as.character(max(years)),
                         summary = "no",
                         .multi = "explode") %>%
    httr2::req_progress()
  # httr2::req_perform()
  # httr2::resp_content_type()

  resp_obs <- req_obs %>%
    httr2::req_perform()

  # resp_encoding(resp_mrr)
  # resp_body_string(resp_mrr) %>% tail()
  resp_obs_list <- resp_obs %>%
    httr2::resp_body_string() %>%
    readr::read_csv() %>%
    # dplyr::filter(!is.na(NFish)) %>%
    dplyr::mutate(Year = as.integer(Year),
                  Length = as.numeric(Length))


  resp_obs_list

  return(resp_obs_list)
}


query_ptagis_snake_river <- function(mmr_sites, obs_sites, years) {
  mmr_sites <- c("BBCTRP","BIGBEC")
  obs_sites <- c("LGR juveniles" = "LWG","LGR adult ladder" = "GRA","GOJ",
                 "Lower Monumental juveniles" = "LMN",
                 "Little Goose juveniles" = "LGS",
                 "Ice Harbor comb juveniles and adults" = "IHA",
                 "McNary juveniles" = "MCN",
                 "McNary adults" = "MCA",
                 "John Day juveniles" = "JDA",
                 "BON juveniles" = "BON",
                 "BON adult fishway" = "B2A",
                 # "Columbia estuary rkm 62" = "PD5",
                 # "Columbia estuary rkm 68" = "PD6",
                 # "Estuary rkm 70" = "PD7",
                 # "Columbia estuary rkm 82" = "PD8",
                 # "Columbia estuary rkm 75 OR" = "PDO",
                 # "Columbia estuary rkm 75 WA" = "PDW",
                 "Estuary towed array" = "TWX"
                 )

  dat_rel <- query_dart_release(mmr_sites = mmr_sites,years = years)

  dat_rec <- query_dart_obs(obs_sites = obs_sites,years = years)
}
