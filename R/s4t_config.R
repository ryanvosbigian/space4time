

#' Create space4time site configuration object (`s4t_config`)
#'
#' @description
#' Create the space4time site configuration object (`s4t_config`) using
#'     custom site configuration. For a simple configuration where
#'     there is only one initial release site use `linear_s4t_config()`,
#'     which assumes that the sites are arranged in a simple sequence.
#'     If there are multiple initial release sites that are immediately
#'     connected to a simple sequence then `simplebranch_s4t_config()`
#'     can be used. Otherwise, use this function to create a custom site
#'     configuration. See details or (LINK VIGNETTE) for more detailed
#'     instructions.
#'
#' @param sites_names a character `vector` of the site names
#' @param sites_config a `matrix` that determine how sites are linked together.
#'     Must be in the same order as sites_names. See LINK VIGNETTE and details.
#' @param holdover_config a `matrix` that determine whether individuals can holdover in between
#'      sites. In the same format as sites_configs. See LINK VIGNETTE and details.
#' @param min_a a `vector` of the minimum ageclass individuals in a site can be.
#'     Must be the same length and order as `sites_names`
#' @param max_a a `vector` of the maximum ageclass individuals in a site can be.
#'     Must be the same length and order as `sites_names`
#' @param sites_to_pool a named list of character vectors that contain the names of sites
#'      to be pooled together and treated as one site. See examples, LINK VIGNETTE, and details.
#'
#' @details
#' Additional details...
#'     `sites_to_pool` is either a named list or a list of character vectors. If
#'     it is unnamed, than the name, min_a, max_a values are used from the first
#'     site listed. If the list is named, than the name, min_a, and max_a are
#'     used from the site in the name.
#'
#'
#'
#' @returns A `s4t_config` object, which contains
#'
#'
#' @examples
#' site_arrangement <- s4t_config(sites_names = c("A","B","C","D"),
#'                                min_a = c(1,1,1,1),
#'                                max_a = c(3,4,4,4),
#'                                sites_config = matrix(
#'                                        c(0,1,0,0,
#'                                          0,0,1,0,
#'                                          0,0,0,1,
#'                                          0,0,0,0),
#'                                          nrow = 4,
#'                                          ncol = 4,
#'                                          byrow = TRUE,
#'                                          dimnames = list(c("A","B","C","D"),
#'                                          c("A","B","C","D"))),
#'                                holdover_config = matrix(
#'                                         c(0,1,0,0,
#'                                           0,0,0,0,
#'                                           0,0,0,0,
#'                                           0,0,0,0),
#'                                           nrow = 4,
#'                                           ncol = 4,
#'                                           byrow = TRUE,
#'                                           dimnames = list(c("A","B","C","D"),
#'                                           c("A","B","C","D")))
#' )
#'
#'
#' @export
s4t_config <- function(sites_names,
                       sites_config,
                       holdover_config,
                       min_a,
                       max_a,
                       sites_to_pool = NULL) {


  max_a <- as.integer(max_a)
  min_a <- as.integer(min_a)

  if (sum(is.na(max_a))>1) {
    stop("max_a must be coercible to an integer")
  }

  if (sum(is.na(min_a))>1) {
    stop("min_a must be coercible to an integer")
  }

  if (length(max_a) != ncol(sites_config)) {
    stop(paste0("max_a should be length = ",ncol(sites_config)))
  }

  if (length(min_a) != ncol(sites_config)) {
    stop(paste0("minx_a should be length = ",ncol(sites_config)))
  }


  sites_names <- as.character(sites_names)
  if (sum(is.na(sites_names))>1) {
    stop("sites_names must be coercible to a character vector")
  }


  if (!is(sites_config,"matrix") | length(dim(sites_config)) != 2) {
    stop("sites_config must be a matrix")
  }


  if (!is(holdover_config,"matrix") | length(dim(holdover_config)) != 2) {
    stop("holdover_config must be a matrix")
  }


  if (!is.null(sites_to_pool)) {
    if (!is(sites_to_pool,"list")) {
      stop("sites_to_pool must be a list")
      # should also check if the names match
    }
    if (length(setdiff(names(sites_to_pool),sites_names)) > 0) {
      stop("Names of sites_to_pool not in sites_names")
    }
    for (i in 1:length(sites_to_pool)) {
      if (length(setdiff(sites_to_pool[[i]],sites_names)) > 0) {
        stop("Names of sites in sites_to_pool not in sites_names")
      }
    }
  }





  # check that sites_names and sites_config, holdover_config are the right dimensions
  if (length(unique(c(nrow(sites_config),ncol(sites_config),
                      nrow(holdover_config ),ncol(holdover_config ),
                      length(sites_names), length(max_a),
                      length(min_a)))) != 1) {
    stop("Dimensions of sites_config, holdover_config")
  }

  # check that if configs are named, that the names match sites_names
  if (!is.null(rownames(sites_config))) {
    if (any(!(rownames(sites_config) == sites_names))) {
      stop("If present, row and column names of sites_config must equal and in the same order as sites_names")
    }
  }
  if (!is.null(colnames(sites_config))) {
    if (any(!(colnames(sites_config) == sites_names))) {
      stop("If present, row and column names of sites_config must equal and in the same order as sites_names")
    }
  }
  if (!is.null(rownames(holdover_config))) {
    if (any(!(rownames(holdover_config) == sites_names))) {
      stop("If present, row and column names of holdover_config must equal and in the same order as sites_names")
    }
  }
  if (!is.null(colnames(holdover_config))) {
    if (any(!(colnames(holdover_config) == sites_names))) {
      stop("If present, row and column names of holdover_config must equal and in the same order as sites_names")
    }
  }

  # set min and max age info

  obs_min_a <- min_a
  obs_max_a <- max_a

  set_min_a <- obs_min_a  - min(obs_min_a) + 1
  set_max_a <- obs_max_a - min(obs_min_a) + 1


  # pool sites

  using_sites <- sites_names



  if (!is.null(sites_to_pool)) {
    # a
    if (is.null(names(sites_to_pool))) {
      names(sites_to_pool) <- 1:length(sites_to_pool)
      for (i in 1:length(sites_to_pool)) {
        names(sites_to_pool)[i] <- as.character(sites_to_pool[[i]][1])
      }
    }

    for (i in 1:length(sites_to_pool)) {

      tmp_sitename <- names(sites_to_pool)[i]

      if (tmp_sitename %in% sites_to_pool[[i]]) sites_to_pool[[i]] <- setdiff(sites_to_pool[[i]],tmp_sitename)


      ## ADD THE CHECK RIGHT HERE
      # tmp_sconfig <- sites_config[which(sites_names %in% sites_to_pool[[i]] | sites_names == tmp_sitename),]
      # tmp_seqal <- apply(tmp_sconfig,MARGIN = 1,FUN = function(x) length(unique(x)))

      tmp_hconfig <- holdover_config[which(sites_names %in% sites_to_pool[[i]] | sites_names == tmp_sitename),]
      tmp_heqal <- apply(tmp_hconfig,MARGIN = 1,FUN = function(x) length(unique(x)))

      if (any(tmp_heqal != 1)) {
        stop(paste0("Cannot pool sites with different holdovers: ",tmp_sitename," with: ",
                    paste0(sites_to_pool[[i]],collapse = ", ")))
      }

      # remove the sites that we are pooling
      using_sites <- setdiff(using_sites,sites_to_pool[[i]])

      using_sites <- ifelse(using_sites == tmp_sitename,tmp_sitename,using_sites)
      # sites_names <- using_sites

      # sites_config and holdovers config to drop:
      keepthesesites <- !(sites_names %in% sites_to_pool[[i]])
      sites_config <- sites_config[keepthesesites,keepthesesites]
      holdover_config <- holdover_config[keepthesesites,keepthesesites]

      sites_names <- using_sites

      # dropping min_a and max_a values
      obs_min_a <- obs_min_a[keepthesesites]
      obs_max_a <- obs_max_a[keepthesesites]
      set_min_a <- set_min_a[keepthesesites]
      set_max_a <- set_max_a[keepthesesites]

      colnames(holdover_config) <- rownames(holdover_config) <- sites_names
      colnames(sites_config) <- rownames(sites_config) <- sites_names

      # sites <- sites_config
      # holdover_config <- holdover_config


    }

  }


  obj <- list(sites_config = sites_config,
              holdover_config = holdover_config,
              sites_names = sites_names,
              sites_to_pool = sites_to_pool,

              max_a_overall = max_a_overall,

              obs_min_a = obs_min_a,
              obs_max_a = obs_max_a,

              set_min_a = set_min_a,
              set_max_a = set_max_a
              )

  class(obj) <- "s4t_config"

  return(obj)
}



#' Create site configuration for simple linear arrangement
#'
#' @description
#' Note that this can be used
#'
#' @param sites_names a character `vector` of the site names. Order indicates the
#'  direction of movement.
#' @param holdover_sites a character `vector` of the sites where after passing,
#'     individuals can holdover (take 1 or more time intervals) before reaching the next site.
#' @param min_a a `vector` of the minimum ageclass individuals in a site can be.
#'     Must be the same length and order as `sites_names`
#' @param max_a a `vector` of the maximum ageclass individuals in a site can be.
#'     Must be the same length and order as `sites_names`
#' @param sites_to_pool a named list of character vectors that contain the names of sites
#'      to be pooled together and treated as one site. See LINK VIGNETTE and details.
#'
#' @export
linear_s4t_config <- function(sites_names,
                                     holdover_sites = NULL,
                                     min_a,
                                     max_a,
                                     sites_to_pool = NULL) {
  sites_names <- as.character(sites_names)


  holdover_config <- sites_config <- matrix(0,nrow = length(sites_names),
                         ncol = length(sites_names),
                         dimnames = list(sites_names, sites_names))

  for (i in 1:(length(sites_names) - 1)) {
    sites_config[i,i + 1] <- 1
  }
  if (!is.null(holdover_sites)) {
   for (j in 1:length(holdover_sites)) {
     i <- which(sites_names == holdover_sites[j])
     holdover_config[i, i + 1] <- 1
   }
  }

  s4t_config(sites_names = sites_names,
             sites_config = sites_config,
             holdover_config = holdover_config,
             min_a = min_a,
             max_a = max_a,
             sites_to_pool = sites_to_pool)
}


#' Create site configuration for simple branching site arrangement
#'
#' @description
#' This is only for when there are two or more initial sites where individuals
#'     can be released, and then the next site is for all sites. After this
#'     shared site, there should be at least one additional site in order to estimate
#'     transitions and apparent survival (because detection probability cannot
#'     be estimated for the last site).
#'
#'
#' @param sites_names a character `vector` of the site names. Order indicates the
#'     direction of movement.
#' @param branch_sites a character `vector` of the initial branching sites
#' @param holdover_sites a character `vector` of the sites where after passing,
#'     individuals can holdover (take 1 or more time intervals) before reaching the next site.
#' @param min_a a `vector` of the minimum ageclass individuals in a site can be.
#'     Must be the same length and order as `sites_names`
#' @param max_a a `vector` of the maximum ageclass individuals in a site can be.
#'     Must be the same length and order as `sites_names`
#' @param sites_to_pool a named list of character vectors that contain the names of sites
#'      to be pooled together and treated as one site. See LINK VIGNETTE and details.
#'
#' @export
simplebranch_s4t_config <- function(sites_names,
                                    branch_sites,
                                    holdover_sites = NULL,
                                    min_a,
                                    max_a,
                                    sites_to_pool = NULL) {

  sites_names <- as.character(sites_names)
  branch_sites <- as.character(branch_sites)


  holdover_config <- sites_config <- matrix(0,nrow = length(sites_names),
                                            ncol = length(sites_names),
                                            dimnames = list(sites_names, sites_names))

  first_intersecting_site <- setdiff(sites_names, branch_sites)[1]

  position_first_intersecting_site <- which(sites_names == first_intersecting_site)

  for (j in 1:length(branch_sites)) {
    tmp_i <- which(sites_names == branch_sites[j])
    sites_config[tmp_i,position_first_intersecting_site] <- 1

    if (!is.null(holdover_sites)) {
      if (branch_sites[j] %in% holdover_sites) {
        holdover_config[tmp_i,position_first_intersecting_site] <- 1
      }

    }
  }



  for (i in position_first_intersecting_site:(length(sites_names) - 1)) {
    sites_config[i,i + 1] <- 1
    if (!is.null(holdover_sites)) {
      if (sites_names[i] %in% holdover_sites) {
        holdover_config[i,i + 1] <- 1
      }
    }
  }


  s4t_config(sites_names = sites_names,
             sites_config = sites_config,
             holdover_config = holdover_config,
             min_a = min_a,
             max_a = max_a,
             sites_to_pool = sites_to_pool)
}

