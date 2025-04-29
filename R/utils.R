
flex_diag <- function(x) {
  if (is.null(dim(x))) return(x)

  diag(x)
}

misc_ind <- function(x) {
  y=vector();
  for (z in 1:x) {
    y <- c(y,z:x)
  }
  y
}

misc_ind_2 <- function(x,y) {
  a=vector();
  for (b in x:y) {
    a <- c(a,b:y)
  }
  a
}

which.min.not.zero <- function(x) {
  # x = as.numeric(x)
  x[x==0] = max(x) + 1
  which.min(x)
}

min_not_zero <- function(x, na.rm = FALSE) {
  x[x==0] = max(x,na.rm = na.rm) + 1
  min(x)
}

rep_crosses <- function(min_a,max_a) {
  rep(min_a:max_a,(max_a:min_a - min_a + 1))
}


# cut this function?
# printable_round <- function(x,digits = 0,...) {
#   y <- round(x=x,digits = digits, ...)
#
#   y_char <- as.character(y)
#
#   contain_period <- grepl("[.]",y_char)
#
#   if (all(contain_period) |
#       all(!contain_period)) {
#     y_char <- y_char
#   } else {
#     y_char[!contain_period] <- paste0(y_char[!contain_period],".")
#   }
#
#
#   # digs <- as.character(y - floor(y))
#   # lens <- sapply(digs, FUN = function(z) nchar(z))
#   #
#
#
#   # if (max(lens) == min(lens)) return(y_char)
#
#
#
#   digits_past_period <- gsub("^.*[.]","",y_char)
#
#   lens2 <- sapply(digits_past_period, FUN = function(z) nchar(z))
#
#   add_zero <- max(lens2) != lens2
#
#   # zeroes <- paste0()
#   for (i in which(add_zero == TRUE)) {
#     y_char[i] <- paste0(y_char[i],paste0(rep("0",times = (max(lens2) - lens2[i])),
#                                          collapse = ""),
#                         collapse = "")
#
#   }
#
#   return(y_char)
# }
