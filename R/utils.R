
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
