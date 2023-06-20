#' Checks whether two IPD datasets can be matched with lpSolve::lp
#'
#' @param ipd1 a dataframe with n1 row and p column, where n1 is number of subjects of the first IPD, and p is the number of variables used in standardization.
#' @param ipd2 a dataframe with n2 row and p column, where n2 is number of subjects of the second IPD, and p is the number of variables used in standardization.
#'
#' @return \item{lp.check}{0 = OS can be conducted; 2 = OS cannot be conducted}
#' @export osLP
#'
#' @author Lillian Yau
## osLP(ipd1, ipd2) ## this would be the example, but ipd1 and ipd2 are not in the package yet
osLP <- function (ipd1, ipd2)
{
  ipd <- as.data.frame(rbind(-1 * ipd1, ipd2))
  oneszeros <- c(rep(1, nrow(ipd1)), rep(0, nrow(ipd2)))
  zerosones <- c(rep(0, nrow(ipd1)), rep(1, nrow(ipd2)))
  ipd <- as.data.frame(cbind(ipd, oneszeros, zerosones))
  ##
  p <- ncol(ipd)
  f.con <- as.matrix(t(ipd))
  f.obj <- rep(0.5, ncol(f.con))
  f.rhs <- as.data.frame(t(c(rep(0, p - 2 ), 1, 1)))
  f.dir <- rep("=", p)
  lp.check <- lpSolve::lp(direction = "max",
                          objective.in = f.obj,
                          const.mat = f.con,
                          const.dir = f.dir,
                          const.rhs = f.rhs,
                          transpose.constraints = TRUE)$status
  return(list(lp.check = lp.check))
}
