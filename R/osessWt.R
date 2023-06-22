#' @titile Optimal Standardization for two IPD's
#'
#' @param ipd1 a dataframe with n row and p column, where n is number of subjects and p is the number of variables used in matching.
#' @param ipd2 the other IPD with the same number of columns
#' @param catigorical.var a list of variable names for the categorical variables in the data
#'
#' @details If dummy variables are already created for the categorical variables in the data set, and are present in \code{ipd1} and \code{ipd2}$, then \code{categorical.var} should be left as NULL.
#' 
#' @return
#' \item{osess.wt }{weights of optimal standardization by maximizing ESS. Scaled to .......}
#' \item{ipd.ess }{effective sample size. It is no smaller than the ESS given by the MAIC weights.}
#' \item{ipd.wtsumm}{weighted summary statistics of the matching variables after matching. they should be identical to the input AD when AD is within the IPD convex hull.}
#' 
#' @export osessWt
#'
#' @examples
osessWt <- function (ipd1, ipd2, catigorical.var = NULL) 
{
  ipd <- as.data.frame(rbind(-1 * ipd1, ipd2))
  oneszeros <- c(rep(1, nrow(ipd1)), rep(0, nrow(ipd2)))
  zerosones <- c(rep(0, nrow(ipd1)), rep(1, nrow(ipd2)))
  ipd <- as.data.frame(cbind(ipd, oneszeros, zerosones))
  ##
  ## means of each variable for ipd1 and ipd2
  ##
  ipd1.bar <- colMeans(ipd1)
  ipd2.bar <- colMeans(ipd2)
  x <- as.data.frame(rbind(ipd1.bar, ipd2.bar))
  bar.min <- apply(x, 2, min)
  bar.max <- apply(x, 2, max) 
  ## rm(x)
  ########
  ipd.n <- nrow(ipd)
  # ones <- rep(1, ipd.n)
  # ipd0 <- data.frame(cbind(ipd, ones))
  # ad0 <- matrix(c(ad, 1, rep(0, nrow(ipd0))), nrow = 1)
  p <- ncol(ipd1) ## number of matching variables
  ad0 <- matrix(c(rep(0, p), 
                  1, 
                  1, 
                  rep(0, nrow(ipd))), 
                nrow = 1)
  # p <- length(ad)
  Amat <- as.matrix(ipd)
  Dmat <- diag(ipd.n)
  ##
  ## ipd blocks
  ##
  x0 <- as.data.frame(rbind(ipd1, ipd2))
  x1 <- as.data.frame(rbind(-1 * ipd1, -1 * ipd2))
  Amat0 <- as.matrix(data.frame(cbind(Amat, Dmat, x0, x1
  )))
  ##rm(x0); rm(x1)
  ## dvec is not needed for this purpose
  dvec <- rep(0, ipd.n)
  bvec <- matrix(c(ad0, 2 * bar.min, 2 * bar.max * (-1)
  ), nrow = 1)
  wts <- quadprog::solve.QP(Dmat = Dmat, 
                            dvec = dvec, 
                            Amat = Amat0, 
                            bvec = bvec, 
                            meq = p+2)
  ## wts.p <- quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat0, 
  ##                         bvec = bvec, meq = p)
  ##
  ipd.wts.me <- wts[["solution"]] ## * ipd.n
  ipd.ess.me <- round(sum(ipd.wts.me)^2/sum(ipd.wts.me^2), 1)
  ipd.wt.mean.me <- colMeans(ipd[, 1:p] * ipd.wts.me)
  ## weights etc. for ipd1
  ipd.1.wts <- ipd.wts.me[1 : nrow(ipd1)]
  ipd.1.wts <- ipd.1.wts * nrow(ipd1) ##/ ipd.n
  ipd.1.ess <- round(sum(ipd.1.wts)^2/sum(ipd.1.wts^2), 1)
  ipd.1.wt.mean <- colMeans(ipd1 * ipd.1.wts) ##* nrow(ipd1)
  ##ipd.1.wt.rs.mean <- colMeans(ipd1 * ipd.1.wts.rs)
  ## weights etc. for ipd2
  ipd.2.wts <- ipd.wts.me[(1 + nrow(ipd1)) : (nrow(ipd1) + nrow(ipd2))]
  ipd.2.wts <- ipd.2.wts * nrow(ipd2) ##/ ipd.n
  ipd.2.ess <- round(sum(ipd.2.wts)^2/sum(ipd.2.wts^2), 1)
  ipd.2.wt.mean <- colMeans(ipd2 * ipd.2.wts) ## * nrow(ipd2)
  ##ipd.2.wt.rs.mean <- colMeans(ipd2 * ipd.2.wts.rs)
  ##
  return(list(osess.wt = ipd.wts.me, 
              ipd.ess = ipd.ess.me, 
              ipd.wt.mean = ipd.wt.mean.me,
              ipd.1.wts = ipd.1.wts,
              ##ipd.1.wts.rs = ipd.1.wts.rs,
              ipd.1.ess = ipd.1.ess,
              ipd.1.wt.mean = ipd.1.wt.mean,
              ##ipd.1.wt.rs.mean = ipd.1.wt.rs.mean,
              ipd.2.wts = ipd.2.wts,
              ipd.2.ess = ipd.2.ess,
              ipd.2.wt.mean = ipd.2.wt.mean ##,
              ##ipd.2.wt.rs.mean = ipd.2.wt.rs.mean ##,
              ##ipd.2.wts.rs = ipd.2.wts.rs
  ))
}