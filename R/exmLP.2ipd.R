#' @title Checks whether two IPD datasets can be matched with lpSolve::lp
#'
#' @param ipd1 a dataframe with n1 row and p column, where n1 is number of subjects of the first IPD, and p is the number of variables used in standardization.
#' @param ipd2 a dataframe with n2 row and p column, where n2 is number of subjects of the second IPD, and p is the number of variables used in standardization.
#' @param vars_to_match variables used for matching. if NULL, use all variables.
#' @param cat_vars_to_01 variable names for the categorical variables that need to be converted to indicator variables.
#' @param mean.constrained whether to restrict the weighted means to be within the ranges of observed means. Default is FALSE. When it is TRUE, there is a higher chance of not having a solution.
#'
#' @details If dummy variables are already created for the categorical variables in the data set, and are present in \code{ipd1} and \code{ipd2}, then \code{cat_vars_to_01} should be left as NULL.
#' 
#' @return \item{lp.check}{0 = OS can be conducted; 2 = OS cannot be conducted}
#' 
#' @export exmLP.2ipd
#'
#' @author Lillian Yau
#' 
## exmLP.2ipd(ipd1, ipd2) ## this would be the example, but ipd1 and ipd2 are not in the package yet

exmLP.2ipd <- function (ipd1, ipd2, vars_to_match = NULL, cat_vars_to_01 = NULL, mean.constrained = FALSE)
{
  ## check vars_to_match
  vars_to_match <- check_data(ipd1, 
                              ipd2, 
                              v.ars_to_match = vars_to_match,
                              c.at_vars_to_01 = cat_vars_to_01
  )
  ## extract only vars_to_match from both ipd's
  ipd1 <- data.frame(ipd1[vars_to_match])
  ipd2 <- data.frame(ipd2[vars_to_match])
  
  if(!is.null(cat_vars_to_01)){

    ## save original input data for later use
    ipd1.o <- ipd1 
    ipd2.o <- ipd2
    
    ## convert categorical variables to indicator variables
    ipd1 <- cat201_minus1(ipd1.o, v.cat = cat_vars_to_01)
    ipd2 <- cat201_minus1(ipd2.o, v.cat = cat_vars_to_01)
  }
  ##
  ## derivation for lpCheck starts here ::::
  ##
  ipd <- as.data.frame(rbind(-1 * ipd1, ipd2))
  oneszeros <- c(rep(1, nrow(ipd1)), rep(0, nrow(ipd2)))
  zerosones <- c(rep(0, nrow(ipd1)), rep(1, nrow(ipd2)))
  ipd <- as.data.frame(cbind(ipd, oneszeros, zerosones))
  p <- ncol(ipd1)
  ## f.con is the A matrix's left 3 colns in the appendix of the paper
  f.con <- as.matrix(t(ipd))  
  f.obj <- rep(0.5, ncol(f.con))
  f.rhs <- as.data.frame(t(c(rep(0, p), 1, 1)))
  f.dir <- rep("=", p + 2)
  
  if (mean.constrained == TRUE) { 
    
    ## re-define ipd1o and ipd2o keeping the reference level for all categorical variables
    if(!is.null(cat_vars_to_01)){
      ipd1 <- cat201(ipd1.o, v.cat = cat_vars_to_01)
      ipd2 <- cat201(ipd2.o, v.cat = cat_vars_to_01)
    }
    ##
    ipd1.bar <- colMeans(ipd1)
    ipd2.bar <- colMeans(ipd2)
    x <- as.data.frame(rbind(ipd1.bar, ipd2.bar))
    bar.min <- apply(x, 2, min)
    bar.max <- apply(x, 2, max)
    f.rhs <- cbind(f.rhs, 2 * t(bar.min), 2 * t(bar.max))
    f.dir <- c(f.dir, rep(">=", ncol(ipd1)), rep("<=", ncol(ipd1)))
    f.con <- data.frame(rbind(f.con, 
                              cbind(t(ipd1), t(ipd2)), 
                              cbind(t(ipd1), t(ipd2))))
  }
  lp.check <- lpSolve::lp(direction = "max", objective.in = f.obj, 
                          const.mat = f.con, const.dir = f.dir, const.rhs = f.rhs, 
                          transpose.constraints = TRUE)$status
  return(list(lp.check = lp.check))
}