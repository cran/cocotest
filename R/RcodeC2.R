#' @name coco2
#' @title Type 2 Rank correlation coefficient based test on positive dependence through stochastic ordering
#' @description This function evaluates the assumption of positive dependence through stochastic ordering in multiple comparison procedures
#' @param cordx a numeric vector
#' @param cordy a numeric vector
#' @param alpha a number of significance level
#' @param Rboot a number of bootstrap replicates
#' @param seed a number of the seed of random number generator
#' @return a vector of three numbers: a lower bound of one-sided confidence interval \code{lower_bound}, a test statistic \code{estimation}, and an indicator whether the PDS condition holds or not \code{PDS_assumption}
#' @export
#' @importFrom boot boot
#' @importFrom boot boot.ci
#' @importFrom stats cor
#' @author Jiangtao Gou
#' @author Fengqing Zhang
#' @details R package \code{boot} is included for computing nonparametric bootstrap confidence intervals
#' @references
#'  Gou, J., Wu, K. and Chen, O. Y. (2024). Rank correlation coefficient based tests on positive dependence through stochastic ordering with application in cancer studies, Technical Report.
#'  Gou, J. (2023). On dependence assumption in p-value based multiple test procedures. \emph{Journal of Biopharmaceutical Statistics}, 33(5), 596-610.
#'  Gou, J. (2024). A test of the dependence assumptions for the Simes-test-based multiple test procedures. \emph{Statistics in Biopharmaceutical Research}, 16(1), 1-7.
#' @examples
#' set.seed(123)
#' cordx <- rnorm(40)
#' cordy <- rnorm(40)
#' coco2(cordx, cordy)
#'
coco2 <- function(cordx, cordy, alpha = 0.05, Rboot = 100, seed = 1) {
  set.seed(seed = seed)
  #
  mytab <- matrix(c(cordx, cordy), ncol = 2)
  #
  sampleStat2 <- function(x, d) {
    sampletau <- stats::cor(x[d,1], x[d,2], method = "kendall")
    samplerho <- stats::cor(x[d,1], x[d,2], method = "spearman")
    samplestat2 <- (1 + 3*sampletau)/(1+samplerho)^2
    return(samplestat2)
  }
  #
  tryCatch({
    bootresult2 <- boot::boot(mytab, sampleStat2, R = Rboot)
    ciresult2 <- boot::boot.ci(bootresult2, conf = 1 - 2*alpha, type = "bca")
    SK2lb <- ciresult2$bca[4]
    SK2est <- ciresult2$t0
    if (SK2lb > 1) {
      SK2conclusion <- FALSE
    } else {
      SK2conclusion <- TRUE
    }
  }, error = function(e){}) # End of tryCatch
  return(list(lower_bound = SK2lb,
              estimation = SK2est,
              PDS_assumption = SK2conclusion))
}
