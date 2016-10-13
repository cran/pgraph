#' Calculate the Projected Covariance of Two Vectors
#'
#' \code{projcov} calculate the projected distance covariance of two vectors given
#' common factors.
#'
#' @export
#' @importFrom stats cor
#' @importFrom stats predict
#' @importFrom stats cov
#' @importFrom energy dcov.test
#' @importFrom glmnet cv.glmnet
#' @importFrom splines ns
#' @importFrom SAM samQL
#' @importFrom stats coef
#' @importFrom stats lm.fit
#' @param x first vector
#' @param y second vector
#' @param b factor matrix
#' @param method projection method. Default = 'linear'.
#' @param one.SE whether to use the 1se rule for glmnet. Default = TRUE.
#' @param refit whether to refit the selected model. Default = TRUE.
#' @param randSeed the random seed for the program. Default = 0.
#' @return a list.
#' \item{test.pearson}{pearson correlection test statistic}
#' \item{test.dcov}{distance covariance test statistic}
#' \item{xeps}{residual of  projection of x on b}
#' \item{yeps}{residual of  projection of y on b}
#' @seealso \code{\link{greg}}, \code{\link{roc}}, \code{\link{pgraph}}
#' @examples
#' library(splines)
#' set.seed(0)
#' K = 3
#' n = 100
#' b = matrix(rnorm(K*n),n,K)
#' bx = 1:3
#' by = c(1,2,2)
#' x = b%*%bx+rnorm(n)
#' y = b%*%by+rnorm(n)
#' fit1 = projcov(x, y, b, method = "lasso")
#' fit2 = projcov(x, y, b, method = "sam")
projcov <- function(x, y, b, method = c("lasso","sam"), one.SE = TRUE, refit = TRUE, randSeed = 0){
  method = match.arg(method)
  set.seed(randSeed)
  xeps = projcore(x, b, method = method, one.SE = one.SE, refit = refit, randSeed = randSeed)
  yeps = projcore(y, b, method = method, one.SE = one.SE, refit = refit, randSeed = randSeed)

  test.pearson = abs(cor(xeps,yeps))
  test.dcov = energy::dcov.test(xeps,yeps)$statistic

return(list=list(test.pearson = test.pearson, test.dcov = test.dcov, xeps=xeps,yeps=yeps))
}
