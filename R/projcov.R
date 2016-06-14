#' Calculate the Projected covariance of two vectors
#'
#' \code{projcov} calculate the projected distance covariance of two vectors given
#' factors.
#'
#' @export
# @importFrom energy dcov.test
# @importFrom glmnet cv.glmnet
# @importFrom splines ns
# @importFrom SAM samQL
#' @importFrom splines ns
#' @param x first vector
#' @param y second vector
#' @param b factor matrix
#' @param method projection method. Default = 'linear'.
#' @return a list.
#' \item{test}{distance covariance test object}
#' \item{xeps}{residual of the projection of x}
#' \item{yeps}{residual of the projection of y}
#' @seealso \code{\link{greg}}, \code{\link{roc}}, \code{\link{pgraph}}, \code{\link{projcore}}
projcov <- function(x, y, b, method = c("lasso","sam")){
  method = match.arg(method)
  if(method == 'lasso'){
xfit = glmnet::cv.glmnet(b,x)
xeps = x - predict(xfit, b)
yfit = glmnet::cv.glmnet(b,y)
yeps = y - predict(yfit, b)
} else if(method == 'sam'){
  xfit = cv.samQL(b,x)
  xeps = x - predict(xfit$sam.fit, b)$values[,xfit$lambda.min]
  yfit = cv.samQL(b,y)
  yeps = y - predict(yfit$sam.fit, b)$values[,yfit$lambda.min]
}

  test.pearson = abs(cor(xeps,yeps))
#  test.dcov = dcor(xeps,yeps)
  test.dcov = energy::dcov.test(xeps,yeps)$statistic

return(list=list(test.pearson = test.pearson, test.dcov = test.dcov, xeps=xeps,yeps=yeps))
}
