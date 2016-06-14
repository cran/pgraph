#' Calculate the residual of a matrix projected on a factor matrix
#'
#' \code{projcore} calculates the residual of a matrix projected on a factor matrix using a method
#' specified.
#'
#' @export
#' @importFrom stats  cor
#' @importFrom stats predict
#' @importFrom stats cov
# @importFrom glmnet cv.glmnet
# @importFrom glmnet predict.glmnet
# @importFrom splines ns
# @importFrom SAM samQL
# @importFrom parcor ridge.cv
#' @param x matrix
#' @param b factor matrix
#' @param method projection method. Default = 'linear'.
#' @return a residual matrix
#' \item{resi}{a residual matrix}
#' @examples
#' library(splines)
#' K = 3
#' n = 100
#' b = matrix(rnorm(K*n),n,K)
#' bx = 1:3
#' by = c(1,2,2)
#' x = b%*%bx+rnorm(n)
#' y = b%*%by+rnorm(n)
#' fit1 = projcore(x, b, method = "lasso")
#' fit2 = projcore(y, b, method = "sam")
#' @seealso \code{\link{greg}}, \code{\link{roc}}, \code{\link{projcov}}, \code{\link{pgraph}}
projcore <- function(x, b, method = c("lasso","sam")){
  method = match.arg(method)
 # cor = match.arg(cor)
p = ncol(x)
resi = x
for (j in 1:p) {
  if(method == 'lasso'){
    xfit = glmnet::cv.glmnet(b,x[,j])
    resi[,j] = x[,j] - predict(xfit, b)
   } else if(method == 'sam'){
    xfit = cv.samQL(b,x[,j])
    resi[,j] = x[,j] - predict(xfit$sam.fit, b)$values[,xfit$lambda.min]
   }
}


return(resi = resi)
}
